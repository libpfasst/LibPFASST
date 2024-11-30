module pf_mod_optimization
  use pfasst
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_mpi
  use pf_mod_ndarray_oc
  use pf_mod_parallel_oc
  use feval
  use solutions
  implicit none

  real(pfdp), parameter, private :: armijoDecrease = 1e-4
  real(pfdp), parameter, private :: minStepSize    = 1e-6
  real(pfdp), parameter, private :: maxStepSize    = 1e6
  real(pfdp), parameter, private :: stepIncFactor  = 10
  real(pfdp), parameter, private :: c2             = 0.1 !0.9 
  integer,    parameter, private :: maxIterZoom    = 20
 
contains

  subroutine restrict_compress_interpolate(pf, t0, dt, coarse_level, flags)
    type(pf_pfasst_t), target, intent(inout) :: pf
    real(pfdp),                intent(in   ) :: t0, dt
    integer,                   intent(in   ) :: coarse_level, flags
    class(pf_level_t), pointer :: c_lev_p    
    class(pf_level_t), pointer :: f_lev_p
    real(pfdp), allocatable :: f_times(:), c_times(:), compressed_values(:)
    class(pf_encap_t), allocatable :: cf_delta(:)   !  Coarse in time but fine in space

    integer                 :: level_index, m, j, nvars, nnodes
    
    ! restrict to coarse level
    do level_index = pf%nlevels,coarse_level+1,-1
      f_lev_p => pf%levels(level_index);
      c_lev_p => pf%levels(level_index-1)
      allocate(f_times(f_lev_p%nnodes))
      f_times = t0 + dt*f_lev_p%nodes
      call restrict_ts(f_lev_p, c_lev_p, f_lev_p%Q, c_lev_p%Q, f_times, flags)
      deallocate(f_times)
    end do

!    if(compress) then
!      nvars = product(pf%levels(coarse_level)%shape)
!      nnodes = pf%levels(coarse_level)%nnodes
!    
!      allocate(compressed_values(nvars))
!    
!      do m=1, nnodes
!        call pf%levels(coarse_level)%Q(m)%pack(compressed_values,flags)
!        compressed_values(:) = floor(compressed_values(:)*N_digits)/dble(N_digits)  
!        call pf%levels(coarse_level)%Q(m)%unpack(compressed_values,flags)
!      end do
!    
!      deallocate(compressed_values)
!    end if
    
    ! and interpolate again   
    do level_index = coarse_level+1,pf%nlevels
      f_lev_p => pf%levels(level_index);
      c_lev_p => pf%levels(level_index-1)
      call f_lev_p%ulevel%factory%create_array(cf_delta, c_lev_p%nnodes, f_lev_p%index,  f_lev_p%lev_shape)
      allocate(c_times(c_lev_p%nnodes))
      c_times = t0 + dt*c_lev_p%nodes

      do m = 1, c_lev_p%nnodes
        call cf_delta(m)%setval(0.0_pfdp,flags)
        call f_lev_p%ulevel%interpolate(f_lev_p, c_lev_p, cf_delta(m), c_lev_p%Q(m), c_times(m), flags)
      end do
      call pf_apply_mat(f_lev_p%Q, 1.0_pfdp, f_lev_p%tmat, cf_delta, .true., flags)
      deallocate(c_times)
      call f_lev_p%ulevel%factory%destroy_array(cf_delta)
    end do
    
  end subroutine restrict_compress_interpolate

  subroutine evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objective, L2NormUSq, savedStates)
    type(pf_pfasst_t),         intent(inout) :: pf
    type(pf_ndarray_oc_t), target,  intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha
    integer,                   intent(in   ) :: nsteps
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(  out) :: objective, L2NormUSq
    real(pfdp),                intent(inout) :: savedStates(:,:,:)    !nsteps, nnodes, nvars
    integer             :: m, step, nnodes, thisstep, ierror, kk, level_index
    real(pfdp) :: stepobjective, roundedValue
    real(pfdp), pointer :: obj(:) !, savedStatesCompressed(:)
    
!     ! for restrict/interpolate (storage)
!     class(pf_level_t), pointer :: c_lev_p    
!     class(pf_level_t), pointer :: f_lev_p
!     real(pfdp), allocatable :: f_times(:)  !!  Simulation time at fine nodes
!     !---

    nnodes = pf%levels(pf%nlevels)%nnodes
    allocate(obj(nnodes))    
    objective = 0
    
!     allocate(savedStatesCompressed(size(savedStates,3)))
    
    pf%state%itcnt   = 0  ! reset iteration count here, not inside pf_pfasst_block_oc

    q1%yflatarray = 0.0_pfdp
    q1%pflatarray = 0.0_pfdp
    if (pf%rank .eq. 0) call initial(q1, 0.0_pfdp, dt) ! sequential run, this is always set at the beginning
    
    if ((pf%rank .eq. 0) .or. predict) then
      call pf%levels(pf%nlevels)%q0%copy(q1, 1) !do this only for rank 0 or predict = true?
    end if
!     call pf%levels(pf%nlevels)%q0%copy(q1, 1) 
    if (pf%rank .eq. 0) print *, ' *********  solve state **************'

    pf%q0_style = 0
    if(.not. predict) pf%q0_style = 2

    do step = 1, nsteps
    
      thisstep = (step-1)*pf%comm%nproc + pf%rank
      pf%state%pfblock = step

      
!     if (nsteps == 1 ) then
!       call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, predict, 1, step=thisstep)      
!     else
      if(.not. predict) then ! warm start, load saved state values
        do m = 1, nnodes
          call pf%levels(pf%nlevels)%Q(m)%unpack(savedStates(step,m,:), 1)
        end do
        if(compress) then
            call restrict_compress_interpolate(pf, thisstep*dt, dt, storage_level, 1)
!             do kk=1, size(savedStates,3)
!               roundedValue = floor(pf%levels(pf%nlevels)%Q(m)%yflatarray(kk)*N_digits)/N_digits
!               pf%levels(pf%nlevels)%Q(m)%yflatarray(kk) = roundedValue
!             end do           
!             savedStatesCompressed(:)= floor(savedStates(step,m,:)*N_digits)/dble(N_digits)           
!             call pf%levels(pf%nlevels)%Q(m)%unpack(savedStatesCompressed(:), 1)
! 
!             ! restrict to coarse level
!             do level_index = pf%nlevels,2,-1
!               f_lev_p => pf%levels(level_index);
!               c_lev_p => pf%levels(level_index-1)
!               allocate(f_times(f_lev_p%nnodes))
!               f_times = thisstep*dt + dt*f_lev_p%nodes
!               call restrict_ts(f_lev_p, c_lev_p, f_lev_p%Q, c_lev_p%Q, f_times, 1)
!               deallocate(f_times)
!             end do
!             ! and interpolate again   
!             do level_index = 2,pf%nlevels
!               call interpolate_time_space(pf, thisstep*dt, dt, level_index, .false., 1)
!             end do
            
          end if
!           else
!             call pf%levels(pf%nlevels)%Q(m)%unpack(savedStates(step,m,:), 1)
!           end if
!           call pf%levels(pf%nlevels)%Q(m)%unpack(savedStates(step,m,:), 1)

!         end do
!         ! use correct initial value at first rank, previous one at other ranks
        if (pf%rank .eq. 0) then
         call pf%levels(pf%nlevels)%Q(1)%copy(pf%levels(pf%nlevels)%q0, 1)
!          additional: get difference between q0 and Q(1) from savedStates and add that to the other Q(m)s? can only be done on rank 0
        else
!         if (pf%rank > 0) then
         call pf%levels(pf%nlevels)%q0%copy(pf%levels(pf%nlevels)%Q(1), 1)
        end if
        ! re-evaluate function values
        call pf%levels(pf%nlevels)%ulevel%sweeper%evaluate_all(pf, pf%state%finest_level, &
                                    thisstep*dt+dt*pf%levels(pf%nlevels)%nodes, flags=1, step=thisstep+1)
!         call restrict_and_evaluate(pf, thisstep*dt, dt, 1, thisstep+1)
      end if
    
      call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., 1, step=thisstep)      
!     end if
    
      ! evaluate objective on current step    
      do m = 1, nnodes
        call objective_function(pf%levels(pf%nlevels)%ulevel%sweeper, pf%levels(pf%nlevels)%Q(m), &
                                    product(pf%levels(pf%nlevels)%lev_shape), m, obj(m), step)

        ! record state solution
        call pf%levels(pf%nlevels)%Q(m)%pack(savedStates(step,m,:), 1)
      end do
      stepobjective = 0
      do m = 1, pf%levels(pf%nlevels)%nnodes-1
        stepobjective = stepobjective + &
                 (obj(m)+obj(m+1))*(pf%levels(pf%nlevels)%nodes(m+1)-pf%levels(pf%nlevels)%nodes(m))*dt
      end do
      objective = objective+stepobjective

      ! copy/broadcase qend to q0 for next step
      if( step < nsteps ) then        
        call pf%levels(pf%nlevels)%qend%pack(pf%levels(pf%nlevels)%send, 1)    !<  Pack away your last solution
        call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, pf%comm%nproc-1)
!         if (pf%rank .eq. 0) then
          call pf%levels(pf%nlevels)%q0%unpack(pf%levels(pf%nlevels)%send, 1)    !<  Everyone resets their q0
!         end if
!         if(pf%rank == pf%comm%nproc-1) then
!           call pf_mpi_send(pf, pf%levels(pf%nlevels), 111, .true., ierror, 1)
!         end if
!         if(pf%rank == 0) then
!           call pf_mpi_recv(pf, pf%levels(pf%nlevels), 111, .true., ierror, 1)
!           call pf%levels(pf%nlevels)%qend%unpack(pf%levels(pf%nlevels)%recv, 1)
!         end if
      end if
    end do
    
    objective = 0.5*objective !0.5 for trapezoidal rule
    call control_L2Q(pf%levels(pf%nlevels)%ulevel%sweeper, dt, pf%levels(pf%nlevels)%nodes, &
                                  product(pf%levels(pf%nlevels)%lev_shape), L2NormUSq)
    objective = 0.5*objective + 0.5*alpha*L2NormUSq
!     deallocate(savedStatesCompressed)
    deallocate(obj)
    
  end subroutine evaluate_objective


  subroutine evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedStates, savedAdjoint)
    type(pf_pfasst_t),        intent(inout) :: pf
    type(pf_ndarray_oc_t), target, intent(inout) :: q1
    real(pfdp),               intent(in   ) :: dt
    integer,                  intent(in   ) :: nsteps
    logical,                  intent(in   ) :: predict
    real(pfdp),               intent(  out) :: gradient(:,:,:), LinftyNormGrad, L2NormGradSq
    real(pfdp),               intent(inout) :: savedStates(:,:,:), savedAdjoint(:,:,:)

    integer :: m, step, nnodes, thisstep, ierror
    real(pfdp), pointer :: obj(:) !, savedAdjointCompressed(:)

!     allocate(savedAdjointCompressed(size(savedAdjoint,3)))
      
    nnodes = pf%levels(pf%nlevels)%nnodes
    q1%pflatarray = 0.0_pfdp
    if ((pf%rank==pf%comm%nproc-1) .or. predict) &
      call pf%levels(pf%nlevels)%qend%copy(q1, 2)  ! this is only done for final time
    

    pf%state%itcnt   = 0  ! reset iteration count here, not inside pf_pfasst_block_oc
    
    pf%q0_style = 1
    if(.not. predict) pf%q0_style = 2

       
    
    if (pf%rank .eq. 0) print *, '*********  solve adjoint *************'
    
     do step = nsteps, 1, -1

        thisstep = (step-1)*pf%comm%nproc + pf%rank
        pf%state%pfblock = step

!       if (nsteps == 1) then
!         call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, predict, 2, step=thisstep) 
!       else
        ! assign savedStates to [Q(m), 1]
        do m = 1, nnodes
          call pf%levels(pf%nlevels)%Q(m)%unpack(savedStates(step,m,:), 1)
        end do
        if(inexact_adjoint) then
           call restrict_compress_interpolate(pf, thisstep*dt, dt, storage_level, 1)
        end if
        call restrict_for_adjoint(pf, thisstep*dt, dt, 1, thisstep+1) ! save states at all levels instead of restricting?      
        
        if(.not. predict) then
          do m = 1, nnodes
            call pf%levels(pf%nlevels)%Q(m)%unpack(savedAdjoint(step,m,:), 2)
           end do
            if(compress) then
              call restrict_compress_interpolate(pf, thisstep*dt, dt, storage_level, 2)
            end if
!             call pf%levels(pf%nlevels)%Q(m)%unpack(savedAdjoint(step,m,:), 2)

! !             do kk=1, size(savedStates,3)
! !               roundedValue = floor(pf%levels(pf%nlevels)%Q(m)%yflatarray(kk)*N_digits)/N_digits
! !               pf%levels(pf%nlevels)%Q(m)%yflatarray(kk) = roundedValue
! !             end do
!             savedAdjointCompressed(:)= floor(savedAdjoint(step,m,:)*N_digits)/dble(N_digits)
!             call pf%levels(pf%nlevels)%Q(m)%unpack(savedAdjointCompressed(:), 2)
! !           end if
!           else
!             call pf%levels(pf%nlevels)%Q(m)%unpack(savedAdjoint(step,m,:), 2)
!             end if
!           end do
          if(pf%rank==pf%comm%nproc-1) then 
              ! this happens in sweep anyway, but not in restriction of q0/qend at start of predictor
            call pf%levels(pf%nlevels)%Q(nnodes)%copy(pf%levels(pf%nlevels)%qend, 2)
          else
!           if(pf%rank < pf%comm%nproc-1) then 
            call pf%levels(pf%nlevels)%qend%copy(pf%levels(pf%nlevels)%Q(nnodes), 2)
          end if
          ! re-evaluate function values
          call pf%levels(pf%nlevels)%ulevel%sweeper%evaluate_all(pf, pf%state%finest_level, &
                                    thisstep*dt+dt*pf%levels(pf%nlevels)%nodes, flags=2, step=thisstep+1)
!           call restrict_and_evaluate(pf, thisstep*dt, dt, 2, thisstep+1)
        end if
       
        call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., 2, step=thisstep) 
!       end if
        
        do m = 1, pf%levels(pf%nlevels)%nnodes
          call pf%levels(pf%nlevels)%Q(m)%pack(gradient(step, m, :), 2)
          savedAdjoint(step,m,:) = gradient(step,m,:)
        end do
                                    
!         call pf%levels(pf%nlevels)%qend%copy(pf%levels(pf%nlevels)%q0, 2)
        if(step > 1) then
          call pf%levels(pf%nlevels)%q0%pack(pf%levels(pf%nlevels)%send, 2)    !<  Pack away your last solution
          call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, 0)
!           if(pf%rank == pf%comm%nproc-1) then
            call pf%levels(pf%nlevels)%qend%unpack(pf%levels(pf%nlevels)%send, 2)    !<  Everyone resets their qend
!           end if
!           if(pf%rank == 0) then
!             call pf_mpi_send(pf, pf%levels(pf%nlevels), 111, .true., ierror, 2)
!           end if
!           if(pf%rank == pf%comm%nproc-1) then
!             call pf_mpi_recv(pf, pf%levels(pf%nlevels), 111, .true., ierror, 2)
!             call pf%levels(pf%nlevels)%q0%unpack(pf%levels(pf%nlevels)%recv, 2)
!           end if
        end if
     end do

    call construct_gradient(pf%levels(pf%nlevels)%ulevel%sweeper, gradient, pf%levels(pf%nlevels)%nodes, &
                                      LinftyNormGrad, L2NormGradSq)

    pf%q0_style = 0     
!     deallocate(savedAdjointCompressed)
  end subroutine evaluate_gradient

!  subroutine armijo_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
!                          globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
!     type(pf_pfasst_t),         intent(inout) :: pf
!     type(ndarray_oc), target,  intent(inout) :: q1
!     real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
!     integer,                   intent(in   ) :: nsteps
!     integer,                   intent(inout) :: itersState, itersAdjoint
!     logical,                   intent(in   ) :: predict
!     real(pfdp),                intent(in   ) :: searchDir(:,:)  
!     real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq
!     real(pfdp),                intent(inout) :: stepSize, LinftyNormGrad, L2NormGradSq, gradient(:,:), savedAdjoint(:,:)
!     logical,                   intent(inout) :: stepTooSmall
!        
!     real(pfdp) :: globObjNew, directionTimesGradient, globDirXGradNew
!     integer    :: l, ierror
!     !real(pfdp) :: armijoDecrease ! this should be in a structure global to the module, and set in
!     !                             ! something like init_optimization, along with other parameters 
!     !real(pfdp) :: minStepSize
!     !armijoDecrease = 1e-4
!     !minStepSize    = 1e-6
! 
!     do
!       call update_control(pf%levels(pf%nlevels)%ctx, searchDir, stepSize)
!       !restrict control
!       do l = pf%nlevels-1,1,-1
! 	call restrict_control(pf%levels(l)%ctx, pf%levels(l+1)%ctx)
!       end do
! 
!       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint)
!       itersState = itersState + pf%state%itcnt
!            
!       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!       if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew
! 
!       if (globObjNew < globObj + armijoDecrease*stepSize*globDirXGrad) then
!         ! evaluate gradient to be consistent with wolfe_powell_step
!         call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
!         itersAdjoint = itersAdjoint + pf%state%itcnt
!         return 
!       end if
! 
!        call update_control(pf%levels(pf%nlevels)%ctx, searchDir, -stepSize)
!        ! no need to restrict here, coarse levels get overwritten above
!        stepSize = 0.5 * stepSize
!        if (stepSize < minStepSize) then
!          stepTooSmall = .true.
!          if (pf%rank .eq. 0) print *, 'no step found, stopping'
!          !call write_control(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
!          return
!        end if
!      end do
!    end subroutine armijo_step
! 
! 
! 
!    subroutine zoom(stepLowIn, stepHighIn, stepSizeIn, stepSizeOut, &
!                    pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedStates, alpha, &
!                    globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepTooSmall, globObjLowIn)
!      type(pf_pfasst_t),         intent(inout) :: pf
!      type(ndarray_oc), target,  intent(inout) :: q1
!      real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
!      integer,                   intent(in   ) :: nsteps
!      integer,                   intent(inout) :: itersState, itersAdjoint
!      logical,                   intent(in   ) :: predict
!      real(pfdp),                intent(in   ) :: searchDir(:,:,:) 
!      real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSizeOut
!      real(pfdp),                intent(in   ) :: stepSizeIn, stepLowIn, stepHighIn, globObjLowIn
!      real(pfdp),                intent(inout) :: gradient(:,:,:), savedStates(:,:,:)
!      logical,                   intent(inout) :: stepTooSmall
! 
!      real(pfdp) :: directionTimesGradient, directionTimesGradientNew
!      real(pfdp) :: stepSize, stepLow, stepHigh, objectiveLow
!      integer    :: l, ierror, iter
!  
!      stepSize = stepSizeIn
!      stepLow = stepLowIn
!      stepHigh = stepHighIn
!      objectiveLow = globObjLowIn
!  
!      iter = 0
!      if (pf%rank .eq. 0) print *, 'in zoom'
! 
!      do
!        call update_control(pf%levels(pf%nlevels)%ctx, searchDir, -stepSize)
!        stepSize = 0.5_pfdp*(stepLow+stepHigh)
!        iter = iter + 1
!        if (stepSize < minStepSize .or. iter > maxIterZoom) then
!          stepTooSmall = .true.
!          if (pf%rank .eq. 0) print *, 'no step found (< min or > maxIter), stopping'
!          !call write_control(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
!          stepSizeOut = stepSize
!          return
!        end if
! 
!        call update_control(pf%levels(pf%nlevels)%ctx, searchDir, stepSize)
!         !restrict control
!         do l = pf%nlevels-1,1,-1
!    	  call restrict_control(pf%levels(l)%ctx, pf%levels(l+1)%ctx)
!         end do
! 
!         call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStates)
!         itersState = itersState + pf%state%itcnt
!         if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', objectiveNew
!       
!         if ( (objectiveNew > globObj + armijoDecrease*stepSize*globDirXGrad) .or. &
!              (objectiveNew >= objectiveLow) )                                     then
!           stepHigh = stepSize
!           if (pf%rank .eq. 0) print *, 'set new stepHigh to stepSize; stepHigh = ', stepHigh, 'stepLow = ', stepLow
!         else          
!           call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedStates)
!           itersAdjoint = itersAdjoint + pf%state%itcnt
!           directionTimesGradientNew = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
!           
!           if (abs(directionTimesGradientNew) <= -c2*globDirXGrad) then
!             stepSizeOut = stepSize
!             return
!           end if
!    
!           if (directionTimesGradientNew*(stepHigh-stepLow) >= 0) then
!             if(pf%rank == 0) print *, stepSize, stepHigh, stepLow
!             stepHigh = stepLow
!             if (pf%rank .eq. 0) print *, 'set new stepHigh to stepLow; stepHigh = ', stepHigh, 'stepLow = ', stepLow, 'stepSize = ', stepSize 
!           end if
!           stepLow = stepSize
!           if (pf%rank .eq. 0) print *, 'set new stepLow; stepHigh = ', stepHigh, 'stepLow = ', stepLow, 'stepSize = ', stepSize 
!           objectiveLow = objectiveNew
!         end if
!       end do
!    end subroutine zoom
!     
!  
!    subroutine wolfe_powell_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedStates, alpha, &
!                                 objective, directionTimesGradient, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
!     type(pf_pfasst_t),         intent(inout) :: pf
!     type(ndarray_oc), target,  intent(inout) :: q1
!     real(pfdp),                intent(in   ) :: dt, alpha, objective, directionTimesGradient
!     integer,                   intent(in   ) :: nsteps
!     integer,                   intent(inout) :: itersState, itersAdjoint
!     logical,                   intent(in   ) :: predict
!     real(pfdp),                intent(in   ) :: searchDir(:,:,:) 
!     real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq
!     real(pfdp),                intent(inout) :: stepSize, gradient(:,:,:), savedStates(:,:,:)
!     logical,                   intent(inout) :: stepTooSmall
! 
!     integer    :: l, ierror, iter
!     real(pfdp) :: prevStepSize, directionTimesGradientNew, prevObjective
!     prevStepSize   = 0.0_pfdp
!     prevObjective = objective
!     iter = 1
! !     
!     do
!       call update_control(pf%levels(pf%nlevels)%ctx, searchDir, stepSize)
!       !restrict control
!       do l = pf%nlevels-1,1,-1
! 	call restrict_control(pf%levels(l)%ctx, pf%levels(l+1)%ctx)
!       end do
!  
!       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStates)
!       itersState = itersState + pf%state%itcnt
!       if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', objectiveNew
!        
!       if ( (objectiveNew > objective + armijoDecrease*stepSize*directionTimesGradient) .or. &
!            ( objectiveNew >= prevObjective .and. iter > 1) )                   then
!         call zoom(prevStepSize, stepSize, stepSize, stepSize, &
!                   pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedStates, alpha, objective, directionTimesGradient, &
!                   objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepTooSmall, prevObjective)
!         return
!       end if                    
! !       
!       if (pf%rank .eq. 0) print *, 'evaluate new gradient'
!       call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedStates)
!       itersAdjoint = itersAdjoint + pf%state%itcnt
!       directionTimesGradientNew = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
!       if (pf%rank .eq. 0) print *, 'globDirXGradNew = ', directionTimesGradientNew
!       if (abs(directionTimesGradientNew) <= -c2*directionTimesGradient) then
!         return
!       end if
! ! 
!       if (directionTimesGradientNew >= 0) then
!         call zoom(stepSize, prevStepSize, stepSize, stepSize, &
!                   pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedStates, alpha, objective, directionTimesGradient, &
!                   objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepTooSmall, objectiveNew)
!         return
!       end if
! 
!       prevObjective  = objectiveNew
!       prevStepSize = stepSize
!       stepSize     = stepIncFactor*stepSize
!       iter = iter + 1
!       if (stepSize > maxStepSize) then
!          stepTooSmall = .true.
!          if (pf%rank .eq. 0) print *, 'no step found (> max), stopping'         
!          return
!        end if
!        call update_control(pf%levels(pf%nlevels)%ctx, searchDir, -prevStepSize)
!     end do
! 
!    end subroutine wolfe_powell_step
   
    subroutine strong_wolfe_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, &
        savedStates, savedAdjoint, &
        alpha, globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, &
        stepSize, stepTooSmall)
     type(pf_pfasst_t),         intent(inout) :: pf
     type(pf_ndarray_oc_t), target,  intent(inout) :: q1
     real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
     integer,                   intent(in   ) :: nsteps
     integer,                   intent(inout) :: itersState, itersAdjoint
     logical,                   intent(in   ) :: predict
     real(pfdp),                intent(in   ) :: searchDir(:,:,:) 
     real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq
     real(pfdp),                intent(inout) :: stepSize, gradient(:,:,:), savedStates(:,:,:), savedAdjoint(:,:,:)
     logical,                   intent(inout) :: stepTooSmall
     integer    :: l, ierror
     real(pfdp) :: low, high
     real(pfdp) :: directionTimesGradient, globObjNew, globDirXGradNew
     logical    :: first
     
     first = .true.
     
     low  = 0.0_pfdp
     high = stepSize
     !high = 1.0_pfdp
     !if (stepSize > high) high = stepSize
     
     do
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, high) !stepSize)
       !restrict control
       do l = pf%nlevels-1,1,-1
       call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
       end do
 !       call dump_control(pf%levels(pf%nlevels)%ctx, pf, 'u1')
 
       !call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStates)
       first = .false.
       
       itersState = itersState + pf%state%itcnt
       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
       if(pf%rank == 0) print *, high, 'objectiveNew (L2) = ', globObjNew ! *, stepSize,
 
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -high) !-stepSize)
       
       if (globObjNew < globObj + armijoDecrease*stepSize*globDirXGrad) then
         high = 2.0*high
       else ! Armijo condition not satisfied, exit loop and reduce stepsize
         exit
       end if
     end do
       
     do
       stepSize = 0.5*(low+high)
       
       if (stepSize < minStepSize ) then
         stepTooSmall = .true.
         if (pf%rank .eq. 0) print *, 'no step found (< min), stopping'
         return
       end if
       
       if (abs(low-high) < minStepSize ) then
         stepTooSmall = .true.
         if (pf%rank .eq. 0) print *, 'no step found (low==high), stopping'
         return
       end if
       
       
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, stepSize)
       do l = pf%nlevels-1,1,-1
        call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
       end do
 
       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStates)
       !(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
       itersState = itersState + pf%state%itcnt
       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
       if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew
       
       if (globObjNew >= globObj + armijoDecrease*stepSize*globDirXGrad) then
         call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
         high = stepSize
         cycle
       end if
 
       ! now Armijo is satisfied again
       !call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
       call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedStates, savedAdjoint)
       itersAdjoint = itersAdjoint + pf%state%itcnt
       directionTimesGradient = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
       call mpi_allreduce(directionTimesGradient, globDirXGradNew, 1, MPI_REAL8, MPI_SUM, &
                          pf%comm%comm, ierror)          
                          
       if (pf%rank .eq. 0) print *, 'globDirXGradNew = ', globDirXGradNew, '-c2*globDirXGrad', -c2*globDirXGrad
 
           
 !       if (abs(globDirXGradNew) <= -c2*globDirXGrad) then  ! strong Wolfe conditions satisfied
 !         return
 !       end if
    
       if (globDirXGradNew < c2*globDirXGrad) then
 !         if(pf%rank == 0) print *, stepSize, high, low
         low = stepSize
         if (pf%rank .eq. 0) print *, 'set new stepLow;  stepHigh = ', high, 'stepLow = ', low, 'stepSize = ', stepSize 
       elseif (globDirXGradNew > -c2*globDirXGrad) then
         high = stepSize
         if (pf%rank .eq. 0) print *, 'set new stepHigh; stepHigh = ', high, 'stepLow = ', low, 'stepSize = ', stepSize 
       else
         exit
       end if
       
       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
     end do
     
   end subroutine strong_wolfe_step
   
end module pf_mod_optimization

