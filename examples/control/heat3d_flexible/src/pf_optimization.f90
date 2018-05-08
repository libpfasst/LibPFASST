module pf_mod_optimization
  use pfasst
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_mpi
  use pf_mod_ndarray_oc
  use feval
  use probin, only: solve_y
  use solutions
  implicit none

  real(pfdp), parameter, private :: armijoDecrease = 1e-4
  real(pfdp), parameter, private :: minStepSize    = 1e-6
  real(pfdp), parameter, private :: maxStepSize    = 1e6
  real(pfdp), parameter, private :: stepIncFactor  = 10
  real(pfdp), parameter, private :: c2             = 0.9 
  integer,    parameter, private :: maxIterZoom    = 20
  
 
contains

  subroutine evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objective, L2NormUSq, savedAdjoint, predictAdj)
    type(pf_pfasst_t),         intent(inout) :: pf
    type(ndarray_oc),          intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha
    integer,                   intent(in   ) :: nsteps
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(  out) :: objective, L2NormUSq
    real(pfdp),                intent(inout) :: savedAdjoint(:,:,:)
    logical, optional,         intent(in   ) :: predictAdj
    integer             :: m, pred_flags(1), l
    real(pfdp), pointer :: obj(:), tmp(:)
    real(pfdp) :: t(pf%levels(pf%nlevels)%nnodes)
    logical    :: predictAdjoint
    
    predictAdjoint = .true.
    if (present(predictAdj) ) predictAdjoint = predictAdj
    allocate(tmp(product(pf%levels(pf%nlevels)%shape)))

   
    solve_y = .true.

    if (pf%rank .eq. 0) then
      call initial(q1, 0.0_pfdp, dt) !dt should be given by rank !initial_rd
    else
      q1%yflatarray = 0.0_pfdp
    end if
    if ((pf%rank .eq. 0) .or. predict) then
      call pf%levels(pf%nlevels)%q0%copy(q1, 1) !do this only for rank 0 or predict = true?
    end if
    
    if (do_mixed .eq. 1 ) then
!       pred_flags(1) = 2
      q1%pflatarray = 0.0_pfdp
      call pf%levels(pf%nlevels)%qend%copy(q1, 2)
!       if( predict .eqv. .false. .and. predictAdjoint .eqv. .true. ) then
        ! call predictor just for adjoint part
!         call pf_predictor(pf, pf%rank*dt, dt, pred_flags)     
!       end if
    end if

    if (pf%rank .eq. 0) print *, ' *********  solve state **************'
    
    if (do_mixed .eq. 1) then
       call pf_pfasst_block_oc(pf, dt, nsteps, predict, flags=0)
!        if (pf%rank .eq. 0) print *, ' *********  now adjoint w/o comm **************'
!        do l = 1, pf%nlevels
!           q1%pflatarray = 0.0_pfdp
!           q1%yflatarray = 0.0_pfdp
!           call pf%levels(l)%q0%copy(q1, 2)
!           call pf%levels(l)%ulevel%sweeper%spreadq0(pf%levels(l), (pf%rank+1)*dt, 2, pf%state%step+1)
!        end do
!        call pf_pfasst_block_oc(pf, dt, nsteps, .false., flags=0)
       
!        ! need to save values in adjoint component of Q(m)!
! !        if(pf%rank .eq. 0) print *, 'saving adjoint'
!        do m = 1, pf%levels(pf%nlevels)%nnodes
!          call pf%levels(pf%nlevels)%Q(m)%pack(tmp, 2)
!          savedAdjoint(m,:,:) = reshape(tmp, (/pf%levels(pf%nlevels)%shape(1), pf%levels(pf%nlevels)%shape(2)/))
!        end do
!        if(pf%rank == 0) print *, "max adjoint val", maxval(savedAdjoint(1,:,:))
    else
       call pf_pfasst_block_oc(pf, dt, nsteps, predict, 1)
    end if
    
!     if (pf%rank .eq. 0) print *, ' *********  compute objective **************'

    allocate(obj(pf%levels(pf%nlevels)%nnodes))    
    do m = 1, pf%levels(pf%nlevels)%nnodes
!       if (pf%rank .eq. 0) print *, m
      call objective_function(pf%levels(pf%nlevels)%ulevel%sweeper, pf%levels(pf%nlevels)%Q(m), &
                                    pf%levels(pf%nlevels)%shape, m, obj(m))
    end do
    objective = 0.0
    do m = 1, pf%levels(pf%nlevels)%nnodes-1
      objective = objective + &
                 (obj(m)+obj(m+1))*(pf%levels(pf%nlevels)%nodes(m+1)-pf%levels(pf%nlevels)%nodes(m))*dt
    end do
    objective = 0.5*objective !0.5 for trapezoidal rule
!     print *, "0.5||y-y_d||^2", 0.5*objective
    call control_L2Q(pf%levels(pf%nlevels)%ulevel%sweeper, dt, pf%levels(pf%nlevels)%nodes, &
                                  pf%levels(pf%nlevels)%shape, L2NormUSq)
    objective = 0.5*objective + 0.5*alpha*L2NormUSq
    deallocate(obj)
    deallocate(tmp)
  end subroutine evaluate_objective



  ! assumes state solution in Q(:), component 1, of pf object
  ! how to do that dimension independent? this is just for 1d because of the gradient array
  subroutine evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
    type(pf_pfasst_t),        intent(inout) :: pf
    type(ndarray_oc), target, intent(inout) :: q1
    real(pfdp),               intent(in   ) :: dt
    integer,                  intent(in   ) :: nsteps
    logical,                  intent(in   ) :: predict
    real(pfdp),               intent(  out) :: gradient(:,:,:), LinftyNormGrad, L2NormGradSq
    real(pfdp),               intent(inout) :: savedAdjoint(:,:,:)
    integer :: m, n, nnodes
    integer :: dest, source, ierror, stat(MPI_STATUS_SIZE)
    real(pfdp), pointer :: tmp(:), terminal(:), adjointAtNode(:,:)
    real(pfdp) :: nodes(pf%levels(pf%nlevels)%nnodes), tend, t
    type(ndarray_oc) :: q_tmp
    
    allocate(tmp(product(pf%levels(pf%nlevels)%shape)))
    allocate(terminal(product(pf%levels(pf%nlevels)%shape)))
    terminal = 0.0_pfdp
    
    ! test: fill with analytic p_tilde on interval, works...
!     if (do_mixed == 1) then
!       nnodes = pf%levels(pf%nlevels)%nnodes
!       nodes = pf%levels(pf%nlevels)%nodes
!       do m=1, nnodes
!         t = pf%rank*dt + dt*(nodes(m)-nodes(1))
!         tend = (pf%rank+1)*dt
!         adjointAtNode => get_array2d_oc(pf%levels(pf%nlevels)%Q(m), 2)
!         call p_tilde(adjointAtNode, shape(adjointAtNode), t, tend)
!       end do
!     end if
    
    
    if (do_mixed .eq. 1 ) then       ! receive terminal value
       if(pf%rank /= pf%comm%nproc-1) then
          source = modulo(pf%rank+1, pf%comm%nproc)
          print *, pf%rank, "receiving from", source, "with buflen", pf%levels(pf%nlevels)%mpibuflen
          call mpi_recv(terminal, pf%levels(pf%nlevels)%mpibuflen, MPI_REAL8, &
                        source, 1, pf%comm%comm, stat, ierror)
          print *, pf%rank, "received ", maxval(terminal)
       end if
       
       ! apply matrix exponential to compute value at Q(1;2)
       call ndarray_oc_build(q_tmp, pf%levels(pf%nlevels)%shape) ! todo: destroy
       call q_tmp%unpack(terminal, 2)
      m = pf%levels(pf%nlevels)%nnodes          
      print *, pf%rank, m, "before exp", q_tmp%norm(2), pf%levels(pf%nlevels)%Q(m)%norm(2)
      call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%copy(q_tmp, 2)
      call pf%levels(pf%nlevels)%qend%copy(pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes), 2)

      print *, pf%rank, m, "after exp", q_tmp%norm(2), pf%levels(pf%nlevels)%Q(m)%norm(2)
       
!        call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%unpack(terminal, 2)
!        call pf%levels(pf%nlevels)%Q(1)%copy(pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes), 2)
       
       ! compute fft of data, multiply with, compute ifft
!        if(pf%rank == 0) print *, "compute Q(1), dt = ", dt
       m = 1          
       print *, pf%rank, m, "before exp", q_tmp%norm(2), pf%levels(pf%nlevels)%Q(m)%norm(2)

       call compute_exp_lap_dt_times_vec(pf%levels(pf%nlevels)%ulevel%sweeper, dt, q_tmp)
       call pf%levels(pf%nlevels)%Q(1)%axpy(1.0_pfdp, q_tmp, 2)
       call pf%levels(pf%nlevels)%q0%copy(pf%levels(pf%nlevels)%Q(1), 2)
       print *, pf%rank, m, "after exp", q_tmp%norm(2), pf%levels(pf%nlevels)%Q(m)%norm(2)

       
       ! send new Q(1;2)
       if (pf%rank /= 0) then
          call pf%levels(pf%nlevels)%Q(1)%pack(tmp, 2)
          dest = modulo(pf%rank-1, pf%comm%nproc)
          
          call mpi_send(tmp, pf%levels(pf%nlevels)%mpibuflen, MPI_REAL8, &   !mpi_send(savedAdjoint(1,:),...
                        dest, 1, pf%comm%comm, stat, ierror)
       end if
       
       tend = (pf%rank+1)*dt
!        print *, pf%rank, "tend", tend
       ! compute remaining values
       if(pf%rank /= pf%comm%nproc-1) then
          nnodes = pf%levels(pf%nlevels)%nnodes
          nodes = pf%levels(pf%nlevels)%nodes
          do m = nnodes-1, 2, -1
!           if(pf%rank == 0) print *, "compute Q(m)", m, " dt = ", dt*(nodes(nnodes)-nodes(m))
             call q_tmp%unpack(terminal, 2)
              t = pf%rank*dt + dt*(nodes(m)-nodes(1))
          print *, pf%rank, m, "before exp", q_tmp%norm(2), pf%levels(pf%nlevels)%Q(m)%norm(2)
             call compute_exp_lap_dt_times_vec(pf%levels(pf%nlevels)%ulevel%sweeper, &
                                            tend-t, q_tmp)
             call pf%levels(pf%nlevels)%Q(m)%axpy(1.0_pfdp, q_tmp, 2)
          print *, pf%rank, m, "after exp", q_tmp%norm(2), pf%levels(pf%nlevels)%Q(m)%norm(2)
          end do
       end if
       return
!        ?restrict to the coarser levels?
       
       ! to write the npy output
!        call call_hooks(pf, -1, PF_POST_STEP)
!         call ndarray_oc_dump_all_hook(pf, pf%levels(pf%nlevels), pf%state)
    else   
      q1%pflatarray = 0.0_pfdp ! only if no final time term in objective, otherwise this is nonzero, and terminal needs to be added to this
    
      call pf%levels(pf%nlevels)%qend%copy(q1, 2)
      ! do this only for final step or predict = true?
      !call restrict_for_adjoint(pf, 1)

    
      if (pf%rank .eq. 0) print *, '*********  solve adjoint *************'
      if(predict) pf%q0_style = 1
      call pf_pfasst_block_oc(pf, dt, nsteps, predict, 2) !predict
      pf%q0_style = 0
    end if ! do_mixed
    
    if (pf%rank .eq. 0) print *, '*********  compute gradient *************'
    
    do m = 1, pf%levels(pf%nlevels)%nnodes
      call pf%levels(pf%nlevels)%Q(m)%pack(tmp, 2)
      gradient(m,:,:) = reshape(tmp, (/pf%levels(pf%nlevels)%shape(1), pf%levels(pf%nlevels)%shape(2)/))
    end do

    call construct_gradient(pf%levels(pf%nlevels)%ulevel%sweeper, gradient, pf%levels(pf%nlevels)%nodes, &
                                      LinftyNormGrad, L2NormGradSq)
    deallocate(terminal)
    deallocate(tmp)
  end subroutine evaluate_gradient



  subroutine armijo_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
                         globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
    type(pf_pfasst_t),         intent(inout) :: pf
    type(ndarray_oc), target,  intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
    integer,                   intent(in   ) :: nsteps
    integer,                   intent(inout) :: itersState, itersAdjoint
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(in   ) :: searchDir(:,:,:)  
    real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq
    real(pfdp),                intent(inout) :: stepSize, LinftyNormGrad, L2NormGradSq, gradient(:,:,:), savedAdjoint(:,:,:)
    logical,                   intent(inout) :: stepTooSmall
       
    real(pfdp) :: globObjNew, directionTimesGradient, globDirXGradNew
    integer    :: l, ierror
    !real(pfdp) :: armijoDecrease ! this should be in a structure global to the module, and set in
    !                             ! something like init_optimization, along with other parameters 
    !real(pfdp) :: minStepSize
    !armijoDecrease = 1e-4
    !minStepSize    = 1e-6

    do
      call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, stepSize)
      !restrict control
      do l = pf%nlevels-1,1,-1
	call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
      end do

      call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint)
      itersState = itersState + pf%state%itcnt
           
      call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
      if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew

      if (globObjNew < globObj + armijoDecrease*stepSize*globDirXGrad) then
        ! evaluate gradient to be consistent with wolfe_powell_step
        call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
        itersAdjoint = itersAdjoint + pf%state%itcnt
        return 
      end if

       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
       ! no need to restrict here, coarse levels get overwritten above
       stepSize = 0.5 * stepSize
       if (stepSize < minStepSize) then
         stepTooSmall = .true.
         if (pf%rank .eq. 0) print *, 'no step found, stopping'
         !call write_control(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
         return
       end if
     end do
   end subroutine armijo_step



!    subroutine zoom(stepLowIn, stepHighIn, stepSizeIn, stepSizeOut, &
!                    pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
!                    globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepTooSmall, globObjLowIn)
!      type(pf_pfasst_t),         intent(inout) :: pf
!      type(ndarray_oc), target,  intent(inout) :: q1
!      real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
!      integer,                   intent(in   ) :: nsteps
!      integer,                   intent(inout) :: itersState, itersAdjoint
!      logical,                   intent(in   ) :: predict
!      real(pfdp),                intent(in   ) :: searchDir(:,:) 
!      real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSizeOut
!      real(pfdp),                intent(in   ) :: stepSizeIn, stepLowIn, stepHighIn, globObjLowIn
!      real(pfdp),                intent(inout) :: gradient(:,:), savedAdjoint(:,:)
!      logical,                   intent(inout) :: stepTooSmall
! 
!      real(pfdp) :: directionTimesGradient, globObjNew, globDirXGradNew
!      real(pfdp) :: stepSize, stepLow, stepHigh, globObjLow
!      integer    :: l, ierror, iter
!  
!      stepSize = stepSizeIn
!      stepLow = stepLowIn
!      stepHigh = stepHighIn
!      globObjLow = globObjLowIn
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
!         call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint)
!         itersState = itersState + pf%state%itcnt
!         call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!         if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew
!       
!         if ( (globObjNew > globObj + armijoDecrease*stepSize*globDirXGrad) .or. &
!              (globObjNew >= globObjLow) )                                       then
!           stepHigh = stepSize
!           if (pf%rank .eq. 0) print *, 'set new stepHigh to stepSize; stepHigh = ', stepHigh, 'stepLow = ', stepLow
!         else          
!           call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
!           itersAdjoint = itersAdjoint + pf%state%itcnt
!           directionTimesGradient = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
!           call mpi_allreduce(directionTimesGradient, globDirXGradNew, 1, MPI_REAL8, MPI_SUM, &
!                              pf%comm%comm, ierror)          
!           if (pf%rank .eq. 0) print *, 'globDirXGradNew = ', globDirXGradNew
! 
!           
!           if (abs(globDirXGradNew) <= -c2*globDirXGrad) then
!             stepSizeOut = stepSize
!             return
!           end if
!    
!           if (globDirXGradNew*(stepHigh-stepLow) >= 0) then
!             if(pf%rank == 0) print *, stepSize, stepHigh, stepLow
!             stepHigh = stepLow
!             if (pf%rank .eq. 0) print *, 'set new stepHigh to stepLow; stepHigh = ', stepHigh, 'stepLow = ', stepLow, 'stepSize = ', stepSize 
!           end if
!           stepLow = stepSize
!           if (pf%rank .eq. 0) print *, 'set new stepLow; stepHigh = ', stepHigh, 'stepLow = ', stepLow, 'stepSize = ', stepSize 
!           globObjLow = globObjNew
!         end if
!       end do
!    end subroutine zoom
! 
! 
!    subroutine wolfe_powell_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
!                                 globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
!     type(pf_pfasst_t),         intent(inout) :: pf
!     type(ndarray_oc), target,  intent(inout) :: q1
!     real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
!     integer,                   intent(in   ) :: nsteps
!     integer,                   intent(inout) :: itersState, itersAdjoint
!     logical,                   intent(in   ) :: predict
!     real(pfdp),                intent(in   ) :: searchDir(:,:) 
!     real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq
!     real(pfdp),                intent(inout) :: stepSize, gradient(:,:), savedAdjoint(:,:)
!     logical,                   intent(inout) :: stepTooSmall
! 
!     integer    :: l, ierror, iter
!     real(pfdp) :: prevStepSize, directionTimesGradient, globDirXGradNew, globObjNew, globObjPrev
!     !, maxStepSize, minStepSize, stepIncFactor, armijoDecrease, c2, 
!     prevStepSize   = 0.0_pfdp
!     !maxStepSize    = 1e6
!     !minStepSize    = 1e-6
!     !stepIncFactor  = 10
!     !armijoDecrease = 1e-4
!     !c2             = 0.9 
!     globObjPrev = globObj
!     iter = 1
!     
!     do
!       call update_control(pf%levels(pf%nlevels)%ctx, searchDir, stepSize)
!       !restrict control
!       do l = pf%nlevels-1,1,-1
! 	call restrict_control(pf%levels(l)%ctx, pf%levels(l+1)%ctx)
!       end do
! !       call dump_control(pf%levels(pf%nlevels)%ctx, pf, 'u1')
! 
!       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint)
!       itersState = itersState + pf%state%itcnt
!       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!       if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew
!       
!       if ( (globObjNew > globObj + armijoDecrease*stepSize*globDirXGrad) .or. &
!            ( globObjNew >= globObjPrev .and. iter > 1) )                   then
!         call zoom(prevStepSize, stepSize, stepSize, stepSize, &
!                   pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, globObj, globDirXGrad, &
!                   objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepTooSmall, globObjPrev)
!         return
!       end if                    
!       
!       if (pf%rank .eq. 0) print *, 'evaluate new gradient'
!       call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
!       itersAdjoint = itersAdjoint + pf%state%itcnt
!       directionTimesGradient = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
!       call mpi_allreduce(directionTimesGradient, globDirXGradNew, 1, MPI_REAL8, MPI_SUM, &
!                          pf%comm%comm, ierror)
!       if (pf%rank .eq. 0) print *, 'globDirXGradNew = ', globDirXGradNew
!       if (abs(globDirXGradNew) <= -c2*globDirXGrad) then
!         return
!       end if
! 
!       if (globDirXGradNew >= 0) then
!         call zoom(stepSize, prevStepSize, stepSize, stepSize, &
!                   pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, globObj, globDirXGrad, &
!                   objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepTooSmall, globObjNew)
!         return
!       end if
! 
!       globObjPrev  = globObjNew
!       prevStepSize = stepSize
!       stepSize     = stepIncFactor*stepSize
!       iter = iter + 1
!       if (stepSize > maxStepSize) then
!          stepTooSmall = .true.
!          if (pf%rank .eq. 0) print *, 'no step found (> max), stopping'
!          !call write_control(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
!          return
!        end if
!        call update_control(pf%levels(pf%nlevels)%ctx, searchDir, -prevStepSize)
!     end do
! 
!    end subroutine wolfe_powell_step
  
  
  subroutine strong_wolfe_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
                          globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
    type(pf_pfasst_t),         intent(inout) :: pf
    type(ndarray_oc), target,  intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
    integer,                   intent(in   ) :: nsteps
    integer,                   intent(inout) :: itersState, itersAdjoint
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(in   ) :: searchDir(:,:,:) 
    real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq
    real(pfdp),                intent(inout) :: stepSize, gradient(:,:,:), savedAdjoint(:,:,:)
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

      call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
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
      
      call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, stepSize)
      do l = pf%nlevels-1,1,-1
	call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
      end do

      call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
      itersState = itersState + pf%state%itcnt
      call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
      if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew
      
      if (globObjNew >= globObj + armijoDecrease*stepSize*globDirXGrad) then
        call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
        high = stepSize
        cycle
      end if

      ! now Armijo is satisfied again
      call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
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

