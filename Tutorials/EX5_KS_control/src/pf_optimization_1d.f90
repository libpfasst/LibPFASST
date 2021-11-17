module pf_mod_optimization
  use pfasst
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_mpi
  use pf_mod_ndarray_oc
  use pf_mod_parallel_oc
  use my_sweeper

  implicit none

  real(pfdp), parameter, private :: armijoDecrease = 1e-4
  real(pfdp), parameter, private :: minStepSize    = 1e-6
  real(pfdp), parameter, private :: maxStepSize    = 1e6
  real(pfdp), parameter, private :: stepIncFactor  = 10
  real(pfdp), parameter, private :: c2             = 0.1 !0.9 
 
contains

  !> Evaluate objective function. For this, solve the state equation with the given control,
  !  evaluate the tracking term and the norm of the control
  subroutine evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objective, L2NormCtrlSq, savedStates, ctrl)
  
    type(pf_pfasst_t),         intent(inout) :: pf
    type(pf_ndarray_oc_t), target,  intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha
    integer,                   intent(in   ) :: nsteps
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(  out) :: objective, L2NormCtrlSq
    real(pfdp),                intent(inout) :: savedStates(:,:,:)    !nsteps, nnodes, nvars
    real(pfdp),                intent(inout) :: ctrl(:) ! initial condition, size=nvars
    integer             :: m, step, nnodes, thisstep, ierror, level_index, nx
    real(pfdp) :: stepobjective, roundedValue
    real(pfdp), pointer :: obj(:) !, savedStatesCompressed(:)
    

    nx = pf%levels(pf%nlevels)%lev_shape(1)
    nnodes = pf%levels(pf%nlevels)%nnodes
    allocate(obj(nnodes))    
    objective = 0
        
    pf%state%itcnt   = 0  ! reset iteration count here, not inside pf_pfasst_block_oc

    q1%yflatarray = ctrl
    q1%pflatarray = 0.0_pfdp
    
    if ((pf%rank .eq. 0) .or. predict) then
      call pf%levels(pf%nlevels)%q0%copy(q1, 1) !do this only for rank 0 or predict = true?
    end if
    if (pf%rank .eq. 0) print *, ' *********  solve state **************'

    pf%q0_style = 0
    if(.not. predict) pf%q0_style = 2

    do step = 1, nsteps
    
      thisstep = (step-1)*pf%comm%nproc + pf%rank
      pf%state%pfblock = step

      
      if(.not. predict) then ! warm start, load saved state values
        do m = 1, nnodes
          call pf%levels(pf%nlevels)%Q(m)%unpack(savedStates(step,m,:), 1)
        end do
        
        if (pf%rank .eq. 0) then
         call pf%levels(pf%nlevels)%Q(1)%copy(pf%levels(pf%nlevels)%q0, 1)
        else
         call pf%levels(pf%nlevels)%q0%copy(pf%levels(pf%nlevels)%Q(1), 1)
        end if
        ! re-evaluate function values
        call pf%levels(pf%nlevels)%ulevel%sweeper%evaluate_all(pf, pf%state%finest_level, &
                                    thisstep*dt+dt*pf%levels(pf%nlevels)%nodes, flags=1, step=thisstep+1)
      end if
    
      call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., 1, step=thisstep)      
!     end if
    
      ! evaluate objective on current step    
      do m = 1, nnodes
        call objective_function(pf%levels(pf%nlevels)%ulevel%sweeper, pf%levels(pf%nlevels)%Q(m), &
                                    pf%levels(pf%nlevels)%lev_shape(1), m, obj(m), step)


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
        call pf%levels(pf%nlevels)%q0%unpack(pf%levels(pf%nlevels)%send, 1)    !<  Everyone resets their q0

      end if
    end do
    
    objective = 0.5*objective !0.5 for trapezoidal rule
    
    ! L2(Omega) control norm
    L2NormCtrlSq = sum(ctrl**2)*Lx/dble(nx)
    
    objective = 0.5*objective + 0.5*alpha*L2NormCtrlSq
!     deallocate(savedStatesCompressed)
    deallocate(obj)
    
  end subroutine evaluate_objective
  
  !> Compute the gradient. For this, use the previously computed state solution to solve the adjoint equation.
  subroutine evaluate_gradient(pf, q1, dt, nsteps, predict, ctrl, alpha, gradient, LinftyNormGrad, L2NormGradSq, savedStates, savedAdjoint)
    type(pf_pfasst_t),        intent(inout) :: pf
    type(pf_ndarray_oc_t), target, intent(inout) :: q1
    real(pfdp),               intent(in   ) :: dt, alpha

    integer,                  intent(in   ) :: nsteps
    logical,                  intent(in   ) :: predict
    real(pfdp),               intent(  out) :: gradient(:), LinftyNormGrad, L2NormGradSq
    real(pfdp),               intent(inout) :: savedStates(:,:,:), savedAdjoint(:,:,:)
    real(pfdp),               intent(inout) :: ctrl(:)

    integer :: m, step, nnodes, thisstep, ierror
    real(pfdp), pointer :: obj(:) 

      
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

        ! load stored state solution
        do m = 1, nnodes
          call pf%levels(pf%nlevels)%Q(m)%unpack(savedStates(step,m,:), 1)
        end do       
        ! restrict saved states to coarse levels
        call restrict_for_adjoint(pf, thisstep*dt, dt, 1) 
        
        ! when using warm starts, load stored adjoint solution from previous optimization iteration

        if(.not. predict) then
          do m = 1, nnodes
            call pf%levels(pf%nlevels)%Q(m)%unpack(savedAdjoint(step,m,:), 2)
           end do
        
          if(pf%rank==pf%comm%nproc-1) then 
              ! this happens in sweep anyway, but not in restriction of q0/qend at start of predictor
            call pf%levels(pf%nlevels)%Q(nnodes)%copy(pf%levels(pf%nlevels)%qend, 2)
          else
            call pf%levels(pf%nlevels)%qend%copy(pf%levels(pf%nlevels)%Q(nnodes), 2)
          end if
          ! re-evaluate function values
          call pf%levels(pf%nlevels)%ulevel%sweeper%evaluate_all(pf, pf%state%finest_level, &
                                    thisstep*dt+dt*pf%levels(pf%nlevels)%nodes, flags=2, step=thisstep+1)
        end if

        ! iterate on current block of time steps
        call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., 2, step=thisstep) 

        ! save computed adjoint
        do m = 1, pf%levels(pf%nlevels)%nnodes
          call pf%levels(pf%nlevels)%Q(m)%pack(savedAdjoint(step, m, :), 2)          
        end do
                                    
        ! if we have more time steps to compute distribute new terminal condition
        if(step > 1) then
          call pf%levels(pf%nlevels)%q0%pack(pf%levels(pf%nlevels)%send, 2)    !<  Pack away your last solution
          call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, 0)
          call pf%levels(pf%nlevels)%qend%unpack(pf%levels(pf%nlevels)%send, 2)    !<  Everyone resets their qend
        end if
     end do

    ! evaluate gradient; as is not time dependent, compute on first rank only (as it requires p(0)) and broadcast.
    if(pf%rank == 0) then
        ! gradient is adjoint at t0 (plus regularization term derivative = alpha*ctrl)
        ! pack away on rank 0
        call pf%levels(pf%nlevels)%q0%pack(pf%levels(pf%nlevels)%send, 2)
    end if
    ! broadcast
    call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, 0)
    ! everybody sets gradient (actually this could be skipped, only rank 0 needs initial condition)
    ! need to check implementation
    gradient = pf%levels(pf%nlevels)%send
    gradient = gradient + alpha * ctrl
    
    pf%q0_style = 0     
  end subroutine evaluate_gradient
  
  !> Compute admissible step size using the Armijo rule. This evaluates the objective, and if a suitable 
  !  step size is found the new gradient is evaluated as well (to be consistend with the strong Wolfe step
  !  routine, where the gradient is required.
   subroutine armijo_step(pf, q1, dt, nsteps, itersState, itersAdjoint, skippedState, predict, searchDir, gradient, &
                        savedStatesFlat, savedAdjointsFlat, ctrl, alpha, globObj, &
                        globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
    type(pf_pfasst_t),         intent(inout) :: pf
    type(pf_ndarray_oc_t), target,  intent(inout) :: q1
    real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
    integer,                   intent(in   ) :: nsteps
    integer,                   intent(inout) :: itersState, itersAdjoint, skippedState
    logical,                   intent(in   ) :: predict
    real(pfdp),                intent(in   ) :: searchDir(:)
    real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq
    real(pfdp),                intent(inout) :: stepSize, LinftyNormGrad, L2NormGradSq, gradient(:), &
                                                savedStatesFlat(:,:,:), savedAdjointsFlat(:,:,:)
    logical,                   intent(inout) :: stepTooSmall
    real(pfdp),                intent(inout) :: ctrl(:)


    real(pfdp) :: globObjNew, directionTimesGradient, globDirXGradNew
    integer    :: l, ierror
    !real(pfdp) :: armijoDecrease ! this should be in a structure global to the module, and set in
    !                             ! something like init_optimization, along with other parameters
    !real(pfdp) :: minStepSize
    !armijoDecrease = 1e-4
    !minStepSize    = 1e-6

    do
      ctrl = ctrl + stepsize*searchDir
!       call update_control(pf%levels(pf%state%finest_level)%ulevel%sweeper, searchDir, stepSize)
!       !restrict control
!       do l = pf%state%finest_level-1,1,-1
! 	      call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
!       end do

      call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStatesFlat, &
                              ctrl)

!       itersState = itersState + pf%state%itcnt
      if(pf%state%finest_level == pf%nlevels) & !count only finest level sweeps
        itersState = itersState + pf%state%itcnt - pf%state%skippedy
      skippedState = skippedState + pf%state%skippedy

      call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!      globObjNew = objectiveNew

      if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew

      if (globObjNew < globObj + armijoDecrease*stepSize*globDirXGrad) then
        ! evaluate gradient to be consistent with wolfe_powell_step
        call evaluate_gradient(pf, q1, dt, nsteps, predict, ctrl, alpha, gradient, LinftyNormGrad, L2NormGradSq, savedStatesFlat, &
                               savedAdjointsFlat)

        if(pf%state%finest_level == pf%nlevels) & !count only finest level sweeps
          itersAdjoint = itersAdjoint + pf%state%itcnt
        return
      end if

      ctrl = ctrl - stepsize*searchDir

!        call update_control(pf%levels(pf%state%finest_level)%ulevel%sweeper, searchDir, -stepSize)
       ! no need to restrict here, coarse levels get overwritten above
       stepSize = 0.5 * stepSize
       if (stepSize < minStepSize) then
         stepTooSmall = .true.
         if (pf%rank .eq. 0) print *, 'no step found, stopping'
         !call write_control(pf%levels(pf%state%finest_level)%ctx, k, "u_sdc_split_final")
         return
       end if
     end do
   end subroutine armijo_step

   
   !> Compute admissible step size using the strong Wolfe criterion. This evaluates the objective, and the new !   gradient. 
   subroutine strong_wolfe_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, &
        savedStates, savedAdjoint, ctrl, &
        alpha, globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, &
        stepSize, stepTooSmall)
     type(pf_pfasst_t),         intent(inout) :: pf
     type(pf_ndarray_oc_t), target,  intent(inout) :: q1
     real(pfdp),                intent(in   ) :: dt, alpha, globObj, globDirXGrad
     integer,                   intent(in   ) :: nsteps
     integer,                   intent(inout) :: itersState, itersAdjoint
     logical,                   intent(in   ) :: predict
     real(pfdp),                intent(in   ) :: searchDir(:) 
     real(pfdp),                intent(  out) :: objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq
     real(pfdp),                intent(inout) :: stepSize, gradient(:), savedStates(:,:,:), savedAdjoint(:,:,:), ctrl(:)
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
       ctrl = ctrl + high*searchDir
!       call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, high) !stepSize)
       !restrict control
!        do l = pf%nlevels-1,1,-1
!        call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
!        end do
 
       !call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStates, ctrl)
       first = .false.
       
       itersState = itersState + pf%state%itcnt
       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
       if(pf%rank == 0) print *, high, 'objectiveNew (L2) = ', globObjNew ! *, stepSize,
 
       ctrl = ctrl - high*searchDir
!        call update_control(pf%levesl(pf%nlevels)%ulevel%sweeper, searchDir, -high) !-stepSize)
       
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
       
       ctrl = ctrl + stepSize*searchDir    
!        call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, stepSize)
!        do l = pf%nlevels-1,1,-1
!         call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
!        end do
 
       call evaluate_objective(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedStates, ctrl)
       !(pf, q1, dt, nsteps, predict, alpha, objectiveNew, L2NormUSq, savedAdjoint, first)
       itersState = itersState + pf%state%itcnt
       call mpi_allreduce(objectiveNew, globObjNew, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
       if(pf%rank == 0) print *, stepSize, 'objectiveNew (L2) = ', globObjNew
       
       if (globObjNew >= globObj + armijoDecrease*stepSize*globDirXGrad) then
         ctrl = ctrl - stepSize*searchDir
         !call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
         high = stepSize
         cycle
       end if
 
       ! now Armijo is satisfied again
       !call evaluate_gradient(pf, q1, dt, nsteps, predict, gradient, LinftyNormGrad, L2NormGradSq, savedAdjoint)
       call evaluate_gradient(pf, q1, dt, nsteps, predict, ctrl, alpha, gradient, LinftyNormGrad, L2NormGradSq, savedStates, savedAdjoint)
       itersAdjoint = itersAdjoint + pf%state%itcnt
       directionTimesGradient = compute_scalar_prod(searchDir, gradient)
       
       ! gradient not time dependent, so scalar prod is same everywhere
       globDirXGradNew = directionTimesGradient
       !call mpi_allreduce(directionTimesGradient, globDirXGradNew, 1, MPI_REAL8, MPI_SUM, &
       !                  pf%comm%comm, ierror)          
                          
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
       
       ctrl = ctrl - stepSize*searchDir     
!        call update_control(pf%levels(pf%nlevels)%ulevel%sweeper, searchDir, -stepSize)
     end do
     
   end subroutine strong_wolfe_step
   
end module pf_mod_optimization

