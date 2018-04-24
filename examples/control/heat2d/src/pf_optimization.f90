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
    integer             :: m
    real(pfdp), pointer :: obj(:)
    real(pfdp) :: t(pf%levels(pf%nlevels)%nnodes)
    logical    :: predictAdjoint
    
    predictAdjoint = .true.
    if (present(predictAdj) ) predictAdjoint = predictAdj
   
    solve_y = .true.

    if (pf%rank .eq. 0) then
      call initial(q1, 0.0_pfdp, dt) !dt should be given by rank !initial_rd
    else
      q1%yflatarray = 0.0_pfdp
    end if
    if ((pf%rank .eq. 0) .or. predict) then
      call pf%levels(pf%nlevels)%q0%copy(q1, 1) !do this only for rank 0 or predict = true?
    end if
    
!     if (do_mixed .eq. 1 ) then
!       q1%pflatarray = 0.0_pfdp
!       call pf%levels(pf%nlevels)%qend%copy(q1, 2)
!       if( predict .eqv. .false. .and. predictAdjoint .eqv. .true. ) then
!         ! call predictor just for adjoint part
!         call pf_predictor(pf, pf%rank*dt, dt, 2)
! !           do m = 1, pf%levels(pf%nlevels)%nnodes
! !             call pf%levels(pf%nlevels)%encap%unpack(pf%levels(pf%nlevels)%Q(m), savedAdjoint(m,:), 2)
! !           end do
! ! 
! !           t = pf%rank*dt + dt*pf%levels(pf%nlevels)%nodes
! !           call pf%levels(pf%nlevels)%sweeper%evaluate_all(pf%levels(pf%nlevels), t, 2)       
!       end if
!     end if

    if (pf%rank .eq. 0) print *, ' *********  solve state **************'
    
!     if (do_mixed .eq. 1) then
!        call pf_pfasst_block(pf, dt, nsteps, predict, 0)
!        ! need to save values in adjoint component of Q(m)!
! !        if(pf%rank .eq. 0) print *, 'saving adjoint'
!        do m = 1, pf%levels(pf%nlevels)%nnodes
!          call pf%levels(pf%nlevels)%Q(m)%pack(savedAdjoint(m,:), 2)
!        end do
! !        print *, pf%rank, maxval(savedAdjoint(1,:))
!     else
       call pf_pfasst_block_oc(pf, dt, nsteps, predict, 1)
!     end if
    
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
    integer :: m
    integer :: dest, source, ierror, stat(MPI_STATUS_SIZE)
    real(pfdp), allocatable :: tmp(:)
    real(pfdp) :: t(pf%levels(pf%nlevels)%nnodes)

    solve_y = .false.      

    allocate(tmp(product(pf%levels(pf%nlevels)%shape)))

    !call initial(q1, pf%rank*dt, (pf%rank+1)*dt)
    q1%pflatarray = 0.0_pfdp ! only if no final time term in objective, otherwise this is nonzero, and terminal needs to be added to this
    
    call pf%levels(pf%nlevels)%qend%copy(q1, 2)
    ! do this only for final step or predict = true?
    !call restrict_for_adjoint(pf, 1)

    
    if (pf%rank .eq. 0) print *, '*********  solve adjoint *************'
    if(predict) pf%q0_style = 1
    call pf_pfasst_block_oc(pf, dt, nsteps, predict, 2) !predict
    pf%q0_style = 0

    if (pf%rank .eq. 0) print *, '*********  compute gradient *************'
    
    do m = 1, pf%levels(pf%nlevels)%nnodes
      call pf%levels(pf%nlevels)%Q(m)%pack(tmp, 2)
      gradient(m,:,:) = reshape(tmp, (/pf%levels(pf%nlevels)%shape(1), pf%levels(pf%nlevels)%shape(2)/))
    end do

    call construct_gradient(pf%levels(pf%nlevels)%ulevel%sweeper, gradient, pf%levels(pf%nlevels)%nodes, &
                                      LinftyNormGrad, L2NormGradSq)
                                      
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

