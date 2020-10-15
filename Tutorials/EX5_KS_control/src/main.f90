program main
  use pfasst
  use pf_mod_mpi
  use pf_mod_ndarray_oc
  use pf_mod_restrict

  use pf_mod_optimization
  use my_sweeper
  use my_level
  use hooks
  use probin

  use fnpy

implicit none

  type(pf_comm_t)          :: comm
  type(pf_pfasst_t)        :: pf
  class(my_sweeper_t), pointer   :: my_sweeper_ptr     !<  pointer to SDC sweeper


  type(pf_ndarray_oc_t) :: q1    !< for setting initial conditions
  integer        :: ierror, l, m, k, i, p
  character(len = 64) :: fout
  character(len =128) :: logfilename
  character(len = 64) :: fname
  character(len=256)     :: npyfname
  character(len = 64) :: probin_fname
  character(len = 64) :: shell_cmd
  integer        :: iout,nout

  real(pfdp), pointer :: gradient(:), prevGrad(:), &    !< gradient of current and previous optimization iter
                         searchDir(:), prevSearchDir(:),& !< current and previous search direction
                         ctrl(:), & !< control (to be determined by the optimization), here: initial condition
                         savedStatesFlat(:,:,:), savedAdjointsFlat(:,:,:) !< to store solutions, e.g., for warm starting
  real(pfdp)          :: LinftyNormGrad, LinftyNormCtrl, L2NormCtrlSq, L2NormGradSq, & !< norms of gradient and control
                         objective, objectiveNew, & !< objective function values (current and after update for step size selection)
                         stepSize, prevStepSize, & !< current initial step size and previously chosen step size for updating the control
                         directionTimesGradient, & !< scalar product of search direction and gradient, used in step size selection
                         globObj, globObjNew, globDirXGrad, prevGlobDirXGrad, globL2NormGradSq, globLinftyNormGrad, & !< accumulated quantities (summed over all ranks)
                         beta, & !< determines the type of optimization method
                         num, denom, & !< helper variables to determine beta
                         L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfEXactCtrl, & !< to determine error in computed control
                         abs_res_tol_save, rel_res_tol_save !< to solve state for data generation with potentially stricter residual tolerance than given in the nml file for the optimization run

  logical             :: stepTooSmall, & !< whether step size determined by alg is too small
                         predict !< use predictor?
  integer             :: itersState, itersAdjoint, skippedState, sumSkippedState, sumItersState, sumItersAdjoint, &  !< counting PFASST iterations on different equations for statistics
                         step, & ! loop counter for time steps
                         nsteps_per_rank, & !< how many timesteps for each CPU?
                         thisstep !< what is the true time step corresponding the the CPUs local step number

  character(8)   :: date  !< wall clock times for timing the optimization run
  character(10)  :: time  !< wall clock times for timing the optimization run
  integer        :: time_start(8), time_end(8)   !< wall clock times for timing the optimization run
  real(pfdp)     :: time_start_sec, time_end_sec !< wall clock times for timing the optimization run

  real(pfdp), allocatable :: targetState(:,:,:), & !< for setting up the tracking type objective function
                             npyOutput(:) !< helper arry for numpy output
  real, allocatable :: rndArr(:)  !< for generating random numbers for perturbations

  !
  ! read options
  !
  if (command_argument_count() >= 1) then
     call get_command_argument(1, value=probin_fname)
  else
     probin_fname = "probin.nml"
  end if
  call probin_init(probin_fname)


  !
  ! initialize mpi
  !
  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."

  !
  ! initialize pfasst
  !
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, comm,  fname=pfasst_nml)

  ! how many time steps for each rank to compute
  nsteps_per_rank = nsteps / comm%nproc
  if(nsteps_per_rank*comm%nproc /= nsteps) stop "ERROR: nsteps not an integer multiple of nproc"

  ! set up levels
  do l = 1, pf%nlevels

       pf%levels(l)%nnodes = nnodes(l)

       !  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)
       allocate(pf_ndarray_oc_factory_t::pf%levels(l)%ulevel%factory)

       !  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,[nvars(l)])
  end do

  call pf_pfasst_setup(pf)

  pf%state%finest_level = pf%nlevels

  pf%outdir = output

  ! do we want to save timings?
  pf%save_timings = 0

  ! saving other stuff is not yet tested for optimal control problems
  pf%save_residuals = .false.
  pf%save_errors=.false.
  pf%save_delta_q0=.false.
  pf%save_solutions=0

  call initialize_results(pf)

  ! output some information each after time step
!   call pf_add_hook(pf,pf%nlevels,PF_POST_BLOCK,echo_error_hook)


  if(pf%rank == 0) print *, "computing ", nsteps_per_rank, "steps per CPU"

  !  Make directory for Data if it does not exist
  call system('if [ ! -e ./Dat ]; then mkdir Dat; fi')
  shell_cmd = 'if [ ! -e ./Dat/'//trim(fbase)//' ]; then mkdir Dat/'//trim(fbase)//'; fi'
  call system(shell_cmd)

  if(write_numpy) then
     call system('if [ ! -e ./npy ]; then mkdir npy; fi')
  end if

  ! open output files
  write (fname, "(A,I0.2,A3,I0.3,A6,I0.3,A6,I0.3)") 'Niter',pf%niters,'_Nx',nvars(pf%nlevels),'_Nstep',nsteps,'_Nproc',comm%nproc
  foutbase = 'Dat/'//trim(fbase)//'/'//trim(fname)

  if (warmstart .eq. 1) then
    predict = .false.
    logfile = trim(logfile)//'_warm'
  else
    predict = .true.
    logfile = trim(logfile)//'_cold'
  endif

  write(logfilename, "(A,'_tol',i0.3,'_optiter',i0.4,'_Nstep',i0.3,'_Nproc',i0.3,'.log')") trim(logfile), test_no, max_opt_iter, nsteps, comm%nproc

  !  Output the run parameters
  if (pf%rank == 0) then
     call pf_print_options(pf, 6,.TRUE.)
     fout = 'Dat/'//trim(fbase)//'/'//trim(fname)//'_params.m'
!     print*, fout
     open(unit=103, file = fout, status = 'unknown', action = 'write')
     do iout=1,2
        if (iout .eq. 1) then
           nout = 103
        else
           nout = 6
        endif

        write(nout,*) 'Nproc=',comm%nproc
        write(nout,*) 'Ndim=',ndim
        write(nout,*) 'Nnodes=',nnodes(1:pf%nlevels)
        write(nout,*) 'Nvars=',nvars(1:pf%nlevels)
        write(nout,*) 'Finterp=',Finterp
        write(nout,*) 'nsteps=',nsteps
        write(nout,*) 'nsteps_per_rank=',nsteps_per_rank
        write(nout,*) 'dt=',dt
        write(nout,*) 'fbase=',fbase
        write(nout,*) 'logfile=',logfilename
        if (do_imex .eq. 1) then
          write(nout,*) 'use IMEX sweeper'
        else
          write(nout,*) 'use MISDC sweeper'
        endif
        if (warmstart .eq. 1) then
          write(nout,*) 'warm start with previous state solution'
        else
          write(nout,*) 'cold start with spread initial value'
        endif
          write(nout,*) 'regularization: alpha=', alpha
          write(nout,*) 'maximum number of optimization iterations=', max_opt_iter
          write(nout,*) 'optimization tolerances (gradient, objective)=', tol_grad, tol_obj
     end do
  end if

  ! create array for initial condition
  call ndarray_oc_build(q1, pf%levels(pf%nlevels)%lev_shape)

  ! set up and initialize optimization specific elements of the sweeper
  do l = 1, pf%nlevels
     call initialize_ocp(pf%levels(l)%ulevel%sweeper, pf, l, nsteps_per_rank)
  end do

  ! get exact control = true initial condition
  allocate(ctrl(nvars(pf%nlevels)))
  call initial_sol(ctrl)

  allocate(npyOutput(nvars(pf%nlevels)))

  ! numpy output
  if (write_numpy) then
     write(npyfname, "(A,'.npy')") "npy/uexact"
     call save_double(npyfname, shape(ctrl), real(ctrl(:),8))
  end if

  ! solve equation with correct initial condition to generate tracking data
  allocate(targetState(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes,nvars(pf%nlevels)))
  targetState = 0.0_pfdp

  pf%q0_style = 0   ! used in the predictor

  ! set initial condition
  if (pf%rank .eq. 0) then
    q1%yflatarray = ctrl
  else
    q1%yflatarray = 0.0_pfdp
  end if
  call pf%levels(pf%nlevels)%q0%copy(q1, 1)


  if (pf%rank .eq. 0)  then
    print *, ' **** solve state with exact initial condition ***'
    if (write_numpy) then
       write(npyfname, "(A,'s',i0.4'.npy')") "npy/y", 0
       call pf%levels(pf%nlevels)%q0%pack(npyOutput,1)
       call save_double(npyfname, shape(npyOutput), real(npyOutput,8))
    end if
  end if

  ! PFASST loop to solve state equation
  ! done manually, as we need to keep track of some quantities (e.g., save state solution values),
  ! so we cannot use the standard PFASST controller magic to do it for us
  do step=1, nsteps_per_rank
      thisstep = (step-1)*pf%comm%nproc + pf%rank
      pf%state%pfblock = step
      call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., 1, step=thisstep+1)

      ! pack true solution to targetState for setting up the tracking type objective function
      do m=1, pf%levels(pf%nlevels)%nnodes
         call pf%levels(pf%nlevels)%Q(m)%pack(targetState(step,m,:),1)

         if (write_numpy) then
            write(npyfname, "(A,'s',i0.4,'m',i0.2,'.npy')") "npy/ytarget", step, m
            call save_double(npyfname, shape(targetState(step,m,:)), real(targetState(step,m,:),8))
         end if

         ! TODO: add some observation noise here
      end do

      ! output state at time step end
      if (write_numpy) then
         write(npyfname, "(A,'s',i0.4'.npy')") "npy/y", thisstep+1
         call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%pack(npyOutput,1)
         call save_double(npyfname, shape(npyOutput), real(npyOutput,8))
      end if

      ! if not done pass on initial condition
      if( step < nsteps_per_rank) then
        call pf%levels(pf%nlevels)%qend%pack(pf%levels(pf%nlevels)%send, 1)    !<  Pack away your last solution
        call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, pf%comm%nproc-1)
        call pf%levels(pf%nlevels)%q0%unpack(pf%levels(pf%nlevels)%send, 1)    !<  Everyone resets their q0
      end if
  end do

  ! set up tracking objective
  call set_ydesired(pf%levels(pf%nlevels)%ulevel%sweeper, targetState)
  do l = pf%nlevels-1,1,-1
    call restrict_ydesired(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
  end do

  ! wait for everybody to finish their setup/data generation phase
  call mpi_barrier(pf%comm%comm, ierror)


  ! set up remaining data structures for the optimization loop
  if(pf%rank == 0) &
     open(unit=105, file = logfilename , &
         status = 'unknown',  action = 'write')

  if(pf%rank == 0) write(105,*) "# iter   L2_grad   objective   stepsize"

  allocate(gradient(nvars(pf%nlevels)))
  allocate(prevGrad(nvars(pf%nlevels)))
  allocate(searchDir(nvars(pf%nlevels)))
  allocate(prevSearchDir(nvars(pf%nlevels)))
  allocate(savedStatesFlat(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, product(pf%levels(pf%nlevels)%lev_shape)))
  allocate(savedAdjointsFlat(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, product(pf%levels(pf%nlevels)%lev_shape)))

  gradient = 0
  prevGrad = 0
  prevGlobDirXGrad = 0.0_pfdp
  globDirXGrad = 0.0_pfdp
  prevStepSize = 1.0_pfdp
  prevSearchDir = 0
  savedStatesFlat = 0.0_pfdp
  savedAdjointsFlat = 0.0_pfdp

  itersState = 0
  itersAdjoint = 0
  skippedState = 0


  ! initial guess for control: some fraction or perturbation of true control
  call initial_sol(ctrl)

  ctrl = 0.75 * ctrl
  ! allocate(rndArr(nvars(pf%nlevels)))
  ! call random_number(rndArr)
  ! ctrl = ctrl + 0.1*rndArr
  ! call pf_broadcast(pf, ctrl, pf%levels(pf%nlevels)%mpibuflen, 0)
  ! deallocate(rndArr)


  call date_and_time(date=date, time=time, values=time_start)
  if (pf%rank .eq. 0) &
    print *, 'start optimization on ', date, ', ',  time

  ! optimization loop
  do k = 1, max_opt_iter

     ! output current control
     if (write_numpy) then
        write(npyfname, "(A,'k',i0.4,'.npy')") "npy/u", k
        call save_double(npyfname, shape(ctrl), real(ctrl(:),8))
     end if

     if(pf%rank == 0) print *, '===============Optimization ITERATION================', k

!      call mpi_barrier(pf%comm%comm, ierror)
!      if ( k .eq. 1 ) then
        call evaluate_objective(pf, q1, dt, nsteps_per_rank, .true., alpha, objective, L2NormCtrlSq, savedStatesFlat, ctrl)
        itersState = itersState + pf%state%itcnt - pf%state%skippedy
        skippedState = skippedState + pf%state%skippedy

     call mpi_allreduce(objective, globObj, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
     if(pf%rank == 0) print *, k, 'objective (L2) = ', globObj
     !exit

!      call mpi_barrier(pf%comm%comm, ierror)
     if ( k .eq. 1 ) then   ! in later iterations, this is done in linesearch
       call evaluate_gradient(pf, q1, dt, nsteps_per_rank, .true., ctrl, alpha, gradient, LinftyNormGrad, L2NormGradSq, savedStatesFlat, &
                              savedAdjointsFlat)
       itersAdjoint = itersAdjoint + pf%state%itcnt
     end if


     if(pf%rank == 0 .and. write_numpy) then
        write(npyfname, "(A,'k',i0.4,'.npy')") "npy/gradient", k
        call save_double(npyfname, shape(gradient), real(gradient(:),8))
     end if


     ! add regularization term derivative to gradient
     globL2NormGradSq = sum(gradient**2)*Lx/dble(nvars(pf%nlevels))
     globLinftyNormGrad = maxval(abs(gradient))

     ! output some information on the optimization iteration
     if(pf%rank == 0) print *, k, 'gradient (L2, Linf) = ', sqrt(globL2NormGradSq), globLinftyNormGrad
     if(pf%rank == 0) write(105,*) k, sqrt(globL2NormGradSq), globObj, prevStepSize

     ! are we done optimizing?
     if (sqrt(globL2NormGradSq) < tol_grad) then
       if(pf%rank == 0) print *, 'optimality condition satisfied (gradient norm small enough), stopping'
       exit
     end if
     if (globObj < tol_obj) then
       if(pf%rank == 0) print *, 'optimality condition satisfied (objective function small enough), stopping'
       exit
     end if

     ! determine search direction depending on the optimization method
     if ( k .eq. 1 ) then
        beta = 0.0
     else
        if (opt_method == 1) then
         !PR:
           num   = compute_scalar_prod(gradient, gradient-prevGrad)
           denom = compute_scalar_prod(prevGrad, prevGrad)
           beta = num/denom
        else if (opt_method == 2) then
!         !DY:
          num   = compute_scalar_prod(gradient, gradient)
          denom = compute_scalar_prod(gradient-prevGrad, prevSearchDir)
          beta = num/denom
        else if (opt_method == 3) then
!         !FR:
          num   = compute_scalar_prod(gradient, gradient)
          denom = compute_scalar_prod(prevGrad, prevGrad)
          beta = num/denom
        else
          ! none of the above, use steepest descent, i.e., beta = 0
          beta = 0
        end if
     end if
     searchDir = -gradient + beta*prevSearchDir

     ! step size control
     globDirXGrad = compute_scalar_prod(searchDir, gradient)
     if(pf%rank == 0) print *, k, 'gradient times searchDir = ', globDirXGrad

     stepSize = min(2.0_pfdp*prevStepSize, 1.0_pfdp)

     ! update for next iteration, as these get re-assigned during linesearch
     prevGlobDirXGrad = globDirXGrad
     prevSearchDir = searchDir
     prevGrad = gradient

     stepTooSmall = .false.

     if(pf%rank == 0) print *, k, 'step size = ', stepSize
     ctrl = ctrl + stepsize*searchDir

     if(use_wolfe) then
        call strong_wolfe_step(pf, q1, dt, nsteps_per_rank, itersState, itersAdjoint, predict, searchDir, gradient, &
        savedStatesFlat, savedAdjointsFlat, ctrl, alpha, &
        globObj, globDirXGrad, objectiveNew, L2NormCtrlSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
     else
        call armijo_step(pf, q1, dt, nsteps_per_rank, itersState, itersAdjoint, skippedState, predict, searchDir, gradient, &
        savedStatesFlat, savedAdjointsFlat, ctrl, alpha, &
        globObj, globDirXGrad, objectiveNew, L2NormCtrlSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
     end if

     if (stepTooSmall) exit

     ! update for next iteration
     prevStepSize = stepSize
   end do

   ! wrap up and output some information
   call mpi_barrier(pf%comm%comm, ierror)
   call date_and_time(date=date, time=time, values=time_end)
   if (pf%rank .eq. 0) then
     print *, 'end optimization on ', date, ', ',  time
     time_start_sec = time_start(3) * 24 * 3600 + time_start(5) * 3600 + time_start(6) * 60 &
           + time_start(7) + 0.001 * time_start(8)
      time_end_sec = time_end(3) * 24 * 3600 + time_end(5) * 3600 + time_end(6) * 60 &
           + time_end(7) + 0.001 * time_end(8)
      print *, 'duration [s]: ', time_end_sec-time_start_sec
      if (time_end(2) /= time_start(2) ) then
          print *, 'run finished in a different month than it started, duration needs to be corrected!'
      end if
   end if

  print *, 'rank:', pf%rank, 'total iterations state, adjoint, skipped state', itersState, itersAdjoint, skippedState

   sumItersState = 0
   sumItersAdjoint = 0
   sumSkippedState = 0
   call mpi_allreduce(itersState,   sumItersState,   1, MPI_INTEGER, MPI_SUM, pf%comm%comm, ierror)
   call mpi_allreduce(itersAdjoint,   sumItersAdjoint,   1, MPI_INTEGER, MPI_SUM, pf%comm%comm, ierror)
   call mpi_allreduce(skippedState,   sumSkippedState,   1, MPI_INTEGER, MPI_SUM, pf%comm%comm, ierror)

   if( pf%rank == 0 ) then
      print *, 'overall sum iterations state, adjoint, skipped', sumItersState, sumItersAdjoint, sumSkippedState
  !    print *, 'absolute error in computed control:  L2 = ', sqrt(globL2errorCtrl)
  !    print *, '                                   Linf = ', globLinfErrorCtrl
   !   print *, 'relative error in computed control:  L2 = ', sqrt(globL2errorCtrl)/sqrt(globL2exactCtrl)
   !   print *, '                                   Linf = ', globLinfErrorCtrl/globLinfExactCtrl
!       print *, 'exact control:                       L2 = ', sqrt(globL2exactCtrl)
!       print *, '                                   Linf = ', globLinfExactCtrl
   end if

  if (write_numpy) then
     write(npyfname, "(A,'.npy')") "npy/ufinal"
     call save_double(npyfname, shape(ctrl), real(ctrl(:),8))

     do step=1, nsteps_per_rank
        thisstep = (step-1)*pf%comm%nproc + pf%rank

        write(npyfname, "(A,'s',i0.4'.npy')") "npy/yfinal", thisstep+1
        call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%pack(npyOutput,1)
        call save_double(npyfname, shape(npyOutput), real(npyOutput,8))

        if(pf%rank == 0) then
           ! for completeness write y at t=0, i.e., the initial condition, as well
           write(npyfname, "(A,'s',i0.4'.npy')") "npy/yfinal", 0
           call save_double(npyfname, shape(ctrl), real(ctrl,8))
        end if
     end do
  end if

  !
  ! cleanup
  !
  deallocate(targetState)
  deallocate(npyOutput)
  deallocate(gradient)
  deallocate(prevGrad)
  deallocate(searchDir)
  deallocate(prevSearchDir)
  deallocate(savedStatesFlat)
  deallocate(savedAdjointsFlat)

  call ndarray_oc_destroy(q1)

  call pf_pfasst_destroy(pf)
  call mpi_finalize(ierror)

end program main
