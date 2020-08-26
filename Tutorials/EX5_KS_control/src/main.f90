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


  type(pf_ndarray_oc_t) :: q1, qend
  integer        :: ierror, iprovided, l, m, k, i, p
  character(len = 64) :: fout
  character(len =128) :: logfilename
  character(len = 64) :: fname
  character(len=256)     :: npyfname
  character(len = 64) :: probin_fname
  character(len = 64) :: shell_cmd
  integer        :: iout,nout

  real(pfdp), pointer :: gradient(:), prevGrad(:), searchDir(:), prevSearchDir(:), ctrl(:), &
                         savedStatesFlat(:,:,:), savedAdjointsFlat(:,:,:)
  real(pfdp)          :: LinftyNormGrad, LinftyNormCtrl, objective, objectiveNew, L2NormCtrlSq, L2NormGradSq, &
                         dx, stepSize, prevStepSize, directionTimesGradient, beta, &
                         globObj, globObjNew, globDirXGrad, prevGlobDirXGrad, globL2NormGradSq, tolGrad, tolObj, &
                         num, denom, globNum, globDenom, globLinftyNormGrad, &
                         L2errorCtrl, LinfErrorCtrl, globL2errorCtrl, globLinfErrorCtrl, &
                         L2exactCtrl, LinfEXactCtrl, globL2exactCtrl, globLinfExactCtrl, &                                                  
                         abs_res_tol_save, rel_res_tol_save
  logical             :: stepTooSmall, predict
  integer             :: itersState, itersAdjoint, skippedState, sumSkippedState, sumItersState, sumItersAdjoint, &
                         step, nsteps_per_rank, root
  character(8)   :: date
  character(10)  :: time
  integer        :: time_start(8), time_end(8)
  real(pfdp)     :: time_start_sec, time_end_sec

  real(pfdp), allocatable :: targetState(:,:,:), npyOutput(:,:,:)

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
!   call mpi_init_thread(mpi_thread_funneled, iprovided, ierror)
  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."

  !
  ! initialize pfasst
  !
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, comm,  fname=pfasst_nml)

  nsteps_per_rank = nsteps / comm%nproc
  if(nsteps_per_rank*comm%nproc /= nsteps) stop "ERROR: nsteps not an integer multiple of nproc"


  pf%save_timings = 0
  pf%save_residuals = .false.

  do l = 1, pf%nlevels

       pf%levels(l)%nnodes = nnodes(l)

       !  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)
       allocate(pf_ndarray_oc_factory_t::pf%levels(l)%ulevel%factory)

       !  Allocate the shape array for level (here just one dimension)

        !  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !>  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,[nvars(l)])

!       call sweeper_setup(pf%levels(l)%ulevel%sweeper, pf%levels(l)%lev_shape, nnodes(l), nsteps_per_rank, comm%nproc)
  end do

  call pf_pfasst_setup(pf)

  pf%state%finest_level = pf%nlevels

!   call initialize_results(pf%results, nsteps, pf%niters, pf%comm%nproc, pf%nlevels,pf%rank)

!   pf%debug = .true.

  pf%outdir = output
!   pf%save_results = .false.


!   call pf_add_hook(pf,pf%nlevels,PF_POST_BLOCK,echo_error_hook)


  if(pf%rank == 0) print *, "computing ", nsteps_per_rank, "steps per CPU"

  !  Make directory for Data if it does not exist


  call system('if [ ! -e ./Dat ]; then mkdir Dat; fi')
  shell_cmd = 'if [ ! -e ./Dat/'//trim(fbase)//' ]; then mkdir Dat/'//trim(fbase)//'; fi'
  call system(shell_cmd)
  ! open output files
  write (fname, "(A,I0.2,A3,I0.3,A6,I0.3,A6,I0.3)") 'Niter',pf%niters,'_Nx',nvars(pf%nlevels),'_Nstep',nsteps,'_Nproc',comm%nproc
  foutbase = 'Dat/'//trim(fbase)//'/'//trim(fname)
!  print *,'foutbase=',foutbase

!   write (logfilename, "(A,I0.2,A3,I0.3,A6,I0.3,A6,I0.3)") 'iter',max_opt_iter,'_Nx',nvars(pf%nlevels),'_Nstep',nsteps,'_Nproc',comm%nproc
  if( do_mixed == 1 ) then
    logfile = trim(logfile)//'_mixed'
  end if

!   logfilename = trim(logfile)
  if (warmstart .eq. 1) then
    predict = .false.
    logfile = trim(logfile)//'_warm'
  else
    predict = .true.
    logfile = trim(logfile)//'_cold'
  endif
!   if (nsteps_per_rank > 1 ) predict = .true. ! for state, on the first block we could in principle warm start as we load the state values
                                             ! for the adjoint solve; for paraexp like integration, however, this will not be the case
                                             ! -> decide in evaluate_objective according to warmstart, step and do_mixed

  write(logfilename, "(A,'_tol',i0.3,'_optiter',i0.4,'_Nstep',i0.3,'_Nproc',i0.3,'.log')") trim(logfile), test_no, max_opt_iter, nsteps, comm%nproc

!  open(unit=101, file = foutbase, status = 'unknown', action = 'write')

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

  call ndarray_oc_build(q1, pf%levels(pf%nlevels)%lev_shape)
  do l = 1, pf%nlevels
     call initialize_ocp(pf%levels(l)%ulevel%sweeper, pf, l, nsteps_per_rank)
  end do

  allocate(ctrl(nvars(pf%nlevels)))
  call initial_sol(ctrl)  ! get exact control

  ! solve equation with correct initial condition to generate tracking data
  allocate(targetState(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes,nvars(pf%nlevels)))
  targetState = 0.0_pfdp

  pf%q0_style = 0

  if (pf%rank .eq. 0) then
    q1%yflatarray = ctrl
  else
    q1%yflatarray = 0.0_pfdp
  end if
  call pf%levels(pf%nlevels)%q0%copy(q1, 1)
  if (pf%rank .eq. 0)  print *, ' **** solve state with exact initial condition ***'


  do step=1, nsteps_per_rank ! nstepsTo25
      pf%state%pfblock = step
      call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., 1, step=(step-1)*pf%comm%nproc + pf%rank)
      
      do m=1, pf%levels(pf%nlevels)%nnodes
         call pf%levels(pf%nlevels)%Q(m)%pack(targetState(step,m,:),1)
         write(npyfname, "(A,'s',i0.2,'m',i0.2,'.npy')") "y", step, m
         call save_double(npyfname, shape(targetState(step,m,:)), real(targetState(step,m,:),8))
      end do
      
      if( step < nsteps_per_rank) then
        call pf%levels(pf%nlevels)%qend%pack(pf%levels(pf%nlevels)%send, 1)    !<  Pack away your last solution
        call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, pf%comm%nproc-1)
        call pf%levels(pf%nlevels)%q0%unpack(pf%levels(pf%nlevels)%send, 1)    !<  Everyone resets their q0
      end if
  end do

  call set_ydesired(pf%levels(pf%nlevels)%ulevel%sweeper, targetState)
  do l = pf%nlevels-1,1,-1
    call restrict_ydesired(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
  end do
  
  call exit(0)

  if(pf%rank == 0) &
     open(unit=105, file = logfilename , &
         status = 'unknown',  action = 'write')

  if(pf%rank == 0) write(105,*) "iter   L2_grad   objective   stepsize"

  allocate(gradient(nvars(pf%nlevels)))
  allocate(prevGrad(nvars(pf%nlevels)))
  allocate(searchDir(nvars(pf%nlevels)))
  allocate(prevSearchDir(nvars(pf%nlevels)))
  allocate(savedStatesFlat(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, product(pf%levels(pf%nlevels)%lev_shape)))
  allocate(savedAdjointsFlat(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, product(pf%levels(pf%nlevels)%lev_shape)))
  
  ctrl = 0.0_pfdp
  gradient = 0
  prevGrad = 0
  prevGlobDirXGrad = 0.0_pfdp
  globDirXGrad = 0.0_pfdp
  prevStepSize = 0.05_pfdp
  prevSearchDir = 0
  savedStatesFlat = 0.0_pfdp
  savedAdjointsFlat = 0.0_pfdp

  itersState = 0
  itersAdjoint = 0
  skippedState = 0

  call date_and_time(date=date, time=time, values=time_start)
  if (pf%rank .eq. 0) &
    print *, 'start optimization on ', date, ', ',  time

  do k = 1, max_opt_iter

     call dump_control(pf%levels(pf%nlevels)%ulevel%sweeper, pf, 'uk', .false., .true.)

     if(pf%rank == 0) print *, '===============Optimization ITERATION================', k

!      call mpi_barrier(pf%comm%comm, ierror)
!      if ( k .eq. 1 ) then
        call evaluate_objective(pf, q1, dt, nsteps_per_rank, .true., alpha, objective, L2NormCtrlSq, savedStatesFlat, ctrl)
        itersState = itersState + pf%state%itcnt - pf%state%skippedy
        skippedState = skippedState + pf%state%skippedy

     call mpi_allreduce(objective, globObj, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
     if(pf%rank == 0) print *, k, 'objective (L2) = ', globObj
     exit

!      call mpi_barrier(pf%comm%comm, ierror)
!      if ( k .eq. 1 ) then   ! in later iterations, this is done in linesearch
       call evaluate_gradient(pf, q1, dt, nsteps_per_rank, .true., gradient, LinftyNormGrad, L2NormGradSq, savedStatesFlat, &
                              savedAdjointsFlat)
       itersAdjoint = itersAdjoint + pf%state%itcnt

     ! done optimizing?
!      call mpi_allreduce(L2NormGradSq, globL2NormGradSq, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
      ! here gradient is not time dependent, every rank has correct norm already
     globL2NormGradSq = L2NormGradSq
     globLinftyNormGrad = LinftyNormGrad
!      call mpi_allreduce(LinftyNormGrad, globLinftyNormGrad, 1, MPI_REAL8, MPI_MAX, pf%comm%comm, ierror)
     if(pf%rank == 0) print *, k, 'gradient (L2, Linf) = ', sqrt(globL2NormGradSq), globLinftyNormGrad
     if(pf%rank == 0) write(105,*) k, sqrt(globL2NormGradSq), globObj, prevStepSize
     if (sqrt(globL2NormGradSq) < tol_grad) then
       if(pf%rank == 0) print *, 'optimality condition satisfied (gradient norm small enough), stopping'
       exit
     end if
     if (globObj < tol_obj) then
       if(pf%rank == 0) print *, 'optimality condition satisfied (objective function small enough), stopping'
       exit
     end if

     if ( k .eq. 1 ) then
        beta = 0.0
     else
!         !PR:
!           globNum   = compute_scalar_prod(gradient, gradient-prevGrad, pf%levels(pf%nlevels)%nodes, dt)
!           globDenom = compute_scalar_prod(prevGrad, prevGrad, pf%levels(pf%nlevels)%nodes, dt)
!         !DY:
          globNum   = compute_scalar_prod(gradient, gradient)
          globDenom = compute_scalar_prod(gradient-prevGrad, prevSearchDir)
!         !FR:
!          globNum   = compute_scalar_prod(gradient, gradient, pf%levels(pf%nlevels)%nodes, dt)
!          globDenom = compute_scalar_prod(prevGrad, prevGrad, pf%levels(pf%nlevels)%nodes, dt)
!           call mpi_allreduce(num,   globNum,   1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!           call mpi_allreduce(denom, globDenom, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
          beta = globNum/globDenom
          !print *,  pf%rank, k, 'beta = ', beta, 'num', globNum, 'denom', globDenom
     end if

     searchDir = -gradient ! + beta*prevSearchDir

     ! loop for step size control
     globDirXGrad = compute_scalar_prod(searchDir, gradient)
!      directionTimesGradient = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
!      call mpi_allreduce(directionTimesGradient, globDirXGrad, 1, MPI_REAL8, MPI_SUM, &
!                          pf%comm%comm, ierror)


     !print *, k, pf%rank, pf%nlevels, 'local gradient times searchDir = ', directionTimesGradient
     if(pf%rank == 0) print *, k, 'gradient times searchDir = ', globDirXGrad


!       stepSize = max_step_size !1.0_pfdp
!       stepSize = 0.0_pfdp

!      stepSize = max(prevStepSize * prevGlobDirXGrad / globDirXGrad, 1.0)
!      stepSize = min(2.0_pfdp*prevStepSize, 1.0_pfdp)

     ! update for next iteration, as these get re-assigned during linesearch
     prevGlobDirXGrad = globDirXGrad
     prevSearchDir = searchDir
     prevGrad = gradient

     stepTooSmall = .false.

     stepSize = prevStepSize
     if(pf%rank == 0) print *, k, 'step size = ', stepSize
     ctrl = ctrl + stepsize*searchDir
     !call update_control(pf%levels(pf%state%finest_level)%ulevel%sweeper, searchDir, stepSize)

    !  call armijo_step(pf, q1, dt, nsteps_per_rank, itersState, itersAdjoint, skippedState, predict, searchDir, gradient, &
!                       savedStatesFlat, savedAdjointsFlat, alpha, &
!                       globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall, &
!                       L2errorState, LinfErrorState, L2exactState, LinfExactState, targetState)
!      call strong_wolfe_step(pf, q1, dt, nsteps_per_rank, itersState, itersAdjoint, predict, searchDir, gradient, &
!                       savedStatesFlat, savedAdjointsFlat, alpha, &
!                       globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall, &
!                       L2errorState, LinfErrorState, L2exactState, LinfExactState, targetState)

     if (stepTooSmall) exit

     ! update for next iteration
     !prevStepSize = 10*stepSize
     prevStepSize = stepSize

   end do

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

  call dump_control(pf%levels(pf%nlevels)%ulevel%sweeper, pf, 'ufinal', .true., .true.)
  if(pf%rank == pf%comm%nproc-1) then ! get final state
      npyOutput = get_array3d_oc(pf%levels(pf%nlevels)%qend, 1)
      call save_double('yfinal.npy', shape(npyOutput), real(npyOutput,8))
      npyOutput = npyOutput - targetState
      call save_double('yfinaldiff.npy', shape(npyOutput), real(npyOutput,8))
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
