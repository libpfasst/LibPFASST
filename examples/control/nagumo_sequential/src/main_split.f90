program main
  use pfasst
  use pf_mod_mpi 
  use pf_mod_ndarray_oc

  use pf_mod_optimization
  use feval
  use hooks
  use probin
  use solutions
!   use transfer

  implicit none

  type(pf_comm_t)          :: comm
  type(pf_pfasst_t)        :: pf
  class(ad_sweeper_t), pointer   :: ad_sweeper_ptr     !<  pointer to SDC sweeper 


  type(ndarray_oc) :: q1, qend
  integer        :: ierror, iprovided, l, m, k, i, p
  character(len = 64) :: fout
  character(len =128) :: logfilename
  character(len = 64) :: fname
  character(len = 64) :: probin_fname
  character(len = 64) :: shell_cmd
  integer        :: iout,nout

  real(pfdp), pointer :: gradient(:,:,:), prevGrad(:,:,:), searchDir(:,:,:), prevSearchDir(:,:,:), solAt25(:), &
                         savedStates(:,:,:), stateSolution(:,:,:), savedAdjoint(:,:,:)
  real(pfdp)          :: LinftyNormGrad, LinftyNormU, objective, objectiveNew, L2NormUSq, L2NormGradSq, &
                         alpha, dx, stepSize, prevStepSize, directionTimesGradient, beta, &
                         globObj, globObjNew, globDirXGrad, prevGlobDirXGrad, globL2NormGradSq, tolGrad, tolObj, &
                         num, denom, globNum, globDenom, globLinftyNormGrad, &
                         L2errorCtrl, LinfErrorCtrl, globL2errorCtrl, globLinfErrorCtrl, &
                         L2exactCtrl, LinfEXactCtrl, globL2exactCtrl, globLinfExactCtrl, &
                         abs_res_tol_save, rel_res_tol_save
  logical             :: stepTooSmall, predict
  integer             :: itersState, itersAdjoint, root, sumItersState, sumItersAdjoint, step, nstepsTo25
  character(8)   :: date
  character(10)  :: time
  integer        :: time_start(8), time_end(8)
  real(pfdp)     :: time_start_sec, time_end_sec
  
  
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
  
  pf%echo_timings = .false.
       
  do l = 1, pf%nlevels
       pf%levels(l)%nsweeps = nsweeps(l)
       pf%levels(l)%nsweeps_pred = nsweeps_pred(l)

       pf%levels(l)%nnodes = nnodes(l)

       !  Allocate the user specific level object
       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(ndarray_oc_factory::pf%levels(l)%ulevel%factory)

       !  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)
       call setup(pf%levels(l)%ulevel%sweeper, nvars(l), nnodes(l), nsteps)

       !  Allocate the shape array for level (here just one dimension)
       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = nvars(l)
       !  Set the size of the send/receive buffer
       pf%levels(l)%mpibuflen  = nvars(l)
  end do

  call pf_pfasst_setup(pf)

  pf%outdir = output

  !
!   ! run
  !

!   call pf_add_hook(pf, 3, PF_PRE_PREDICTOR, echo_error_hook)
!   call pf_add_hook(pf, 3, PF_POST_PREDICTOR, echo_error_hook)
!   call pf_add_hook(pf,3,PF_POST_ITERATION,echo_error_hook)
  call pf_add_hook(pf,-1,PF_POST_STEP,echo_error_hook)
!   call pf_add_hook(pf,-1,PF_POST_SWEEP,echo_residual_hook)
!   call pf_add_hook(pf,-1,PF_POST_ITERATION,echo_residual_hook)


  if (len_trim(output) > 0) then
     call ndarray_mkdir(output, len_trim(output))
     call pf_add_hook(pf, -1, PF_POST_STEP, ndarray_oc_dump_all_hook)
  end if

  
  if (comm%nproc > 1 ) stop 'this is the sequential version, can only be run on one rank!'

  allocate(stateSolution(nsteps, nnodes(pf%nlevels), nvars(pf%nlevels)))

  !  Make directory for Data if it does not exist

  !if (len_trim(output) > 0) then
  !   call ndarray_mkdir(output, len_trim(output))
  !   call pf_add_hook(pf, -1, PF_POST_SWEEP, ndarray_oc_dump_hook)
  !end if


  call system('if [ ! -e ./Dat ]; then mkdir Dat; fi')
  shell_cmd = 'if [ ! -e ./Dat/'//trim(fbase)//' ]; then mkdir Dat/'//trim(fbase)//'; fi'
  call system(shell_cmd)
  ! open output files
  write (fname, "(A,I0.2,A3,I0.3,A6,I0.3,A6,I0.3)") 'Niter',pf%niters,'_Nx',nvars(pf%nlevels),'_Nstep',nsteps,'_Nproc',comm%nproc
  foutbase = 'Dat/'//trim(fbase)//'/'//trim(fname)
!  print *,'foutbase=',foutbase
  
  logfilename = trim(logfile)
  if (warmstart .eq. 1) then
    predict = .false.
  else
    predict = .true.
  endif

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
        write(nout,*) 'nprob=',nprob
        write(nout,*) 'nsteps=',nsteps
        write(nout,*) 'dt=',dt
        write(nout,*) 'nu=',nu
        write(nout,*) 'v=',v
        write(nout,*) 'kfreq=',kfreq
        write(nout,*) 'do_spec',do_spec
        write(nout,*) 'fbase=',fbase
        if (do_spec .eq. 0) then
           write(nout,*) 'N_Vcycles=',N_Vcycles
           write(nout,*) 'Nrelax',Nrelax
           write(nout,*) 'mg_interp_order=',mg_interp_order
        endif
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

     end do
  end if

  alpha = 1e-6 !0.0_pfdp

  call ndarray_oc_build(q1, pf%levels(pf%nlevels)%shape)
  do l = 1, pf%nlevels
     call initialize_oc(pf%levels(l)%ulevel%sweeper, pf%rank*dt, dt, pf%levels(l)%nodes, nvars(l), alpha)
  end do

  ! for Nagumo model: solve state equation with zero right hand side to determine natural
  ! solution used in objective function
  ! solve up to high accuracy
  abs_res_tol_save = abs_res_tol
  rel_res_tol_save = rel_res_tol
  abs_res_tol = 1e-11
  rel_res_tol = 1e-11
  
  pf%q0_style = 0
 
  if (pf%rank .eq. 0) then
    call initial(q1, pf%rank*dt, pf%rank*dt+dt) !dt should be given by rank !initial_rd
  else
    q1%yflatarray = 0.0_pfdp
  end if
  call pf%levels(pf%nlevels)%q0%copy(q1, 1)
  
  if (pf%rank .eq. 0)  print *, ' **** solve state with zero control until T=2.5 ****'
  nstepsTo25 = floor(2.5/dt)
  do step = 1, nstepsTo25
    call pf_pfasst_block_oc(pf, dt, step, .true., 1, step-1)
    call fill_ydesired_nagumo(pf%levels(pf%nlevels)%ulevel%sweeper, pf, step)
    call pf%levels(pf%nlevels)%q0%copy(pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes),1)
  end do
    
  print *, ' **** use solution at 2.5 to fill remaining ydesired ****'
  allocate(solAt25(nvars(pf%nlevels)))    
  call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%pack(solAt25, 1)
  do step = nstepsTo25+1, nsteps
    call fill_ydesired_nagumo_25(pf%levels(pf%nlevels)%ulevel%sweeper, pf, step, solAt25)
  end do
  
  print *, ' **** initializing control ****'
  call initialize_control(pf%levels(pf%nlevels)%ulevel%sweeper, solAt25, pf%rank*dt, dt, pf%levels(pf%nlevels)%nodes)
  print *, ' **** restricting control and ydesired ****'
  do l = pf%nlevels-1,1,-1
    call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
    call dump_control(pf%levels(l)%ulevel%sweeper, pf, 'u', l)
    call restrict_ydesired(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
  end do

  abs_res_tol = abs_res_tol_save
  rel_res_tol = rel_res_tol_save      

  tolGrad = 1.0e-6
  tolObj  = 1.0e-6


  
  if(pf%rank == 0) &
     open(unit=105, file = logfilename , & 
         status = 'unknown',  action = 'write')

         
  allocate(gradient(nsteps, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(prevGrad(nsteps, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(searchDir(nsteps, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(prevSearchDir(nsteps, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
!   allocate(savedAdjoint(pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(savedStates(nsteps, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  gradient = 0
  prevGrad = 0
  prevGlobDirXGrad = 0.0_pfdp
  globDirXGrad = 0.0_pfdp
  prevStepSize = 0.1_pfdp
  prevSearchDir = 0
  savedStates = 0.0_pfdp
  

  itersState = 0
  itersAdjoint = 0

  
  call date_and_time(date=date, time=time, values=time_start)
  if (pf%rank .eq. 0) &
    print *, 'start optimization on ', date, ', ',  time

  do k = 1, 200
  
     if(pf%rank == 0) print *, '===============Optimization ITERATION================', k
     
!      call mpi_barrier(pf%comm%comm, ierror)
     if ( k .eq. 1 ) then
        call evaluate_objective(pf, q1, dt, nsteps, .true., alpha, objective, L2NormUSq, savedStates)
        itersState = itersState + pf%state%itcnt
!         exit
     else
        objective = objectiveNew  !can we return globObjNew from linesearch, avoids communication; is objective needed somewhere?
     end if 
!      call mpi_allreduce(objective, globObj, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
     if(pf%rank == 0) print *, k, 'objective (L2) = ', objective 
!      exit 
     
!      call mpi_barrier(pf%comm%comm, ierror)
     if ( k .eq. 1 ) then   ! in later iterations, this is done in linesearch
       call evaluate_gradient(pf, q1, dt, nsteps, .true., gradient, LinftyNormGrad, L2NormGradSq, savedStates)
       itersAdjoint = itersAdjoint + pf%state%itcnt
     end if
     !print *, k, pf%rank, pf%nlevels, 'objective (L2) = ', objective
     !print *, k, pf%rank, pf%nlevels, 'gradient (Linfty, L2) = ', LinftyNormGrad, sqrt(L2NormGradSq)
!      exit
     ! done optimizing?
!      call mpi_allreduce(L2NormGradSq, globL2NormGradSq, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!      call mpi_allreduce(LinftyNormGrad, globLinftyNormGrad, 1, MPI_REAL8, MPI_MAX, pf%comm%comm, ierror)
     if(pf%rank == 0) print *, k, 'gradient (L2, Linf) = ', sqrt(L2NormGradSq), LinftyNormGrad
     if(pf%rank == 0) write(105,*) k, sqrt(L2NormGradSq), objective, prevStepSize
     if (sqrt(L2NormGradSq) < tolGrad) then
       if(pf%rank == 0) print *, 'optimality condition satisfied (gradient norm small enough), stopping'
       !call write_control_work1(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
       exit
     end if
     if (objective < tolObj) then
       if(pf%rank == 0) print *, 'optimality condition satisfied (objective function small enough), stopping'
       !call write_control_work1(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
       exit
     end if
     
! exit

     
     if ( k .eq. 1 ) then
        beta = 0.0
     else
!         !PR:
          !num   = compute_scalar_prod(gradient, gradient-prevGrad, pf%levels(pf%nlevels)%nodes, dt)
          !denom = compute_scalar_prod(prevGrad, prevGrad, pf%levels(pf%nlevels)%nodes, dt)
!         !DY:
          num   = compute_scalar_prod(gradient, gradient, pf%levels(pf%nlevels)%nodes, dt) 
          denom = compute_scalar_prod(gradient-prevGrad, prevSearchDir, pf%levels(pf%nlevels)%nodes, dt)
!         !FR:
!         !num   = compute_scalar_prod(gradient, gradient, pf%levels(pf%nlevels)%nodes, dt) 
!         !denom = compute_scalar_prod(prevGrad, prevGrad, pf%levels(pf%nlevels)%nodes, dt)
!           call mpi_allreduce(num,   globNum,   1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!           call mpi_allreduce(denom, globDenom, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
          beta = num/denom
          !print *,  pf%rank, k, 'beta = ', beta, 'num', globNum, 'denom', globDenom
     end if

     searchDir = -gradient + beta*prevSearchDir

     ! loop for step size control
     directionTimesGradient = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
!      call mpi_allreduce(directionTimesGradient, globDirXGrad, 1, MPI_REAL8, MPI_SUM, &
!                          pf%comm%comm, ierror)
     
     !print *, k, pf%rank, pf%nlevels, 'local gradient times searchDir = ', directionTimesGradient
     if(pf%rank == 0) print *, k, 'gradient times searchDir = ', directionTimesGradient


!       stepSize = 1.0_pfdp
!      stepSize = max(prevStepSize * prevGlobDirXGrad / globDirXGrad, 1.0)
     stepSize = 10*prevStepSize
     
     ! update for next iteration, as these get re-assigned during linesearch
     prevGlobDirXGrad = directionTimesGradient
     prevSearchDir = searchDir
     prevGrad = gradient
    
     stepTooSmall = .false.
!      call armijo_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
!                       globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
!      call wolfe_powell_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
!                     globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)   
     call strong_wolfe_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedStates, alpha, &
                    objective, directionTimesGradient, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, &
                    stepTooSmall)    
                      
     if (stepTooSmall) exit

     ! update for next iteration
     !prevStepSize = 10*stepSize
     prevStepSize = stepSize
     
   end do
   
!    call mpi_barrier(pf%comm%comm, ierror)
   call date_and_time(date=date, time=time, values=time_end)
   if (pf%rank .eq. 0) then
     print *, 'end optimization on ', date, ', ',  time
     time_start_sec = time_start(5) * 3600 + time_start(6) * 60 &
           + time_start(7) + 0.001 * time_start(8)
      time_end_sec = time_end(5) * 3600 + time_end(6) * 60 &
           + time_end(7) + 0.001 * time_end(8)
      print *, 'duration [s]: ', time_end_sec-time_start_sec
   end if

  call dump_control(pf%levels(pf%nlevels)%ulevel%sweeper, pf, 'u')
  call dump_exact_control(pf%levels(pf%nlevels)%ulevel%sweeper, pf, solAt25, pf%rank*dt, dt, pf%levels(pf%nlevels)%nodes, &
                                  L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl)
!   call mpi_allreduce(L2errorCtrl, globL2errorCtrl, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!   call mpi_allreduce(LinfErrorCtrl, globLinfErrorCtrl, 1, MPI_REAL8, MPI_MAX, pf%comm%comm, ierror)
!   call mpi_allreduce(L2exactCtrl, globL2exactCtrl, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
!   call mpi_allreduce(LinfExactCtrl, globLinfExactCtrl, 1, MPI_REAL8, MPI_MAX, pf%comm%comm, ierror) 
   
!   print *, 'rank:', pf%rank, 'total iterations state, adjoint', itersState, itersAdjoint

   sumItersState = 0
   sumItersAdjoint = 0
!    call mpi_allreduce(itersState,   sumItersState,   1, MPI_INTEGER, MPI_SUM, pf%comm%comm, ierror)
!    call mpi_allreduce(itersAdjoint,   sumItersAdjoint,   1, MPI_INTEGER, MPI_SUM, pf%comm%comm, ierror)
  
   if( pf%rank == 0 ) then
      print *, 'overall sum iterations state, adjoint', itersState, itersAdjoint
      print *, 'absolute error in computed control: L2 = ', sqrt(L2errorCtrl)
      print *, '                                  Linf = ', LinfErrorCtrl
      print *, 'relative error in computed control: L2 = ', sqrt(L2errorCtrl)/sqrt(L2exactCtrl)
      print *, '                                  Linf = ', LinfErrorCtrl/LinfExactCtrl   
      print *, 'exact control:                      L2 = ', sqrt(L2exactCtrl)
      print *, '                                  Linf = ', LinfExactCtrl
   end if

   
  !
  ! cleanup
  !

  deallocate(gradient)
  deallocate(prevGrad)
  deallocate(searchDir)
  deallocate(prevSearchDir)
  deallocate(savedStates)

  call ndarray_oc_destroy(q1)

  call pf_pfasst_destroy(pf)
  call mpi_finalize(ierror)
!   call fftw_cleanup()

end program main
