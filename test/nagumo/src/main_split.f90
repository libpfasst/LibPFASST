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


  type(pf_ndarray_oc_t) :: q1, qend
  integer        :: ierror, iprovided, l, m, k, i, p
  character(len = 64) :: fout
  character(len =128) :: logfilename
  character(len = 64) :: fname
  character(len = 64) :: probin_fname
  character(len = 64) :: shell_cmd
  integer        :: iout,nout

  real(pfdp), pointer :: gradient(:,:,:), prevGrad(:,:,:), searchDir(:,:,:), prevSearchDir(:,:,:), &
                         savedStates(:,:,:), savedAdjoint(:,:,:), solAt25(:)
  real(pfdp)          :: LinftyNormGrad, LinftyNormU, objective, objectiveNew, L2NormUSq, L2NormGradSq, &
                         dx, stepSize, prevStepSize, directionTimesGradient, beta, &
                         globObj, globObjNew, globDirXGrad, prevGlobDirXGrad, globL2NormGradSq, &
                         num, denom, globNum, globDenom, globLinftyNormGrad, &
                         L2errorCtrl, LinfErrorCtrl, globL2errorCtrl, globLinfErrorCtrl, &
                         L2exactCtrl, LinfEXactCtrl, globL2exactCtrl, globLinfExactCtrl, &
                         abs_res_tol_save, rel_res_tol_save, gamma_save, &
                         initialGradNorm
  logical             :: stepTooSmall, predict, retry
  integer             :: itersState, itersAdjoint, root, sumItersState, sumItersAdjoint
  character(8)   :: date
  character(10)  :: time
  integer        :: time_start(8), time_end(8)
  real(pfdp)     :: time_start_sec, time_end_sec
  real(pfdp), pointer :: tmp(:,:), tmpnodes(:), tmprmat(:,:)
  integer ::  stat(MPI_STATUS_SIZE)
  integer ::  nsteps_per_rank, nstepsTo25, step

  
  
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
  
  pf%state%nsteps = nsteps
  
  pf%save_timings = 0
       
  do l = 1, pf%nlevels
!        pf%levels(l)%nsweeps = nsweeps(l)
!        pf%levels(l)%nsweeps_pred = nsweeps_pred(l)

       pf%levels(l)%nnodes = nnodes(l)

       !  Allocate the user specific level object
       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(pf_ndarray_oc_factory_t::pf%levels(l)%ulevel%factory)

       !  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)
       
       call pf_level_set_size(pf,l,[nvars(l)])

       call setup(pf%levels(l)%ulevel%sweeper, nvars(l), nnodes(l), nsteps_per_rank, comm%nproc)

       !  Allocate the shape array for level (here just one dimension)
!        allocate(pf%levels(l)%shape(1))
!        pf%levels(l)%shape(1) = nvars(l)
       !  Set the size of the send/receive buffer
!        pf%levels(l)%mpibuflen  = nvars(l)
  end do

  call pf_pfasst_setup(pf)

  pf%state%finest_level = pf%nlevels


  pf%debug = .false.
  
  do l = 1, pf%nlevels
    call setrank(pf%levels(l)%ulevel%sweeper, pf%rank)
  end do
  
  !
!   ! run
  !

!   call pf_add_hook(pf, 3, PF_PRE_PREDICTOR, echo_error_hook)
!   call pf_add_hook(pf, 3, PF_POST_PREDICTOR, echo_error_hook)
!   call pf_add_hook(pf,3,PF_POST_ITERATION,echo_error_hook)
 call pf_add_hook(pf,pf%nlevels,PF_POST_BLOCK,echo_error_hook) !STEP
!   call pf_add_hook(pf,-1,PF_POST_SWEEP,echo_residual_hook)
!   call pf_add_hook(pf,-1,PF_POST_ITERATION,echo_residual_hook)

  if(do_imex /= 1 .and. lagging .eqv. .false.) &
    write(output, "(A,'_nl')") trim(output)
    
  write(output, "(A,'_gamma',F0.2)") trim(output), gamma
  if(pf%rank == 0) print *, 'output directory for npy: ', output
  pf%outdir = output
  
  if (len_trim(output) > 0) then
    !  Make directory for Data if it does not exist
!     call ndarray_mkdir(output, len_trim(output))
!      call pf_add_hook(pf, pf%nlevels, PF_POST_STEP, ndarray_oc_dump_all_hook)
!      call pf_add_hook(pf, pf%nlevels, PF_POST_STEP, ndarray_oc_dump_hook)
  end if

  if (nsteps < comm%nproc) then
     nsteps = comm%nproc
  end if

  nsteps_per_rank = nsteps / pf%comm%nproc
  if(nsteps_per_rank*pf%comm%nproc /= nsteps) stop "ERROR: nsteps not an integer multiple of nproc"
  if(pf%rank == 0) print *, "computing ", nsteps_per_rank, "steps per CPU"
  
  call initialize_results(pf)

  !call system('if [ ! -e ./Dat ]; then mkdir Dat; fi')
  !shell_cmd = 'if [ ! -e ./Dat/'//trim(fbase)//' ]; then mkdir Dat/'//trim(fbase)//'; fi'
  !call system(shell_cmd)
  ! open output files
  write (fname, "(A,I0.2,A3,I0.3,A6,I0.3,A6,I0.3)") 'Niter',pf%niters,'_Nx',nvars(pf%nlevels),'_Nstep',nsteps,'_Nproc',comm%nproc
  foutbase = 'Dat/'//trim(fbase)//'/'//trim(fname)
!  print *,'foutbase=',foutbase
  
!   logfilename = trim(logfile)
!   if (warmstart .eq. 1) then
!     predict = .false.
!   else
!     predict = .true.
!   endif
  if (warmstart .eq. 1) then
    predict = .false.
    logfile = trim(logfile)//'_warm'
    if(restart_interval < max_opt_iter) then
      write (logfilename, "(A,'_rest',I0.4)") trim(logfile), restart_interval
      logfile = logfilename
    end if
    if(compress) then
      write (logfilename, "(A,'_compress',I0.12)") trim(logfile), N_digits
      logfile = logfilename
      write (logfilename, "(A,'_storage',I0.12)") trim(logfile), storage_level
      logfile = logfilename
      if(inexact_adjoint) then
logfile = trim(logfile)//'_inexAdj'
end if
end if

else
predict = .true.
logfile = trim(logfile)//'_cold'
endif    
if(adapt_res_tol) then
logfile = trim(logfile)//'_adapt'
end if

write(logfilename, "(A,'_tol',i0.3,'_optiter',i0.4,'_Nstep',i0.3,'_Nproc',i0.3,'.log')") trim(logfile), test_no, max_opt_iter, nsteps, comm%nproc


!  open(unit=101, file = foutbase, status = 'unknown', action = 'write')

  !  Output the run parameters
  if (pf%rank == 0) then
     call pf_print_options(pf, 6,.TRUE.)
!     fout = 'Dat/'//trim(fbase)//'/'//trim(fname)//'_params.m'
!     print*, fout
!     open(unit=103, file = fout, status = 'unknown', action = 'write')
     do iout=2,2 !1,2
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
        write(nout,*) 'nsteps_per_rank=',nsteps_per_rank
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
        write(nout,*) 'gamma=',gamma
        write(nout,*) 'logfile=',logfilename
        if (do_imex .eq. 1) then
          write(nout,*) 'use IMEX sweeper'
        else  
          if( lagging .eqv. .false.) then
            write(nout,*) 'use MISDC sweeper (nonlinear)'            
          else
            write(nout,*) 'use MISDC sweeper (w/ lagging)'
          endif
        endif
        if (warmstart .eq. 1) then
          write(nout,*) 'warm start with previous state solution'
        else
          write(nout,*) 'cold start with spread initial value'
        endif
        if(adapt_res_tol) then
           write(nout, *) 'adapt residual tolerance, (grad_decrease, factor)=', res_tol_graddec, res_tol_factor
        end if
        write(nout,*) 'regularization: alpha=', alpha
        write(nout,*) 'maximum number of optimization iterations=', max_opt_iter
        write(nout,*) 'optimization tolerances (gradient, objective)=', tol_grad, tol_obj
     end do
  end if

  !alpha = 1e-6 !0.0_pfdp set in probin

  call ndarray_oc_build(q1, pf%levels(pf%nlevels)%lev_shape)
  do l = 1, pf%nlevels
     call initialize_oc(pf%levels(l)%ulevel%sweeper, alpha)
  end do

  ! for Nagumo model: solve state equation with zero right hand side to determine natural
  ! solution used in objective function
  ! solve up to high accuracy
  gamma_save = gamma
!   gamma = 1.0_pfdp
  abs_res_tol_save = pf%abs_res_tol
  rel_res_tol_save = pf%rel_res_tol
  if(pf%rank == 0) &
    print *, abs_res_tol_save, rel_res_tol_save, gamma_save
  
  
!   print *, n_digits, 1.23456789, floor(1.23456789*N_digits)/dble(N_digits)
  
  pf%abs_res_tol = 1e-11
  pf%rel_res_tol = 1e-11
  
  pf%q0_style = 0
  
  if (pf%rank .eq. 0) then
    call initial(q1, 0.0_pfdp, dt) 
  else
    q1%yflatarray = 0.0_pfdp
  end if
  call pf%levels(pf%nlevels)%q0%copy(q1, 1)
  if (pf%rank .eq. 0)  print *, ' **** solve state with zero control ***'

  ! solution at t=2.5 has to be send to all later ranks
  allocate(solAt25(nvars(pf%nlevels)))   
  solAt25 = 0.0_pfdp

  ! only compute up to t=2.5
  nstepsTo25=floor(2.5/dt)
!   if (pf%rank .eq. 0)  print *, nstepsTo25
  do step=1, nsteps_per_rank ! nstepsTo25
      pf%state%pfblock = step
      call pf_pfasst_block_oc(pf, dt, step*pf%comm%nproc, .true., 1, step=(step-1)*pf%comm%nproc + pf%rank)     
      call fill_ydesired_nagumo(pf%levels(pf%nlevels)%ulevel%sweeper, pf, step)
      if ( abs(((step-1)*pf%comm%nproc + pf%rank+1)*dt-2.5) < 1e-6) then 
        print *, pf%rank, ((step-1)*pf%comm%nproc + pf%rank+1)*dt
        call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%pack(solAt25, 1)
        
!         print *, "----state at t=2.5----"
!         call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%eprint(1)
!         print *, "---------------------------"              
      end if
      
!     call pf_pfasst_block_oc(pf, dt, nsteps, .true., flags=1)
      if( step < nsteps_per_rank) then
        call pf%levels(pf%nlevels)%qend%pack(pf%levels(pf%nlevels)%send, 1)    !<  Pack away your last solution
        call pf_broadcast(pf, pf%levels(pf%nlevels)%send, pf%levels(pf%nlevels)%mpibuflen, pf%comm%nproc-1)
        call pf%levels(pf%nlevels)%q0%unpack(pf%levels(pf%nlevels)%send, 1)    !<  Everyone resets their q0
      end if
  end do
  
  ! determine which cpu has 2.5 as final time
  ! smarter way?
  do step=1,nsteps_per_rank
    do m=0, pf%comm%nproc-1
      if (abs(((step-1)*pf%comm%nproc + m + 1)*dt -2.5) < 1e-6 ) then
        root = m
      end if
    end do
  end do
!   if ( pf%rank==pf%comm%nproc-1) then
!     call pf%levels(pf%nlevels)%Q(pf%levels(pf%nlevels)%nnodes)%pack(solAt25, 1)
!   end if 
!   root = pf%comm%nproc-1
!   print *, pf%rank, root

  call mpi_bcast(solAt25, nvars(pf%nlevels), MPI_REAL8, root, pf%comm%comm, ierror)
  do m = 1, pf%levels(pf%nlevels)%nnodes
    call pf%levels(pf%nlevels)%Q(m)%unpack(solAt25, 1)
  end do
  ! now every rank has values for the objective in Q(m) -> fill remaining y_desired with this
  do step=1, nsteps_per_rank
    if( ((step-1)*pf%comm%nproc + pf%rank+1)*dt > 2.5) then
      call fill_ydesired_nagumo(pf%levels(pf%nlevels)%ulevel%sweeper, pf, step)
    end if
  end do
  call initialize_control(pf%levels(pf%nlevels)%ulevel%sweeper, solAt25, 0.0_pfdp, dt, pf%levels(pf%nlevels)%nodes)
  do l = pf%nlevels-1,1,-1
    call restrict_control(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
    call restrict_ydesired(pf%levels(l)%ulevel%sweeper, pf%levels(l+1)%ulevel%sweeper)
  end do

  pf%abs_res_tol = abs_res_tol_save
  pf%rel_res_tol = rel_res_tol_save     
  gamma = gamma_save
  
  if(pf%rank == 0) &
    print *, pf%abs_res_tol, pf%rel_res_tol, gamma
  
 ! if(pf%rank == 0) &
 !    open(unit=105, file = logfilename , & 
 !        status = 'unknown',  action = 'write')
 ! if(pf%rank == 0) write(105,*) 'iter ', 'L2_grad ', 'objective ', 'stepsize '
         
  allocate(gradient(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(prevGrad(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(searchDir(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(prevSearchDir(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(savedStates(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  allocate(savedAdjoint(nsteps_per_rank, pf%levels(pf%nlevels)%nnodes, nvars(pf%nlevels)))
  gradient = 0
  prevGrad = 0
  prevGlobDirXGrad = 0.0_pfdp
  globDirXGrad = 0.0_pfdp
  prevStepSize = 0.1_pfdp
  prevSearchDir = 0
  savedAdjoint = 0.0_pfdp
  

  itersState = 0
  itersAdjoint = 0
  
  retry = .true.
  
  call date_and_time(date=date, time=time, values=time_start)
  if (pf%rank .eq. 0) &
    print *, 'start optimization on ', date, ', ',  time

  do k = 1, max_opt_iter

     if(pf%rank == 0) print *, '===============Optimization ITERATION================', k
     
!      call mpi_barrier(pf%comm%comm, ierror)
     if ( k .eq. 1 ) then
        call evaluate_objective(pf, q1, dt, nsteps_per_rank, .true., alpha, objective, L2NormUSq, savedStates)
        itersState = itersState + pf%state%itcnt

!         if(pf%rank == pf%comm%nproc-1) then
!           print *, "----state at final time----"
!           call pf%levels(pf%nlevels)%qend%eprint(1)
!           print *, "---------------------------"
!           call dump_stuff(pf, savedStates, "y")
!         end if
!         exit
     else
        objective = objectiveNew  !can we return globObjNew from linesearch, avoids communication; is objective needed somewhere?
     end if 
     call mpi_allreduce(objective, globObj, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
     if(pf%rank == 0) print *, k, 'objective (L2) = ', globObj 
     !exit 
     
!      call mpi_barrier(pf%comm%comm, ierror)
     if ( k .eq. 1 ) then   ! in later iterations, this is done in linesearch
       call evaluate_gradient(pf, q1, dt, nsteps_per_rank, .true., gradient, LinftyNormGrad, L2NormGradSq, & 
                              savedStates, savedAdjoint)
            
       itersAdjoint = itersAdjoint + pf%state%itcnt
       
!        if(pf%rank == 0) then
!           print *, "----adjoint at t=0----"
!           call pf%levels(pf%nlevels)%q0%eprint(2)
!           print *, "---------------------------"
!           call dump_stuff(pf, savedAdjoint, "p")
!         end if

     end if
     !print *, k, pf%rank, pf%nlevels, 'objective (L2) = ', objective
!      print *, k, pf%rank, pf%nlevels, 'gradient (Linfty, L2) = ', LinftyNormGrad, sqrt(L2NormGradSq)
!     exit
     ! done optimizing?
     call mpi_allreduce(L2NormGradSq, globL2NormGradSq, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
     call mpi_allreduce(LinftyNormGrad, globLinftyNormGrad, 1, MPI_REAL8, MPI_MAX, pf%comm%comm, ierror)
     if(pf%rank == 0) print *, k, 'gradient (L2, Linf) = ', sqrt(globL2NormGradSq), globLinftyNormGrad
!     if(pf%rank == 0) write(105,*) k, sqrt(globL2NormGradSq), globObj, prevStepSize
     if(k == 1) then 
        initialGradNorm = sqrt(globL2NormGradSq)
     else
        if(adapt_res_tol .and. sqrt(globL2NormGradSq) < res_tol_graddec*initialGradNorm) then
           initialGradNorm = sqrt(globL2NormGradSq)
           pf%abs_res_tol = pf%abs_res_tol * res_tol_factor
           pf%rel_res_tol = pf%rel_res_tol * res_tol_factor
           if(pf%rank == 0) print *, 'adapted residual tolerances: ', pf%abs_res_tol, pf%rel_res_tol
           predict = .true.
        end if
     end if
     if (sqrt(globL2NormGradSq) < tol_grad) then
       if(pf%rank == 0) print *, 'optimality condition satisfied (gradient norm small enough), stopping'
       !call write_control_work1(pf%levels(pf%nlevels)%ctx, k, "u_sdc_split_final")
       exit
     end if
     if (globObj < tol_obj) then
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
          call mpi_allreduce(num,   globNum,   1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
          call mpi_allreduce(denom, globDenom, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
          beta = globNum/globDenom
          !print *,  pf%rank, k, 'beta = ', beta, 'num', globNum, 'denom', globDenom
     end if

     searchDir = -gradient + beta*prevSearchDir

     ! loop for step size control
     directionTimesGradient = compute_scalar_prod(searchDir, gradient, pf%levels(pf%nlevels)%nodes, dt)
     call mpi_allreduce(directionTimesGradient, globDirXGrad, 1, MPI_REAL8, MPI_SUM, &
                         pf%comm%comm, ierror)
     
     !print *, k, pf%rank, pf%nlevels, 'local gradient times searchDir = ', directionTimesGradient
     if(pf%rank == 0) print *, k, 'gradient times searchDir = ', globDirXGrad


!       stepSize = 1.0_pfdp
!      stepSize = max(prevStepSize * prevGlobDirXGrad / globDirXGrad, 1.0)
     stepSize = 10*prevStepSize
     
     ! update for next iteration, as these get re-assigned during linesearch
     prevGlobDirXGrad = globDirXGrad
     prevSearchDir = searchDir
     prevGrad = gradient

! if (warmstart .eq. 0 ) then
!   ! set Q, F to zero
!   do l = 1, pf%nlevels-1
!     do m = 1, pf%levels(l)%nnodes
!       call pf%levels(l)%encap%setval(pf%levels(l)%Q(m), 0.0_pfdp, 0) 
!       call pf%levels(l)%encap%setval(pf%levels(l)%pQ(m), 0.0_pfdp, 0) ! only on coarser levels, also tauQ
!       if(m < pf%levels(l)%nnodes) then
!         call pf%levels(l)%encap%setval(pf%levels(l)%tauQ(m), 0.0_pfdp, 0)
!         call pf%levels(l)%encap%setval(pf%levels(l)%tau(m), 0.0_pfdp, 0)
!         call pf%levels(l)%encap%setval(pf%levels(l)%I(m), 0.0_pfdp, 0) 
!         call pf%levels(l)%encap%setval(pf%levels(l)%S(m), 0.0_pfdp, 0) 
!       end if
!       do p = 1,size(pf%levels(l)%F(1,:))
!         call pf%levels(l)%encap%setval(pf%levels(l)%F(m,p), 0.0_pfdp, 0)
!         !call pf%levels(l)%encap%setval(pf%levels(l)%pF(m,p), 0.0_pfdp, 0)
!       end do
!     end do
!   end do
!   
!   do m = 1, pf%levels(pf%nlevels)%nnodes
!     call pf%levels(pf%nlevels)%encap%setval(pf%levels(pf%nlevels)%Q(m), 0.0_pfdp, 0) 
!     do p = 1,size(pf%levels(pf%nlevels)%F(1,:))
!         call pf%levels(pf%nlevels)%encap%setval(pf%levels(pf%nlevels)%F(m,p), 0.0_pfdp, 0)
!     end do
!   end do
! end if

!      if ((.not. predict) .and. mod(k,restart_interval) == 0 ) then ! use cold start every n iterations
!        predict = .true.
!      else if(warmstart == 1) then
!        predict = .false.
!      end if
      
     stepTooSmall = .false.
!      call armijo_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
!                       globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)
!      call wolfe_powell_step(pf, q1, dt, nsteps, itersState, itersAdjoint, predict, searchDir, gradient, savedAdjoint, alpha, &
!                     globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)   
     call strong_wolfe_step(pf, q1, dt, nsteps_per_rank, itersState, itersAdjoint, predict, searchDir, gradient, &
                    savedStates, savedAdjoint, alpha, &
                    globObj, globDirXGrad, objectiveNew, L2NormUSq, LinftyNormGrad, L2NormGradSq, stepSize, stepTooSmall)    
                      
     if (stepTooSmall) then
        if (adapt_res_tol .and. retry) then
           initialGradNorm = sqrt(globL2NormGradSq)
           pf%abs_res_tol = pf%abs_res_tol * res_tol_factor
           pf%rel_res_tol = pf%rel_res_tol * res_tol_factor
           if(pf%rank == 0) print *, 'adapted residual tolerances: ', pf%abs_res_tol, pf%rel_res_tol
           predict = .true.
        else
           exit
        end if
     end if

     ! update for next iteration
     !prevStepSize = 10*stepSize
     prevStepSize = stepSize
          
   end do
   
   call mpi_barrier(pf%comm%comm, ierror)
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
  call dump_exact_control(pf%levels(pf%nlevels)%ulevel%sweeper, pf, solAt25, 0.0_pfdp, dt, pf%levels(pf%nlevels)%nodes, &
                                  L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl)
  call mpi_allreduce(L2errorCtrl, globL2errorCtrl, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
  call mpi_allreduce(LinfErrorCtrl, globLinfErrorCtrl, 1, MPI_REAL8, MPI_MAX, pf%comm%comm, ierror)
  call mpi_allreduce(L2exactCtrl, globL2exactCtrl, 1, MPI_REAL8, MPI_SUM, pf%comm%comm, ierror)
  call mpi_allreduce(LinfExactCtrl, globLinfExactCtrl, 1, MPI_REAL8, MPI_MAX, pf%comm%comm, ierror) 
   
  print *, 'rank:', pf%rank, 'total iterations state, adjoint', itersState, itersAdjoint

   sumItersState = 0
   sumItersAdjoint = 0
   call mpi_allreduce(itersState,   sumItersState,   1, MPI_INTEGER, MPI_SUM, pf%comm%comm, ierror)
   call mpi_allreduce(itersAdjoint,   sumItersAdjoint,   1, MPI_INTEGER, MPI_SUM, pf%comm%comm, ierror)
  
   if( pf%rank == 0 ) then
      print *, 'overall sum iterations state, adjoint', sumItersState, sumItersAdjoint
      print *, 'absolute error in computed control: L2 = ', sqrt(globL2errorCtrl)
      print *, '                                  Linf = ', globLinfErrorCtrl
      print *, 'relative error in computed control: L2 = ', sqrt(globL2errorCtrl)/sqrt(globL2exactCtrl)
      print *, '                                  Linf = ', globLinfErrorCtrl/globLinfExactCtrl   
      print *, 'exact control:                      L2 = ', sqrt(globL2exactCtrl)
      print *, '                                  Linf = ', globLinfExactCtrl

   end if

   
  !
  ! cleanup
  !

!   print *, pf%rank, "cleanup"
  deallocate(gradient)
  deallocate(prevGrad)
  deallocate(searchDir)
  deallocate(prevSearchDir)
  deallocate(savedAdjoint)

!   print *, pf%rank, "destroy q1"
  call ndarray_oc_destroy(q1)
  
!   print *, pf%rank, "destroy levels"
   do l = 1, pf%nlevels
    call destroy(pf%levels(l)%ulevel%sweeper)
  end do

!   print *, pf%rank, "destroy pf"
  call pf_pfasst_destroy(pf)
  
!   print *, pf%rank, "mpi_finalize"
  call mpi_finalize(ierror)
!   call fftw_cleanup()

end program main
