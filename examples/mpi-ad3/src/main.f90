program main
  use pfasst
  use pf_mod_mpi, only: MPI_THREAD_FUNNELED, MPI_COMM_WORLD
  use pf_mod_ndarray

  use feval
  use hooks
  use probin
  use solutions
  use transfer

  implicit none

  type(pf_comm_t)          :: comm
  type(pf_pfasst_t)        :: pf
  type(pf_encap_t), target :: encap

  type(ndarray),    target :: q1
  integer        :: ierror, iprovided, l
  character(len = 64) :: fout
  character(len = 64) :: fname
  character(256) :: probin_fname
  integer        :: iout,nout


  !
  ! read options
  !

  if (command_argument_count() == 1) then
     call get_command_argument(1, value=probin_fname)
  else
     probin_fname = "probin.nml"
  end if
  call probin_init(probin_fname)


  !
  ! initialize mpi
  !
  call mpi_init_thread(mpi_thread_funneled, iprovided, ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst
  !

  call ndarray_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, comm,  fname=pfasst_nml)

  pf%niters = niters
!  pf%qtype  = SDC_GAUSS_LOBATTO ! + SDC_PROPER_NODES
  pf%qtype  = SDC_UNIFORM + SDC_NO_LEFT
  pf%window = wtype

  pf%abs_res_tol = abs_tol
  pf%rel_res_tol = rel_tol

  pf%echo_timings = .true.

  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     allocate(pf%levels(l)%shape(ndim))
     pf%levels(l)%shape(ndim) = nvars(l)

     pf%levels(l)%nvars  = product(pf%levels(l)%shape)
     pf%levels(l)%nnodes = nnodes(l)
     pf%levels(l)%Finterp = Finterp

     pf%levels(l)%encap       => encap
     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict

     select case(ndim)
     case(1)
        call create_work1(pf%levels(l)%levelctx, pf%levels(l)%shape(1))
        call pf_imex_create(pf%levels(l)%sweeper, f1eval1, f2eval1, f2comp1)
     case(2)
        call create_work2(pf%levels(l)%levelctx, pf%levels(l)%shape(1))
        call pf_imex_create(pf%levels(l)%sweeper, f1eval2, f2eval2, f2comp2)
     case(3)
        call create_work3(pf%levels(l)%levelctx, pf%levels(l)%shape(1))
        call pf_imex_create(pf%levels(l)%sweeper, f1eval3, f2eval3, f2comp3)
     end select
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  pf%outdir = output

  !
  ! run
  !

!  call pf_add_hook(pf, nlevs, PF_POST_SWEEP, echo_error_hook)
  call pf_add_hook(pf,-1,PF_POST_ITERATION,echo_error_hook)

  call ndarray_create_simple(q1, pf%levels(nlevs)%shape)
  call initial(q1)

  if (len_trim(output) > 0) then
     call ndarray_mkdir(output, len_trim(output))
     ! call pf_add_hook(pf, nlevs, PF_POST_SWEEP, ndarray_dump_hook)
!     call pf_add_hook(pf, -1, PF_POST_SWEEP, ndarray_dump_hook)
  end if

  if (nsteps < comm%nproc) then
     nsteps = comm%nproc
  end if


  !  Output the run parameters
  if (pf%rank == 0) then
     call pf_print_options(pf, 6,.TRUE.)
     fout = 'Dat/'//trim(fbase)//'_'//trim(fname)//'_params.m'
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
        write(nout,*) 'do_spec',do_spec
        if (do_spec .eq. 0) then
           write(nout,*) 'N_Vcycles=',N_Vcycles
           write(nout,*) 'mg_interp_order=',mg_interp_order
        endif

     end do
  end if


  call mpi_barrier(mpi_comm_world, ierror)
  call pf_pfasst_run(pf, c_loc(q1), dt, 0.0_pfdp, nsteps=nsteps)

  !
  ! cleanup
  !

  deallocate(q1%flatarray)

  ! do l = 1, nlevs
  !    if (ndim == 1) then
  !       call destroy_work1(pf%levels(l)%levelctx)
  !    else
  !       ! call destroy_work2(pf%levels(l)%levelctx)
  !    end if
  ! end do

  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(ierror)
  call fftw_cleanup()

end program main
