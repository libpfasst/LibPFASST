program fpfasst
  use pfasst
  use feval
  use initial
  use hooks
  use encap
  use transfer
  use probin
  use pf_mod_ndarray
  use pf_mod_version
  use pf_mod_fimex
  use pf_mod_mpi, only: MPI_THREAD_MULTIPLE, MPI_COMM_WORLD

  implicit none

  type(pf_pfasst_t)     :: pf
  type(pf_comm_t)       :: tcomm
  type(carray4), target :: q0
  integer               :: nthreads, first

  integer               :: ierror, iprovided, l
  character(len=32)     :: arg
  character(len=256)    :: probin_fname

  type(pf_encap_t),   target :: encaps

  ! initialize mpi
  call mpi_init_thread(mpi_thread_multiple, iprovided, ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."

  ! set nthreads from environment
  nthreads = -1
  if (nthreads < 0) then
     call getenv("OMP_NUM_THREADS", arg)
     if (len_trim(arg) > 0) then
        read (arg, *) nthreads
     else
        nthreads = 1
     end if
  end if

  ! read probin
  if (command_argument_count() == 1) then
     call get_command_argument(1, value=probin_fname)
  else
     probin_fname = "probin.nml"
  end if
  call probin_init(probin_fname)

  ! init pfasst
  call carray4_encap_create(encaps)
  call pf_mpi_create(tcomm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, tcomm, nlevs)

  first = size(nx) - nlevs + 1

  pf%niters       = niters
  pf%qtype        = SDC_GAUSS_LOBATTO
  pf%echo_timings = .true.

  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, pf%nlevels
     pf%levels(l)%nvars  = nvars(first-1+l)
     pf%levels(l)%nnodes = nnodes(first-1+l)

     allocate(pf%levels(l)%shape(4))
     pf%levels(l)%shape = [ nx(first-1+l), nx(first-1+l), nx(first-1+l), 3 ]
     call feval_create(nx(first-1+l), 1.0d0, nu, nthreads, pf%levels(l)%levelctx)

     pf%levels(l)%encap       => encaps
     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     call pf_fimex_create(pf%levels(l)%sweeper, eval_f1, eval_f2, comp_f2, eval_force)
  end do

  call pf_mpi_setup(tcomm, pf)
  call pf_pfasst_setup(pf)

  pf%outdir = output

  ! load initial condition, set hooks
  call carray4_create(q0, pf%levels(nlevs)%shape)
  !call load(q0, input)
  call random_full(q0)

  call dump(output, "initial.npy", q0)
  call pf_add_hook(pf, -1, PF_POST_SWEEP, project_hook)
  ! call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error_hook)
  call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_energy_hook)

  ! run
  if (pf%rank == 0) then
     print *, 'ns run'
     print *, '------'
     print *, 'nx:       ', nx(first:)
     print *, 'nlevs:    ', nlevs
     print *, 'nnodes:   ', nnodes(first:)
     print *, 'nthreads: ', nthreads
     print *, 'nsteps:   ', nsteps
     print *, 'niters:   ', niters
     print *, 'nprocs:   ', pf%comm%nproc
     print *, 'nu:       ', nu
     print *, 'dt:       ', dt
     print *, 'input:    ', trim(input)
     print *, 'output:   ', trim(pf%outdir)
     print *, 'version:  ', pf_git_version
     print *, ''
  end if

  if (len_trim(pf%outdir) > 0) then
     call ndarray_mkdir(pf%outdir, len_trim(pf%outdir))
     call pf_add_hook(pf, -1, PF_POST_SWEEP, dump_hook)
  end if

  ! call pf_logger_attach(pf)
  call mpi_barrier(pf%comm%comm, ierror)
  call pf_pfasst_run(pf, c_loc(q0), dt, 0.0d0, nsteps=nsteps)

  ! done
  call mpi_barrier(pf%comm%comm, ierror)

  do l = 1, nlevs
     call feval_destroy(pf%levels(l)%levelctx)
  end do

  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(tcomm)
  call feval_finalize()

  call mpi_finalize(ierror)

end program fpfasst
