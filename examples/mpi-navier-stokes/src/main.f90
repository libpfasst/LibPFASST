program fpfasst
  use pfasst
  use feval
  use hooks
  use encap
  use transfer
  use pf_mod_mpi, only: MPI_COMM_WORLD

  implicit none

  type(pf_pfasst_t)     :: pf
  type(pf_comm_t)       :: tcomm
  type(carray4), target :: q0
  integer               :: nlevs, nthreads, nsteps, first
  integer               :: nx(3), nvars(3), nnodes(3)
  integer               :: ierror, l
  double precision      :: dt
  character(len=32)     :: arg

  type(pf_sweeper_t), target :: sweeper
  type(pf_encap_t),   target :: encaps

  ! initialize mpi
  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."

  nthreads = -1
  nsteps   = 32
  nlevs    = 3

  if (nthreads < 0) then
     call getenv("OMP_NUM_THREADS", arg)
     if (len_trim(arg) > 0) then
        read (arg, *) nthreads
     else
        nthreads = 1
     end if
  end if


  ! initialize pfasst
  nx     = [ 32, 64, 128 ]
  nvars  = 2 * 3 * nx**3
  nnodes = [ 2, 3, 5 ]
  dt     = 0.001d0
  first  = size(nx) - nlevs + 1

  call carray4_encap_create(encaps)
  call pf_mpi_create(tcomm, MPI_COMM_WORLD)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)
  call pf_pfasst_create(pf, tcomm, nlevs)

  pf%niters = 5
  pf%qtype  = 1

  pf%echo_timings = .true.
  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     pf%levels(l)%nvars  = nvars(first-1+l)
     pf%levels(l)%nnodes = nnodes(first-1+l)

     allocate(pf%levels(l)%shape(4))
     pf%levels(l)%shape = [ nx(first-1+l), nx(first-1+l), nx(first-1+l), 3 ]
     call feval_create(nx(first-1+l), 1.0d0, 2.0d-3, nthreads, pf%levels(l)%levelctx)

     pf%levels(l)%encap       => encaps
     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%sweeper     => sweeper
  end do

  call pf_mpi_setup(tcomm, pf)
  call pf_pfasst_setup(pf)

  !!!! initialize advection/diffusion
  call carray4_create(q0, pf%levels(nlevs)%shape)
  call vortex_sheets(q0)

  do l = 1, nlevs
     call pf_add_hook(pf, l, PF_POST_SWEEP, project_hook)
  end do


  !!!! run
  if (pf%rank == 0) then
     print *, 'NX:       ', nx(first:)
     print *, 'NLEVS:    ', nlevs
     print *, 'NNODES:   ', nnodes(first:)
     print *, 'NTHREADS: ', nthreads
     print *, 'NSTEPS:   ', nsteps
     print *, 'NPROCS:   ', pf%comm%nproc
  end if

  call pf_pfasst_run(pf, c_loc(q0), dt, 0.0d0, nsteps=nsteps)


  !!!! done
  call mpi_barrier(pf%comm%comm, ierror)

  do l = 1, nlevs
     call feval_destroy(pf%levels(l)%levelctx)
  end do

  ! call destroy(q0)
  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(tcomm)
  call feval_finalize()

  call mpi_finalize(ierror)

end program fpfasst
