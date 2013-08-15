!
! Simple example of using LIBPFASST.
!

program main
  use pfasst
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use feval
  use hooks
  use transfer

  implicit none

  integer, parameter :: nlevs = 3

  type(pf_pfasst_t)           :: pf
  type(pf_comm_t)             :: comm
  type(pf_sweeper_t), target  :: sweeper
  type(pf_encap_t),   target  :: encap
  type(ndarray),      pointer :: q0
  integer                     :: argc, err, nvars(nlevs), nnodes(nlevs), l
  double precision            :: dt
  character(len=128)          :: arg


  !
  ! initialize mpi
  !

  call mpi_init(err)
  if (err .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst using three levels
  !

  nvars  = [ 32, 64, 128 ]   ! number of dofs on the time/space levels
  nnodes = [ 2, 3, 5 ]       ! number of sdc nodes on time/space levels
  dt     = 0.1_pfdp

  call ndarray_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)
  call pf_pfasst_create(pf, comm, nlevs)

  pf%qtype  = SDC_GAUSS_LOBATTO
  pf%niters = 12

  do argc = 1, command_argument_count()
     call get_command_argument(argc, arg)
     select case(arg)
     case ("--ring")
        pf%window      = PF_WINDOW_RING
        pf%abs_res_tol = 1.d-10
     case default
        stop "Usage: main.exe [--ring]"
     end select
  end do
  
  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     pf%levels(l)%nvars  = nvars(l)
     pf%levels(l)%nnodes = nnodes(l)

     call feval_create_workspace(pf%levels(l)%ctx, pf%levels(l)%nvars)

     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%encap       => encap
     pf%levels(l)%sweeper     => sweeper

     allocate(pf%levels(l)%shape(1))
     pf%levels(l)%shape(1)    = nvars(l)
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  if (pf%rank == 0) then
     print *, 'nvars: ', pf%levels(:)%nvars
     print *, 'nnodes:', pf%levels(:)%nnodes
  end if


  !
  ! compute initial condition, add hooks, run
  !

  allocate(q0)
  call ndarray_create_simple(q0, [ nvars(3) ])
  call initial(q0)

  pf%abs_res_tol = 1.d-10

  call pf_add_hook(pf, nlevs, PF_POST_ITERATION, echo_error)
  call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
  call pf_pfasst_run(pf, c_loc(q0), dt, tend=0.d0, nsteps=4*comm%nproc)


  !
  ! cleanup
  !

  call ndarray_destroy(c_loc(q0))

  do l = 1, nlevs
     call feval_destroy_workspace(pf%levels(l)%ctx)
  end do

  call pf_imex_destroy(sweeper)
  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(err)
  call fftw_cleanup()

end program main
