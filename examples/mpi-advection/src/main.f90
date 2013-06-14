!
! Simple example of using LIBPFASST.
!

program main
  use pf_mod_dtype
  use pf_mod_parallel
  use pf_mod_pfasst
  use pf_mod_mpi
  use pf_mod_comm_mpi

  use feval
  use hooks
  use transfer
  use encap_array1d

  implicit none

  type(pf_pfasst_t)  :: pf
  type(pf_comm_t)    :: comm
  type(pf_sweeper_t), target :: sweeper
  type(pf_encap_t),   target :: encap
  integer            :: ierror, nlevs, nvars(3), nnodes(3), l
  double precision   :: dt

  type(array1d), target :: q0


  !
  ! initialize mpi
  !
  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst using three levels
  !

  nvars  = [ 32, 64, 128 ]      ! number of dofs on the time/space levels
  nnodes = [ 2, 3, 5 ]          ! number of sdc nodes on time/space levels
  dt     = 0.1_pfdp
  nlevs  = 3

  call array1d_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)
  call pf_pfasst_create(pf, comm, nlevs)

  pf%niters = 12
  pf%qtype  = SDC_GAUSS_LOBATTO + SDC_PROPER_NODES

  pf%echo_timings = .false.

  pf%window      = PF_WINDOW_RING
  pf%abs_res_tol = 1.d-8

  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     pf%levels(l)%nvars  = nvars(3-nlevs+l)
     pf%levels(l)%nnodes = nnodes(3-nlevs+l)

     call feval_create_workspace(pf%levels(l)%ctx, pf%levels(l)%nvars)

     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%encap       => encap
     pf%levels(l)%sweeper     => sweeper
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  if (pf%rank == 0) then
     print *, 'nvars: ', pf%levels(:)%nvars
     print *, 'nnodes:', pf%levels(:)%nnodes
  end if


  !
  ! run
  !
  allocate(q0%array(pf%levels(nlevs)%nvars))
  call initial(q0)

  ! call pf_logger_attach(pf)
  call pf_add_hook(pf, nlevs, PF_POST_ITERATION, echo_error)
  ! call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
  call pf_pfasst_run(pf, c_loc(q0), dt, 0.0_pfdp, nsteps=2*comm%nproc)


  !
  ! cleanup
  !
  deallocate(q0%array)

  do l = 1, nlevs
     call feval_destroy_workspace(pf%levels(l)%ctx)
  end do

  call pf_imex_destroy(sweeper)
  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(ierror)
  call fftw_cleanup()

end program main
