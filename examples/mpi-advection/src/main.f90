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

  integer, parameter :: maxlevs = 3

  type(pf_pfasst_t)        :: pf
  type(pf_comm_t)          :: comm
  type(pf_encap_t), target :: encap
  type(ndarray), pointer   :: q0
  integer                  :: err, nvars(maxlevs), nnodes(maxlevs), l
  double precision         :: dt


  !
  ! initialize mpi
  !

  call mpi_init(err)
  if (err .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst
  !

  nvars  = [ 32, 64, 128 ]   ! number of dofs on the time/space levels
!  nvars  = [ 128, 128, 128 ]   ! number of dofs on the time/space levels
  nnodes = [ 3, 5, 9 ]       ! number of sdc nodes on time/space levels
!  nvars  = [  16 ]   ! number of dofs on the time/space levels
!  nnodes = [ 2 ]       ! number of sdc nodes on time/space levels
  dt     = 0.005_pfdp

  call ndarray_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, comm, maxlevs)

  pf%qtype  = SDC_GAUSS_LOBATTO
  pf%niters = 32

  if (pf%nlevels > 1) then
     pf%levels(1)%nsweeps = 3
  end if

  do l = 1, pf%nlevels
     pf%levels(l)%nsweeps = 1

     pf%levels(l)%nvars  = nvars(maxlevs-pf%nlevels+l)
     pf%levels(l)%nnodes = nnodes(maxlevs-pf%nlevels+l)

     call feval_create_workspace(pf%levels(l)%levelctx, pf%levels(l)%nvars)

     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%encap       => encap
    call pf_imexQ_create(pf%levels(l)%sweeper, eval_f1, eval_f2, comp_f2)
!     call pf_implicitQ_create(pf%levels(l)%sweeper,  eval_f2, comp_f2)

     allocate(pf%levels(l)%shape(1))
     pf%levels(l)%shape(1)    = pf%levels(l)%nvars
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  if (pf%rank == 0) then
     ! print *, 'nvars: ', pf%levels(:)%nvars
     ! print *, 'nnodes:', pf%levels(:)%nnodes
     call pf_print_options(pf)
  end if


  !
  ! compute initial condition, add hooks, run
  !

  allocate(q0)
  call ndarray_create_simple(q0, [ pf%levels(pf%nlevels)%nvars ])
  call initial(q0)

  if (pf%window == PF_WINDOW_RING) pf%abs_res_tol = 1.d-9

  ! call pf_cycle_print(pf)

  call pf_add_hook(pf, pf%nlevels, PF_POST_ITERATION, echo_error)
  ! call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
  call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
  call pf_pfasst_run(pf, c_loc(q0), dt, tend=0.d0, nsteps=1*comm%nproc)

  !
  ! cleanup
  !

  call ndarray_destroy(c_loc(q0))

  do l = 1, pf%nlevels
     call feval_destroy_workspace(pf%levels(l)%levelctx)
  end do


  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(err)
  call fftw_cleanup()

end program main
