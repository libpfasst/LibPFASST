!
! Simple example of using LIBPFASST.
!

program main

  call ad

contains

  subroutine ad()
  use pfasst
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use feval
  use hooks
  use transfer

  implicit none

  integer, parameter :: maxlevs = 3

  type(pf_pfasst_t)             :: pf
  type(pf_comm_t)               :: comm
  type(ndarray_factory), target :: factory
  type(ndarray), allocatable    :: q0
  type(ad_sweeper_t), target    :: sweepers(maxlevs)
  integer                       :: err, nvars(maxlevs), nnodes(maxlevs), l
  double precision              :: dt


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
  nnodes = [ 3, 5, 9 ]       ! number of sdc nodes on time/space levels
  dt     = 0.005_pfdp

  ! call ndarray_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, comm, maxlevs)

  pf%qtype  = SDC_GAUSS_LOBATTO
  pf%niters = 4

  if (pf%nlevels > 1) then
     pf%levels(1)%nsweeps = 3
  end if

  do l = 1, pf%nlevels
     pf%levels(l)%nsweeps = 1

     pf%levels(l)%nvars  = nvars(maxlevs-pf%nlevels+l)
     pf%levels(l)%nnodes = nnodes(maxlevs-pf%nlevels+l)

     pf%levels(l)%factory     => factory
     pf%levels(l)%sweeper     => sweepers(l)
     call sweepers(l)%setup(pf%levels(l)%nvars)

     allocate(pf%levels(l)%shape(1))
     pf%levels(l)%shape(1)    = pf%levels(l)%nvars
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  !
  ! compute initial condition, add hooks, run
  !

  allocate(q0)
  call ndarray_build(q0, [ pf%levels(pf%nlevels)%nvars ])
  call initial(q0)

  call pf_add_hook(pf, pf%nlevels, PF_POST_ITERATION, echo_error)
  call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
  call pf_pfasst_run(pf, q0, dt, tend=0.d0, nsteps=1*comm%nproc)

  !
  ! cleanup
  !
  call pf_pfasst_destroy(pf)    ! XXX
  call pf_mpi_destroy(comm)     ! XXX
  call mpi_finalize(err)
  call fftw_cleanup()

end subroutine ad

end program main
