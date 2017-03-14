!
! Simple example of using LIBPFASST.
!

program main
  use pf_mod_mpi!, only: mpi_init, mpi_finalize
  integer :: err

  call mpi_init(err)
  if (err /= 0) &
       stop "ERROR: Can't initialize MPI."
  call ad()

  call mpi_finalize(err)

contains

  subroutine ad()
    use pfasst
    use feval
    use hooks
    use pf_mod_mpi

    implicit none

    integer, parameter :: maxlevs = 1

    type(pf_pfasst_t)              :: pf
    type(pf_comm_t)                :: comm
    type(ndarray), allocatable     :: q0
    integer                        :: nvars(maxlevs), nnodes(maxlevs), l
    double precision               :: dt

    real :: val

    !
    ! initialize pfasst
    !

    nvars  = [ 1 ]   ! number of dofs on the time/space levels
    nnodes = [ 9 ]       ! number of sdc nodes on time/space levels
    dt     = 0.1_pfdp

    call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_pfasst_create(pf, comm, maxlevs)

    pf%qtype  = SDC_GAUSS_LOBATTO
    pf%niters = 30

    do l = 1, pf%nlevels
       pf%levels(l)%nsweeps = 1

       pf%levels(l)%nvars  = nvars(maxlevs-pf%nlevels+l)
       pf%levels(l)%nnodes = nnodes(maxlevs-pf%nlevels+l)

       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)

       call setup(pf%levels(l)%ulevel%sweeper, pf%levels(l)%nvars)

       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = pf%levels(l)%nvars
    end do

    if (pf%nlevels > 1) then
       pf%levels(1)%nsweeps = 2 
       pf%levels(1)%nsweeps_pred = 2
    end if

    call pf_mpi_setup(comm, pf) ! XXX: move this into pf_pfasst_setup
    call pf_pfasst_setup(pf)

    !
    ! compute initial condition, add hooks, run
    !
    allocate(q0)
    call ndarray_build(q0, [ pf%levels(pf%nlevels)%nvars ])
    call initial(q0)

    val = q0%flatarray(1)
    print *, val

    call pf_add_hook(pf, pf%nlevels, PF_POST_ITERATION, echo_error)
    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
    call pf_pfasst_run(pf, q0, dt, tend=0.d0, nsteps=1*comm%nproc)

    val = q0%flatarray(1)
    print *, val

    !
    ! cleanup
    !
    call pf_mpi_destroy(comm)     ! XXX

  end subroutine ad

end program
