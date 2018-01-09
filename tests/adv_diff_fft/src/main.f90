!
! Simple example of using LIBPFASST.
!

program main
  use pf_mod_mpi!, only: mpi_init, mpi_finalize
!  use feval, only: fftw_cleanup
  integer ::  ierror

  call mpi_init(ierror)
  if (ierror /= 0) &
       stop "ERROR: Can't initialize MPI."
  call ad()
!  call fftw_cleanup()
  call mpi_finalize(ierror)

contains

  subroutine ad()
    use pfasst
    use feval
    use hooks
    use pf_mod_mpi
    use probin      !< should be library module for reading/parsing problem parameters
    implicit none


    type(pf_pfasst_t)              :: pf
    type(pf_comm_t)                :: comm
    type(ndarray), allocatable     :: q0
    integer                        ::  l   !  Loop variable
    character(256) :: probin_fname       !<  file name for input
    class(ad_sweeper_t), pointer   :: ad_sweeper_ptr

    !
    ! initialize pfasst
    !

!    nvars  = [ 128 ]   ! number of dofs on the time/space levels
!    nnodes = [ 5 ]       ! number of sdc nodes on time/space levels

!    nvars  = [ 32, 64,  128]   ! number of dofs on the time/space levels
!    nnodes = [ 3,5,9 ]       ! number of sdc nodes on time/space levels
!    dt     = 0.05_pfdp

      probin_fname = "probin.nml"
      if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
      call probin_init(probin_fname)
      print *, nvars
      print *, nnodes
      print *, nsteps,Tfin,dt

    call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_pfasst_create(pf, comm, fname=probin_fname)

!    pf%qtype  = SDC_GAUSS_LOBATTO
!    pf%niters = 20
!    pf%abs_res_tol=1.0D-15    
!    pf%rel_res_tol=1.0D-12

    do l = 1, pf%nlevels
       pf%levels(l)%nsweeps = nsweeps(l)
       pf%levels(l)%nsweeps_pred = nsweeps_pred(l)

       pf%levels(l)%nvars  = nvars(l)
       pf%levels(l)%nnodes = nnodes(l)

       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)

       call setup(pf%levels(l)%ulevel%sweeper, pf%levels(l)%nvars)

       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = pf%levels(l)%nvars
    end do

    if (pf%nlevels > 1) then
       pf%levels(1)%nsweeps = 1 
       pf%levels(1)%nsweeps_pred = 1
    end if

    call pf_mpi_setup(comm, pf,ierror) ! XXX: move this into pf_pfasst_setup
    call pf_pfasst_setup(pf)

    !
    ! compute initial condition, add hooks, run
    !
    allocate(q0)
    call ndarray_build(q0, [ pf%levels(pf%nlevels)%nvars ])
    call initial(q0)

    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_residual)
    call pf_pfasst_run(pf, q0, dt, 0.d0, nsteps)

    deallocate(q0%flatarray)
    deallocate(q0%shape)
    deallocate(q0)

    !
    ! cleanup
    !
    call pf_mpi_destroy(comm)     ! XXX
    call pf_pfasst_destroy(pf)

  end subroutine ad

end program
