!
! This file is part of LIBPFASST.
!
!> Simple example of using LIBPFASST.
!!
!!  This program solves the 1-d advection diffusion problem on a periodic domain

!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
program main
  use pf_mod_mpi

  integer ::  ierror

  !> initialize MPI
  call mpi_init(ierror)
  if (ierror /= 0) &
       stop "ERROR: Can't initialize MPI."

  !> call the advection-diffusion solver 
  call ad()

  !> close mpi
  call mpi_finalize(ierror)

contains
  !>  This subroutine implements pfasst to solve the advection diffusion equation
  subroutine ad()  
    use pfasst  !<  This module has include statements for the main pfasst routines
    use feval   !<  Local module for function evaluations
    use hooks   !<  Local module for diagnostics and i/o
    use probin      !< Local module reading/parsing problem parameters

    implicit none

    !>  Local variables
    type(pf_pfasst_t)              :: pf       !<  the main pfasst structure
    type(pf_comm_t)                :: comm     !<  the communicator (here is is mpi)
    type(ndarray)                  :: q0       !<  the initial condition
    type(ndarray)                  :: qend     !<  the solution at the final time
    character(256)                 :: probin_fname       !<  file name for input parameters

    integer                        ::  l   !  loop variable over levels


    !> read problem parameters
    probin_fname = "probin.nml"
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
    call probin_init(probin_fname)

    !>  set up communicator
    call pf_mpi_create(comm, MPI_COMM_WORLD)

    !>  create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=probin_fname)


    !> loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       pf%levels(l)%nsweeps = nsweeps(l)
       pf%levels(l)%nsweeps_pred = nsweeps_pred(l)

       pf%levels(l)%nnodes = nnodes(l)

       !  Allocate the user specific level object
       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)

       !  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)
       call sweeper_setup(pf%levels(l)%ulevel%sweeper, nx(l))


       !  Allocate the shape array for level (here just one dimension)
       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = nx(l)
       !  Set the size of the send/receive buffer
       pf%levels(l)%mpibuflen  = nx(l)

    end do

    !>  Set up some parameters
    call pf_pfasst_setup(pf)

    call pf_print_options(pf,6,.false.)

    !> allocate starting and end solutions
    call ndarray_build(q0, [ nx(pf%nlevels) ])
    call ndarray_build(qend, [ nx(pf%nlevels) ])    

    !> compute initial condition
    call initial(q0)

    !> add some hooks
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_residual)
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

    !> do the PFASST stepping
    call pf_pfasst_run(pf, q0, dt, 0.d0, nsteps,qend)

    !>  deallocate initial condition
    print *,'cleaning up'
    call ndarray_destroy(q0)
    call ndarray_destroy(qend)
    
    !>  deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine ad

end program
