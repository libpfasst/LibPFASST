!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the 1-d advection diffusion problem on a periodic domain

!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
program main
  use pf_mod_mpi

  integer ::  ierror

  !> Initialize MPI
  call mpi_init(ierror)
  if (ierror /= 0) &
       stop "ERROR: Can't initialize MPI."

  !> Call the advection-diffusion solver 
  call run_pfasst()

  !> Close mpi
  call mpi_finalize(ierror)

contains
  !>  This subroutine implements pfasst to solve the advection diffusion equation
  subroutine run_pfasst()  
    use pfasst  !< This module has include statements for the main pfasst routines
    use pf_my_sweeper  !< Local module for sweeper and function evaluations
    use pf_my_level    !< Local module for the levels
    use hooks   !< Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters

    implicit none

    !>  Local variables
    type(pf_pfasst_t) :: pf       !<  the main pfasst structure
    type(pf_comm_t)   :: comm     !<  the communicator (here it is mpi)
    type(ndarray)     :: y_0      !<  the initial condition
    type(ndarray)     :: y_end    !<  the solution at the final time
    character(256)    :: probin_fname   !<  file name for input parameters

    integer           ::  l   !  loop variable over levels

    !> Set the name of the input file
    probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
    
    !> Read problem parameters
    call probin_init(probin_fname)

    !>  Set up communicator
    call pf_mpi_create(comm, MPI_COMM_WORLD)

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=probin_fname)

    !> Loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(ad_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data factory
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !>  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,[nx(l)])
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)


    !> Add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
    call pf_add_hook(pf, -1, PF_POST_SWEEP, pf_echo_residual)

    !>  Output the run options 
    call pf_print_options(pf,un_opt=6)

    !>  Output local parameters
    call print_loc_options(pf,un_opt=6)
    
    !> Allocate initial and final solutions
    call ndarray_build(y_0, [ nx(pf%nlevels) ])
    call ndarray_build(y_end, [ nx(pf%nlevels) ])    

    !> compute initial condition
    call initial_sol(y_0)

    !> Do the PFASST stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)
    
    !>  Wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)

    !>  Deallocate initial condition and final solution
    call ndarray_destroy(y_0)
    call ndarray_destroy(y_end)

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst

end program
