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

  !> initialize MPI
  call mpi_init(ierror)
  if (ierror /= 0) &
       stop "ERROR: Can't initialize MPI."

  !> call the advection-diffusion solver 
  call run_pfasst()

  !> close mpi
  call mpi_finalize(ierror)

contains
  !>  This subroutine implements pfasst to solve the advection diffusion equation
  subroutine run_pfasst()  
    use pfasst  !<  This module has include statements for the main pfasst routines
    use feval   !<  Local module for function evaluations
    use hooks   !<  Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters

    implicit none

    !>  Local variables
    type(pf_pfasst_t) :: pf       !<  the main pfasst structure
    type(pf_comm_t)   :: comm     !<  the communicator (here it is mpi)
    type(ndarray)     :: y_0      !<  the initial condition
    type(ndarray)     :: y_end    !<  the solution at the final time
    character(256)    :: probin_fname   !<  file name for input parameters

    integer           ::  l   !  loop variable over levels

    !> set the name of the input file
    probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
    
    !> read problem parameters
    call probin_init(probin_fname)

    !>  set up communicator
    call pf_mpi_create(comm, MPI_COMM_WORLD)

    !>  create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=probin_fname)

    !> loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)

       !>  Allocate the shape array for level (here just one dimension)
       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = nx(l)

       !>  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)
       call sweeper_setup(pf%levels(l)%ulevel%sweeper, pf%levels(l)%shape)


       !>  Set the size of the send/receive buffer
       pf%levels(l)%mpibuflen  = nx(l)
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)


    !> add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)

    !>  output the run options 
    call pf_print_options(pf,un_opt=6)

    !>  Output local parameters
    call ad_print_options(pf,un_opt=6)
    
    !> allocate initial and final solutions
    call ndarray_build(y_0, [ nx(pf%nlevels) ])
    call ndarray_build(y_end, [ nx(pf%nlevels) ])    

    !> compute initial condition
    call initial(y_0)

    !> do the PFASST stepping
    call pf_pfasst_run(pf, y_0, dt, 0.d0, nsteps,y_end)
    
    !>  deallocate initial condition and final solution
    call ndarray_destroy(y_0)
    call ndarray_destroy(y_end)
    
    !>  deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst

end program
