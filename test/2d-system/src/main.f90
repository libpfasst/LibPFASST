!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the 2-d advection diffusion problem on a periodic domain

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
    type(ndsysarray)     :: y_0      !<  the initial condition
    type(ndsysarray)     :: y_end    !<  the solution at the final time
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
       allocate(ndsysarray_factory::pf%levels(l)%ulevel%factory)


       !>  Allocate the arr_shape array for level (here just two dimension)
       allocate(pf%levels(l)%shape(3))
       pf%levels(l)%shape(1) = nx(l)
       pf%levels(l)%shape(2) = ny(l)
       pf%levels(l)%shape(3) = 2

       !>  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)
       print *,'calling setup',pf%levels(l)%shape
       call sweeper_setup(pf%levels(l)%ulevel%sweeper, pf%levels(l)%shape)

       !>  Set the size of the send/receive buffer
       pf%levels(l)%mpibuflen  = product(pf%levels(l)%shape)
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    !> add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

    !>  output the run options 
    call pf_print_options(pf,un_opt=6)

    !>  output the run options 
    call ad_print_options(pf,un_opt=6)

    !> allocate initial and final solutions
    call ndsysarray_build(y_0,  pf%levels(pf%nlevels)%shape)
    call ndsysarray_build(y_end, pf%levels(pf%nlevels)%shape)

    !> compute initial condition
    call initial(y_0)

    !> do the PFASST stepping
    call pf_pfasst_run(pf, y_0, dt, 0.d0, nsteps,y_end)

    !>  wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)
    
    !>  deallocate initial condition and final solution
    call ndsysarray_destroy(y_0)
    call ndsysarray_destroy(y_end)

    !> deallocate sweeper stuff
    do l = 1, pf%nlevels
       call sweeper_destroy(pf%levels(l)%ulevel%sweeper)
    end do

    !>  deallocate pfasst structure
    call pf_pfasst_destroy(pf)
  end subroutine run_pfasst

end program
