!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the linear scalar equation
!!
!!    y'=lambda*y
!!
!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
program main
  use pf_mod_mpi

  integer ::  ierror

  !> Initialize MPI
  call mpi_init(ierror)
  if (ierror /= 0) &
       stop "ERROR: Can't initialize MPI."

  !> Call the  solver 
  call run_pfasst()

  !> Close mpi
  call mpi_finalize(ierror)

contains
  !>  This subroutine set ups and calls libpfasst 
  subroutine run_pfasst()  
    use pfasst            !< This module has include statements for the main pfasst routines
    use pf_my_sweeper     !< Local module for sweeper
    use pf_my_level       !< Local module for level
    use hooks             !< Local module for diagnostics and i/o
    use probin            !< Local module reading/parsing problem parameters
    use encap             !< Local module defining the encapsulation

    implicit none

    !>  Local variables
    type(pf_pfasst_t) :: pf        !<  the main pfasst structure
    type(pf_comm_t)   :: comm      !<  the communicator (here it is mpi)
    type(scalar_encap) :: y_0      !<  the initial condition
    type(scalar_encap) :: y_end    !<  the solution at the final time
    character(256)    :: pf_fname   !<  file name for input of PFASST parameters

    integer           ::  l   !  loop variable over levels

    print *,'2'
    !> Read problem parameters
    call probin_init(pf_fname)

    print *,'3'
    !>  Set up communicator
    call pf_mpi_create(comm, MPI_COMM_WORLD)

    print *,'4'    
    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)

    print *,'5'    
    !> Loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       print *,'6'       
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       print *,'7'
       !>  Allocate the user specific data constructor
       allocate(scalar_factory::pf%levels(l)%ulevel%factory)

       print *,'8'
       !>  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       print *,'9'
       !>  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,[1])
    end do

       print *,'10'    
    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    !>  Add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

    !>  Output run parameters
    call print_loc_options(pf,un_opt=6)

    print *,'11'    
    !> Set the initial condition
    call y_0%setval(1.0_pfdp)

    print *,'12'    
    !> Do the PFASST time stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)

    !>  Wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)
    
    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst

end program
