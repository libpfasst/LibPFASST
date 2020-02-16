!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the PDEs in periodic domains using a
!!  spectral representation of the solution
!! 
!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
program main
  use pfasst         !<  This module has include statements for the main pfasst routines
  use pf_mod_mpi

  integer ::  ierror

  !> initialize MPI
  call mpi_init(ierror)
  if (ierror /= 0) &
       call pf_stop(__FILE__,__LINE__,'Can not initialize mpi ',ierror)


  !> call the routine to do PFASST
  call run_pfasst()

  !> close mpi
  call mpi_finalize(ierror)

contains
  !>  This subroutine implements pfasst to PDEs in spectral space
  subroutine run_pfasst()  
    use pf_my_sweeper  !<  Local module for defining sweeper and function evaluations
    use pf_my_level    !< Local module defining the levels and restriction/interpolation
    use pf_mod_zndarray  !< Libpfasst encapsulation for complex N-dim arrays
    use hooks   !<  Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters
    use pf_mod_zutils  !< module with a few helpful subroutines
    implicit none

    !>  Local variables
    type(pf_pfasst_t)   :: pf       !<  the main pfasst structure
    type(pf_comm_t)     :: comm     !<  the communicator (here it is mpi)
    type(pf_zndarray_t) :: y_0      !<  the initial condition
    type(pf_zndarray_t) :: y_end    !<  the solution at the final time
    character(256)      :: pf_fname   !<  file name for input parameters
    class(my_sweeper_t), pointer :: sweeper

    integer           ::  l   !  loop variable over levels

    !> read problem parameters
    call probin_init(pf_fname)

    !>  set up communicator
    call pf_mpi_create(comm, MPI_COMM_WORLD)

    !>  create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)

    !> loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Assign the factory for making the solution encapsulation
       allocate(pf_zndarray_factory_t::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !>  Allocate the shape array for level 
       grid_size=nx(l)
       call pf_level_set_size(pf,l,grid_size,2*PRODUCT(grid_size)) ! 2 because it is complex
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    !> add some hooks for output
!    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
!    call pf_add_hook(pf, -1, PF_POST_ITERATION, set_error)
!    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
!    call pf_add_hook(pf, -1, PF_POST_SWEEP, set_error)
!    call pf_add_hook(pf, -1, PF_POST_ALL, echo_error)    
    call pf_add_hook(pf, -1, PF_POST_ALL, set_error)    
    call pf_add_hook(pf, -1, PF_POST_ALL, echo_error)    

    !>  output the run options 
    call pf_print_options(pf,un_opt=6)

    !>  output local parameters
    call print_loc_options(pf,un_opt=6)
    
    !> allocate initial and final solutions
    grid_size=nx(pf%nlevels)
    call zndarray_build(y_0, grid_size)
    call zndarray_build(y_end, grid_size)

    !> compute initial condition on finest level (owned by sweeper)
    sweeper => as_my_sweeper(pf%levels(pf%nlevels)%ulevel%sweeper)        
    call set_ic(sweeper,y_0)

    !> do the PFASST stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)
    
    !>  wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)

    !>  deallocate initial condition and final solution
    call zndarray_destroy(y_0)
    call zndarray_destroy(y_end)

    !>  deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst

  !> Routine to set initial condition.
  subroutine set_ic(this,y_0)
    use pf_my_sweeper  !<  Local module for defining sweeper and function evaluations
    use pf_mod_zutils    

    class(my_sweeper_t), intent(inout) :: this
    type(pf_zndarray_t), intent(inout) :: y_0

    call exact(this%fft_tool,0.0_pfdp, y_0)    
    
  end subroutine set_ic
end program
