!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the 1-d advection diffusion problem on a periodic domain

!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
program main
  use pfasst  !<  This module has include statements for the main pfasst routines
  use pf_mod_mpi

  integer ::  ierror

  !> initialize MPI
  call mpi_init(ierror)
  if (ierror /= 0) &
       call pf_stop(__FILE__,__LINE__,'Can not initialize MPI, ierrer',ierror)

  !> call the advection-diffusion solver 
  call run_pfasst()

  !> close mpi
  call mpi_finalize(ierror)

contains
  !>  This subroutine implements pfasst to solve the advection diffusion equation
  subroutine run_pfasst()  
    use, intrinsic :: iso_c_binding
    use pf_mod_parareal  !<  Parareal routines
    use pf_my_stepper   !<  Local module for function evaluations
    use pf_my_level    !< Local module for the levels
    use hooks   !<  Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters

    implicit none

    !>  Local variables
    type(pf_pfasst_t)  :: pf       !<  the main pfasst structure
    type(pf_comm_t)    :: comm     !<  the communicator (here it is mpi)
    type(pf_ndarray_t):: y_0      !<  the initial condition
    type(pf_ndarray_t):: y_end    !<  the solution at the final time
    character(256)     :: pf_fname   !<  file name for input parameters
    class(my_stepper_t), pointer :: stepper
    real(c_double),        pointer :: yvec(:)     ! underlying vector
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
       allocate(pf_ndarray_factory_t::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(my_stepper_t::pf%levels(l)%ulevel%stepper)

       !>  Allocate the shape array for level (here just one dimension)
       grid_size=nx(l)
       call pf_level_set_size(pf,l,grid_size,PRODUCT(grid_size))
    end do

    !>  Set up some pfasst stuff
    pf%use_sdc_sweeper=.FALSE.
    call pf_pfasst_setup(pf)

    !> add some hooks for output
        call pf_add_hook(pf, -1, PF_POST_PREDICTOR, echo_error)
        !        call pf_add_hook(pf, 2, PF_POST_ITERATION, echo_error)
        call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)        
!    call pf_add_hook(pf, -1, PF_POST_CONVERGENCE, echo_error)    

    !>  output the run options 
    call pf_print_options(pf,un_opt=6)

    !>  Output local parameters
    call print_loc_options(pf,un_opt=6)
    
    !> allocate initial and final solutions
    call ndarray_build(y_0, grid_size)
    call ndarray_build(y_end, grid_size)    

    !> compute initial condition
    stepper => as_my_stepper(pf%levels(pf%nlevels)%ulevel%stepper)            
    call set_ic(stepper,y_0)

    !> do the PFASST stepping
    call pf_parareal_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)    
    
    !>  wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)

    !>  deallocate initial condition and final solution
    call ndarray_destroy(y_0)
    call ndarray_destroy(y_end)

    !>  deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst

  !> Routine to set initial condition.
  subroutine set_ic(this,y_0)
    use pf_my_stepper   !<  Local module for function evaluations
    use pf_mod_rutils    

    class(my_stepper_t), intent(inout) :: this
    type(pf_ndarray_t), intent(inout) :: y_0

    call exact(this%fft_tool,0.0_pfdp, y_0)    
    
  end subroutine set_ic

  
end program
