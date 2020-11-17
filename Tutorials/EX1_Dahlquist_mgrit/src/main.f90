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
  !>  This subroutine setups and calls libpfasst 
  subroutine run_pfasst()  
    use pfasst        !<  This module has include statements for the main pfasst routines
    use pf_my_stepper !<  Local module stepper
    use pf_my_level   !<  Local module for level
    use probin       !<  Local module reading/parsing problem parameters
    use pf_mod_MGRIT
    use pf_mod_parareal
    use hooks

    implicit none

    !>  Local variables
    type(pf_pfasst_t) :: pf       !<  the main pfasst structure
    type(pf_comm_t)   :: comm     !<  the communicator (here it is mpi)
    type(pf_ndarray_t):: y_0 , y_end   !<  the initial condition and solution
    character(256)    :: pf_fname   !<  file name for input of PFASST parameters

    integer           ::  l   !  loop variable over levels
    real(pfdp) :: T0
    type(mgrit_level_data), allocatable :: mg_ld(:)
    integer :: n_coarse, refine_factor

    !> Read problem parameters
    call probin_init(pf_fname)

    !>  Set up communicator
    call pf_mpi_create(comm, MPI_COMM_WORLD)

    pf%use_rk_stepper = .true.
    pf%use_sdc_sweeper = .false.

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)

    !> Loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data constructor
       allocate(pf_ndarray_factory_t::pf%levels(l)%ulevel%factory)

       !>  Allocate the stepper at this level
       allocate(my_stepper_t::pf%levels(l)%ulevel%stepper)

       !>  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,[1])

       pf%levels(l)%ulevel%stepper%order = 1
       pf%levels(l)%ulevel%stepper%nsteps = nsteps_rk(l)
       !pf%rk_order(l) = 1
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    T0 = 0.0_pfdp
    n_coarse = max(1, mgrit_n_coarse/pf%comm%nproc)
    refine_factor = mgrit_refine_factor
    call mgrit_initialize(pf, mg_ld, T0, Tfin, n_coarse, refine_factor)

    do l = 1, pf%nlevels
       !mg_ld(l)%FCF_flag = .false.
       !pf%levels(l)%ulevel%stepper%nsteps = mg_ld(l)%Nt
    end do

    !> add some hooks for output  (using a LibPFASST hook here)
    !call pf_add_hook(pf, -1, PF_POST_ITERATION, pf_echo_residual)
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

    !>  Output run parameters to screen
    call print_loc_options(pf,un_opt=6)
    
    !>  Allocate initial consdition
    call ndarray_build(y_0, [ 1 ])
    call ndarray_build(y_end, [ 1 ])

    !> Set the initial condition 
    call y_0%setval(1.0_pfdp)
  
    if (mgrit_flag .eqv. .true.) then 
       call pf_MGRIT_run(pf, mg_ld, y_0, y_end)
       if (pf%rank .eq. pf%comm%nproc-1) call y_end%eprint()
    else
       call pf_parareal_run(pf, y_0, dt, Tfin, nsteps, y_end)
       if (pf%rank .eq. pf%comm%nproc-1) call y_end%eprint()
    end if
    
    !>  Wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)

    !>  Deallocate initial condition and final solution
    call ndarray_destroy(y_0)
    call ndarray_destroy(y_end)
    
    !>  Deallocate pfasst structure
    !call pf_pfasst_destroy(pf)

  end subroutine run_pfasst

end program
