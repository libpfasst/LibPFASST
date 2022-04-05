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
    use my_sweeper  !< Local module for sweeper and function evaluations
    use my_stepper
    use my_level    !< Local module for the levels
    use hooks   !< Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters
    use pf_mod_MGRIT
    use pf_mod_parareal

    implicit none

    !>  Local variables
    type(pf_pfasst_t) :: pf       !<  the main pfasst structure
    type(pf_comm_t)   :: comm     !<  the communicator (here it is mpi)
    type(pf_ndarray_t):: y_0      !<  the initial condition
    type(pf_ndarray_t):: y_end    !<  the solution at the final time
    character(256)    :: pf_fname   !<  file name for input of PFASST parameters

    integer           ::  l   !  loop variable over levels
    type(mgrit_level_data), allocatable :: mg_ld(:)
    integer :: mgrit_nsteps, coarsen_factor
    real(pfdp) :: T0
    logical :: FAS_flag, FCF_flag, setup_start_coarse_flag

    real(pfdp), allocatable :: yexact_1d(:), yexact_2d(:,:), yexact_3d(:,:,:)
    real(pfdp), pointer :: y_0_1d(:), y_0_2d(:,:), y_0_3d(:,:,:)
    real(pfdp) :: error
    integer :: temp_nsweeps, temp_nsweeps_pred, l_nx_not_zero

    integer :: nx_old(PF_MAXLEVS)

    integer, save :: nx_coarse

    !> Read problem parameters
    call probin_init(pf_fname)

    v(1) = vx
    v(2) = vx
    v(3) = vx
    kfreq(1) = kfreqx
    kfreq(2) = kfreqx
    kfreq(3) = kfreqx
    space_len(1) = Lx
    space_len(2) = Lx
    space_len(3) = Lx

    !>  Set up communicator
    call pf_mpi_create(comm, MPI_COMM_WORLD)

    if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3) .or. (solver_type .eq. 4)) then
       pf%use_rk_stepper = .true.
       pf%use_sdc_sweeper = .false.
    else
       pf%use_rk_stepper = .false.
       pf%use_sdc_sweeper = .true.
    end if

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)
    if (((solver_type .eq. 2)  .or. (solver_type .eq. 4)) .and. (pf%nlevels .ne. 2)) then
       print *, 'ERROR: nlevels must be 2 for Parareal.'
       return
    end if

    do l = 1,pf%nlevels 
       if (nx(l) .gt. 0) then
          l_nx_not_zero = l
       end if
    end do
    nx_old = nx
    if (l_nx_not_zero .lt. pf%nlevels) then
       do l = pf%nlevels,1,-1
          if (l_nx_not_zero .gt. 0) then
             nx(l) = nx_old(l_nx_not_zero)
             l_nx_not_zero = l_nx_not_zero - 1
          else
             nx(l) = 0
          end if
       end do
    end if

    do l = pf%nlevels,1,-1
       if (nx(l) .le. 0) then
          nx(l) = nx_coarse
       else
          nx_coarse = nx(l)
       end if
       !print *,nx(l)
    end do

    !> Loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data factory
       allocate(pf_ndarray_factory_t::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3) .or. (solver_type .eq. 4)) then
          allocate(my_stepper_t::pf%levels(l)%ulevel%stepper)
       else
          allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)
       end if

       !>  Set the size of the data on this level (here just one)
       if (ndim == 1) then
          call pf_level_set_size(pf,l,[nx(l)])
       else if (ndim == 2) then
          call pf_level_set_size(pf,l,[nx(l), nx(l)])
       else if (ndim == 3) then
          call pf_level_set_size(pf,l,[nx(l), nx(l), nx(l)])
       else

       end if

       if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3) .or. (solver_type .eq. 4)) then
          pf%levels(l)%ulevel%stepper%order = rk_order
          if (solver_type .eq. 2) then
             pf%levels(l)%ulevel%stepper%nsteps = nsteps_rk(l)
          end if
       end if
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    if ((solver_type .eq. 1) .or. (solver_type .eq. 4)) then
       if (solver_type .eq. 4) then
          FCF_flag = .false.
       else
          FCF_flag = .true.
       end if

       FAS_flag = .false.
       T0 = 0.0_pfdp
       setup_start_coarse_flag = .false.
       mgrit_nsteps = max(1, nsteps/pf%comm%nproc)
       call mgrit_initialize(pf, mg_ld, T0, Tfin, mgrit_nsteps, mgrit_coarsen_factor, FAS_flag, FCF_flag, setup_start_coarse_flag)
    end if

    !> Add some hooks for output
    if (solver_type .eq. 1) then
       call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
    else
       !call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
       call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
    end if

    !>  Output the run options 
    call pf_print_options(pf,un_opt=6)

    !>  Output local parameters
    call print_loc_options(pf,un_opt=6)
    
    !> Allocate initial and final solutions

    if (ndim == 1) then
       call ndarray_build(y_0, [nx(pf%nlevels)])
       call ndarray_build(y_end, [nx(pf%nlevels)])    
    else if (ndim == 2) then
       call ndarray_build(y_0, [nx(pf%nlevels), nx(pf%nlevels)])
       call ndarray_build(y_end, [nx(pf%nlevels), nx(pf%nlevels)])
    else if (ndim == 3) then
       call ndarray_build(y_0, [nx(pf%nlevels), nx(pf%nlevels), nx(pf%nlevels)])
       call ndarray_build(y_end, [nx(pf%nlevels), nx(pf%nlevels), nx(pf%nlevels)])
    else

    end if

    !> compute initial condition
    call initial_sol(y_0)

    !if (ndim == 1) then
    !   allocate(yexact_1d(pf%levels(pf%nlevels)%lev_shape(1)))
    !   y_0_1d => get_array1d(y_0)
    !   call exact(pf%state%t0+pf%state%dt, yexact_1d)
    !   error = maxval(abs(y_0_1d-yexact_1d))
    !   deallocate(yexact_1d)
    !else if (ndim == 2) then
    !   allocate(yexact_2d(pf%levels(pf%nlevels)%lev_shape(1), &
    !                      pf%levels(pf%nlevels)%lev_shape(2)))
    !   y_0_2d => get_array2d(y_0)
    !   call exact(pf%state%t0+pf%state%dt, yexact_2d)
    !   error = maxval(abs(y_0_2d-yexact_2d))
    !   deallocate(yexact_2d)
    !else if (ndim == 3) then
    !   allocate(yexact_3d(pf%levels(pf%nlevels)%lev_shape(1), &
    !                      pf%levels(pf%nlevels)%lev_shape(2), &
    !                      pf%levels(pf%nlevels)%lev_shape(3)))
    !   y_0_3d => get_array3d(y_0)
    !   call exact(pf%state%t0+pf%state%dt, yexact_3d)
    !   error = maxval(abs(y_0_3d-yexact_3d))
    !   deallocate(yexact_3d)
    !else
    !end if

    if ((solver_type .eq. 1) .or. (solver_type .eq. 4)) then !> MGRIT and MGRIT-Parareal
       call pf_MGRIT_run(pf, mg_ld, y_0, y_end)
    else if (solver_type .eq. 2) then !> Parareal
       call pf_parareal_run(pf, y_0, dt, Tfin, nsteps, y_end)
    else if (solver_type .eq. 3) then !> Sequential solver
       pf%state%nsteps = nsteps
       pf%state%step = nsteps-1
       pf%state%iter = 1
       pf%state%sweep = 1;
       call initialize_results(pf)
       if (pf%save_timings > 0) call pf_start_timer(pf, T_TOTAL)
       call pf%levels(1)%ulevel%stepper%do_n_steps(pf, 1, T0, y_0, y_end, dt, nsteps)
       call pf%levels(pf%nlevels)%qend%copy(y_end)
       call echo_error(pf, pf%nlevels)
       if (pf%save_timings > 0) call pf_stop_timer(pf, T_TOTAL)
       call pf_dump_stats(pf)
    else
       call pf_pfasst_run(pf, y_0, dt, Tfin, nsteps, y_end)
    end if
    !if (pf%rank .eq. pf%comm%nproc-1) call y_end%eprint()

    !>  Wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)
    
    !>  Deallocate initial condition and final solution
    call ndarray_destroy(y_0)
    call ndarray_destroy(y_end)

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst

end program
