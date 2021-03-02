!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the linear hypre_vector equation
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
    use pf_space_comm
    use pfasst_hypre
    use pf_mod_parareal

    implicit none

    !>  Local variables
    type(pf_pfasst_t) :: pf        !<  the main pfasst structure
    type(pf_comm_t) :: comm      !<  the communicator (here it is mpi)
    type(hypre_vector_encap) :: y_0      !<  the initial condition
    type(hypre_vector_encap) :: y_end    !<  the solution at the final time
    type(hypre_vector_encap) :: u
    class(pf_encap_t), allocatable :: y_0_base, y_end_base, u_base
    character(256) :: pf_fname   !<  file name for input of PFASST parameters
    integer, allocatable :: lev_shape(:,:)
    type(hypre_vector_factory) :: hvf
    type(my_sweeper_t) :: s_finest, s
    type(my_level_t) :: my_lev 

    integer :: l, l_finest   !  loop variable over levels
    integer :: n, m
    integer :: space_comm, time_comm, space_color, time_color
    integer :: level_index

    integer :: nproc, rank, error
    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1
    integer :: spacial_coarsen_flag
    type(mgrit_level_data), allocatable :: mg_ld(:)

    ! check size
    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

    !> Read problem parameters
    call probin_init(pf_fname)

    n = num_grid_points * num_grid_points
    call create_simple_communicators(nspace, ntime, space_comm, time_comm, space_color, time_color, space_dim)

    !>  Set up communicator
    call pf_mpi_create(comm, time_comm)
    

    if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3)) then
       pf%use_rk_stepper = .true.
       pf%use_sdc_sweeper = .false.
    else
       pf%use_rk_stepper = .false.
       pf%use_sdc_sweeper = .true.
    end if

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)
    if ((solver_type .eq. 2) .and. (pf%nlevels .ne. 2)) then
       print *, 'ERROR: nlevels must be 2 for Parareal.'
       return
    end if

    spacial_coarsen_flag = 0
    call PfasstHypreInit(pf, mg_ld, lev_shape, space_color, time_color, spacial_coarsen_flag)
    !print *,time_color,space_color,pf%rank
    
    !>  Add some hooks for output
    if (solver_type .eq. 1) then
       call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
    else
       !call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
       call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
    end if
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

    !>  Output run parameters
    call print_loc_options(pf,un_opt=6)

    level_index = pf%nlevels
    !> Set the initial condition
    call hvf%create_single(y_0_base, level_index, lev_shape(pf%nlevels,:))
    call hvf%create_single(y_end_base, level_index, lev_shape(pf%nlevels,:))
    y_0 = cast_as_hypre_vector(y_0_base)
    y_end = cast_as_hypre_vector(y_end_base)
    ! call y_0%setval(init_cond)
    ! call y_end%setval(init_cond)
    call initial(y_0)

    !> Do the PFASST time stepping
    if (solver_type .eq. 1) then
       call pf_MGRIT_run(pf, mg_ld, y_0, y_end)
    else if (solver_type .eq. 2) then
       call pf_parareal_run(pf, y_0, dt, Tfin, nsteps, y_end)
    else if (solver_type .eq. 3) then
       call initialize_results(pf)
       if (pf%save_timings > 0) call pf_start_timer(pf, T_TOTAL)
       call pf%levels(1)%ulevel%stepper%do_n_steps(pf, 1, T0, y_0, y_end, dt, nsteps)
       if (pf%save_timings > 0) call pf_stop_timer(pf, T_TOTAL)
       call pf_dump_stats(pf)
    else
       call pf_pfasst_run(pf, y_0, dt, Tfin, nsteps, y_end)
    end if
    !if (pf%rank .eq. pf%comm%nproc-1) call y_end%eprint()
    !call y_end%eprint()

    call mpi_comm_size(pf%comm%comm, nproc, error)
    call mpi_comm_rank(pf%comm%comm, rank, error)

    !>  Wait for everyone to be done
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    
    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst
end program
