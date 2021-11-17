module pfasst_hypre
  use encap
  use pfasst
  use pf_my_level
  use probin
  use pf_my_sweeper
  use pf_my_stepper
  use pf_mod_mgrit
  implicit none
contains

  subroutine PfasstHypreInit(pf, mg_ld, lev_shape, space_color, time_color, spacial_coarsen_flag) 
    type(pf_pfasst_t), intent(inout) :: pf
    type(mgrit_level_data), allocatable, intent(inout) :: mg_ld(:)
    integer, allocatable, intent(inout) :: lev_shape(:,:)
    integer, intent(in) :: space_color, time_color
    integer, intent(in) :: spacial_coarsen_flag
    
    type(pf_comm_t) :: comm
    type(my_sweeper_t) :: sw_finest, sw_lev
    type(my_stepper_t) :: st_finest, st_lev
    type(my_level_t) :: my_lev

    integer :: l, l_finest   !  loop variable over levels
    integer :: n, m
    integer :: level_index

    integer :: nproc, rank, error
    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1
    integer :: n_init, refine_factor, FComp_setup_flag
    logical :: setup_start_coarse_flag

    if ((solver_type .eq. 1) .and. (rk_order .eq. 1)) then
       FComp_setup_flag = 0
    else
       FComp_setup_flag = 1
    end if
    

    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)


    allocate(lev_shape(pf%nlevels,10))
    !> Loop over levels and set some level specific parameters
    do l = 1,pf%nlevels
       lev_shape(l,1) = num_grid_points
       lev_shape(l,2) = space_color
       lev_shape(l,3) = space_dim
       lev_shape(l,4) = max_space_v_cycles
       lev_shape(l,10) = spacial_coarsen_flag
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data constructor
       allocate(hypre_vector_factory::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3)) then
          allocate(my_stepper_t::pf%levels(l)%ulevel%stepper)
       else
          allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)
       end if

       call pf_level_set_size(pf, l, lev_shape(l,:), 0)

       if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3)) then
          pf%levels(l)%ulevel%stepper%order = rk_order
          if (solver_type .eq. 2) then
             pf%levels(l)%ulevel%stepper%nsteps = nsteps_rk(l)
          end if
       end if
    end do

    l_finest = pf%nlevels
    if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3)) then
       if (spacial_coarsen_flag .eq. 1) then
          st_finest = cast_as_my_stepper_t(pf%levels(l_finest)%ulevel%stepper)
          call HypreSolverInit(st_finest%c_hypre_solver_ptr, &
                               l_finest, &
                               num_grid_points, &
                               space_color, &
                               space_dim, &
                               max_space_v_cycles, &
                               pf%nlevels, &
                               spacial_coarsen_flag)
       end if
       do l = pf%nlevels,1,-1 
          if (spacial_coarsen_flag .eq. 1) then
             st_lev = cast_as_my_stepper_t(pf%levels(l_finest)%ulevel%stepper)
          else
             st_lev = cast_as_my_stepper_t(pf%levels(l)%ulevel%stepper)
             call HypreSolverInit(st_lev%c_hypre_solver_ptr, &
                                  l, &
                                  num_grid_points, &
                                  space_color, &
                                  space_dim, &
                                  max_space_v_cycles, &
                                  pf%nlevels, &
                                  spacial_coarsen_flag)
          end if
 
          n = st_lev%get_nrows(l)
          ilower0 = st_lev%get_extent(l, 0)
          ilower1 = st_lev%get_extent(l, 1)
          iupper0 = st_lev%get_extent(l, 2)
          iupper1 = st_lev%get_extent(l, 3)
   
          lev_shape(l,5) = n
          lev_shape(l,6) = ilower0
          lev_shape(l,7) = ilower1
          lev_shape(l,8) = iupper0
          lev_shape(l,9) = iupper1
   
          call pf_level_set_size(pf, l, lev_shape(l,:), n)
       end do
    else
       if (spacial_coarsen_flag .eq. 1) then
          sw_finest = cast_as_my_sweeper_t(pf%levels(l_finest)%ulevel%sweeper)
          call HypreSolverInit(sw_finest%c_hypre_solver_ptr, &
                               l_finest, &
                               num_grid_points, &
                               space_color, &
                               space_dim, &
                               max_space_v_cycles, &
                               pf%nlevels, &
                               spacial_coarsen_flag)
       end if

       do l = pf%nlevels,1,-1
          if (spacial_coarsen_flag .eq. 1) then
             sw_lev = cast_as_my_sweeper_t(pf%levels(l_finest)%ulevel%sweeper)
          else
             sw_lev = cast_as_my_sweeper_t(pf%levels(l)%ulevel%sweeper)
             call HypreSolverInit(sw_lev%c_hypre_solver_ptr, &
                                  l, &
                                  num_grid_points, &
                                  space_color, &
                                  space_dim, &
                                  max_space_v_cycles, &
                                  pf%nlevels, &
                                  spacial_coarsen_flag)
          end if

          n = sw_lev%get_nrows(l)
          ilower0 = sw_lev%get_extent(l, 0)
          ilower1 = sw_lev%get_extent(l, 1)
          iupper0 = sw_lev%get_extent(l, 2)
          iupper1 = sw_lev%get_extent(l, 3)

          lev_shape(l,5) = n
          lev_shape(l,6) = ilower0
          lev_shape(l,7) = ilower1
          lev_shape(l,8) = iupper0
          lev_shape(l,9) = iupper1

          call pf_level_set_size(pf, l, lev_shape(l,:), n)
       end do
    end if

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3)) then
       call stepper_hypre_set_level_data(pf)
    else
       call sweeper_hypre_set_level_data(pf)
    end if

    if (solver_type .eq. 1) then
       T0 = 0.0_pfdp
       setup_start_coarse_flag = .false.
       if (setup_start_coarse_flag .eqv. .true.) then
          n_init = max(1, mgrit_n_init/pf%comm%nproc)
       else
          n_init = mgrit_n_init / pf%comm%nproc
          !n_init = pf%state%nsteps / pf%comm%nproc
       end if
       refine_factor = mgrit_refine_factor
       call mgrit_initialize(pf, mg_ld, T0, Tfin, n_init, refine_factor, FAS_flag, FCF_flag, setup_start_coarse_flag)
    end if

    if (FComp_setup_flag .eq. 0) then
       if (solver_type .eq. 1) then
          do l = pf%nlevels,1,-1
             st_lev = cast_as_my_stepper_t(pf%levels(l)%ulevel%stepper)
             call HypreImplicitSolverInit(st_lev%c_hypre_solver_ptr, &
                                          l, &
                                          num_grid_points, &
                                          space_color, &
                                          space_dim, &
                                          max_space_v_cycles, &
                                          pf%nlevels, &
                                          mg_ld(l)%dt)
          end do
       else
          do l = pf%nlevels,1,-1
             sw_lev = cast_as_my_sweeper_t(pf%levels(l)%ulevel%sweeper)
             call HypreImplicitSolverInit(sw_lev%c_hypre_solver_ptr, &
                                          l, &
                                          num_grid_points, &
                                          space_color, &
                                          space_dim, &
                                          max_space_v_cycles, &
                                          pf%nlevels, &
                                          mg_ld(l)%dt)
          end do
       end if
    end if
  end subroutine PfasstHypreInit

end module pfasst_hypre
