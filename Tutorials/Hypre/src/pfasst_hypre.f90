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

  subroutine PfasstHypreInit(pf, mg_ld, lev_shape, space_color, time_color) 
    type(pf_pfasst_t), intent(inout) :: pf
    integer, allocatable, intent(inout) :: lev_shape(:,:)
    integer, intent(in) :: space_color, time_color
    type(mgrit_level_data), allocatable, intent(inout) :: mg_ld(:)
    
    type(pf_comm_t) :: comm
    type(my_sweeper_t) :: s_finest, s
    type(my_level_t) :: my_lev

    integer :: l, l_finest   !  loop variable over levels
    integer :: n, m
    integer :: level_index

    integer :: nproc, rank, error
    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1

    integer :: n_coarse, refine_factor
    real(pfdp) :: T0
    logical :: FAS_flag, FCF_flag

    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

    allocate(lev_shape(pf%nlevels,9))

    !> Loop over levels and set some level specific parameters
    do l = 1,pf%nlevels
       lev_shape(l,1) = num_grid_points
       lev_shape(l,2) = space_color
       lev_shape(l,3) = space_dim
       lev_shape(l,4) = max_space_v_cycles
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data constructor
       allocate(hypre_vector_factory::pf%levels(l)%ulevel%factory)

       if (use_mgrit .eqv. .true.) then
          allocate(my_stepper_t::pf%levels(l)%ulevel%stepper)
       else
          allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)
       end if

       call pf_level_set_size(pf, l, lev_shape(l,:), n)

       if (use_mgrit .eqv. .true.) then
          pf%levels(l)%ulevel%stepper%order = rk_order
          pf%levels(l)%ulevel%stepper%nsteps = nsteps_rk(l)
       end if
    end do

    l_finest = pf%nlevels
    s_finest = cast_as_my_sweeper_t(pf%levels(l_finest)%ulevel%sweeper)
    call HypreSolverInit(s_finest%c_hypre_solver_ptr, &
                         l_finest, &
                         num_grid_points, &
                         space_color, &
                         space_dim, &
                         max_space_v_cycles, &
                         pf%nlevels)
    do l = 1,pf%nlevels
       n = s_finest%get_nrows(l)
       call pf_level_set_size(pf, l, lev_shape(l,:), n)
    end do

    do l = 1,pf%nlevels
       n = s_finest%get_nrows(l)
       ilower0 = s_finest%get_extent(l, 0)
       ilower1 = s_finest%get_extent(l, 1)
       iupper0 = s_finest%get_extent(l, 2)
       iupper1 = s_finest%get_extent(l, 3)

       lev_shape(l,5) = n
       lev_shape(l,6) = ilower0
       lev_shape(l,7) = ilower1
       lev_shape(l,8) = iupper0
       lev_shape(l,9) = iupper1

       call pf_level_set_size(pf, l, lev_shape(l,:), n)
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)
    call hypre_set_level_data(pf)

    if (use_mgrit .eqv. .true.) then
       FAS_flag = .true.
       FCF_flag = .false.
       T0 = 0.0_pfdp
       n_coarse = max(1, mgrit_n_coarse/pf%comm%nproc)
       refine_factor = mgrit_refine_factor
       call mgrit_initialize(pf, mg_ld, T0, Tfin, n_coarse, refine_factor, FAS_flag, FCF_flag)
    end if
  end subroutine PfasstHypreInit

end module pfasst_hypre
