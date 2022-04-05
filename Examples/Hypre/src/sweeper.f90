module pf_my_sweeper
  use encap
  use pf_mod_imex_sweeper
  use pfasst
  use pf_my_level
  use probin
  implicit none

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: my_sweeper_t
   type(c_ptr) :: c_hypre_solver_ptr = c_null_ptr
   contains

     procedure :: f_eval    !  Computes the explicit rhs terms
     procedure :: f_comp    !  Does implicit solves 

     procedure :: initialize !  Overwrites imex sweeper initialize
     procedure :: destroy    !  Overwrites imex sweeper destroy

     procedure :: get_nrows
     procedure :: get_extent

  end type my_sweeper_t

  interface

     subroutine HypreSolverInit(hypre_solver_ptr, level_index, nx, comm_color, space_dim, max_iter, num_levels, spatial_coarsen_flag) bind(c, name="HypreSolverInit")
        use iso_c_binding
        type(c_ptr) :: hypre_solver_ptr
        integer, value :: nx, level_index, comm_color, space_dim, max_iter, num_levels, spatial_coarsen_flag
     end subroutine HypreSolverInit

     subroutine HypreImplicitSolverInit(hypre_solver_ptr, level_index, nx, comm_color, space_dim, max_iter, num_levels, dtq) bind(c, name="HypreImplicitSolverInit")
        use iso_c_binding
        type(c_ptr) :: hypre_solver_ptr
        integer, value :: nx, level_index, comm_color, space_dim, max_iter, num_levels
        real(c_double), value :: dtq
     end subroutine HypreImplicitSolverInit
   
     subroutine HypreSolverDestroy(hypre_solver, level_index) bind(c, name="HypreSolverDestroy")
        use iso_c_binding
        type(c_ptr), value :: hypre_solver
        integer, value :: level_index
     end subroutine HypreSolverDestroy
     
     subroutine HypreSolverFEval(hypre_solver, y, t, level_index, f, piece) bind(c, name="HypreSolverFEval")
        use iso_c_binding
        type(c_ptr), value :: hypre_solver
        type(c_ptr), value :: y
        real(c_double), value :: t
        integer, value :: level_index
        type(c_ptr), value :: f
        integer, value :: piece
     end subroutine HypreSolverFEval
   
     subroutine HypreSolverFComp(hypre_solver, y, t, dtq, rhs, level_index, f, piece) bind(c, name="HypreSolverFComp")
        use iso_c_binding
        type(c_ptr), value :: hypre_solver
        type(c_ptr), value :: y
        real(c_double), value :: t
        real(c_double), value :: dtq
        type(c_ptr), value :: rhs
        integer, value :: level_index
        type(c_ptr), value :: f
        integer, value :: piece
     end subroutine HypreSolverFComp

     subroutine HypreVectorSetSinInitCond(y_0) bind(c, name="HypreVectorSetSinInitCond")
        use iso_c_binding
        type(c_ptr), value :: y_0
     end subroutine HypreVectorSetSinInitCond

     function HypreSolverGetNumLevels(hypre_solver) result(num_levels) bind(c, name="HypreSolverGetNumLevels")
        use iso_c_binding
        type(c_ptr), value :: hypre_solver
        integer :: num_levels
     end function HypreSolverGetNumLevels

     subroutine HypreSolverSetLevelData(y, x, level_index) bind(c, name="HypreSolverSetLevelData")
        use iso_c_binding
        type(c_ptr) :: y
        type(c_ptr), value :: x
        integer, value :: level_index
     end subroutine HypreSolverSetLevelData

     function HypreSolverGetNumRowsLevel(hypre_solver, level_index) result(nrows) bind(c, name="HypreSolverGetNumRowsLevel")
        use iso_c_binding
        type(c_ptr), value :: hypre_solver
        integer, value :: level_index
        integer :: nrows
     end function

     function HypreSolverGetExtentLevel(hypre_solver, level_index, i) result(extent) bind(c, name="HypreSolverGetExtentLevel")
        use iso_c_binding
        type(c_ptr), value :: hypre_solver
        integer, value :: level_index, i
        integer :: extent
     end function

   end interface

contains


  !>  Routine to set up sweeper variables and operators
  subroutine initialize(this, pf, level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    integer, intent(in) :: level_index
    integer :: nx, comm_color, space_dim, max_space_v_cycles, spatial_coarsen_flag

    if (imex_stat .eq. 0 ) then
       this%explicit=.TRUE.
       this%implicit=.FALSE.
    elseif (imex_stat .eq. 1 ) then
       this%implicit=.TRUE.
       this%explicit=.FALSE.
    else
       this%implicit=.TRUE.
       this%explicit=.TRUE.
    end if
 
    !>  Call the imex sweeper initialization
    call this%imex_initialize(pf,level_index)

    ! Space variables
    nx = pf%levels(level_index)%lev_shape(1)
    comm_color = pf%levels(level_index)%lev_shape(2)
    space_dim = pf%levels(level_index)%lev_shape(3)
    max_space_v_cycles = pf%levels(level_index)%lev_shape(4)
    spatial_coarsen_flag = pf%levels(level_index)%lev_shape(10)

    !> Call the Hypre solver initialization
    call HypreSolverInit(this%c_hypre_solver_ptr, &
                         level_index, &
                         nx, &
                         comm_color, &
                         space_dim, &
                         max_space_v_cycles, &
                         pf%nlevels, &
                         spatial_coarsen_flag)
  end subroutine initialize

  !>  destroy the sweeper type
  subroutine destroy(this, pf,level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    integer, intent(in) :: level_index

    !> Destroy Hypre solver variables
    call HypreSolverDestroy(this%c_hypre_solver_ptr, level_index)

    !>  Call the imex sweeper destroy
    call this%imex_destroy(pf,level_index)
  end subroutine destroy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    real(pfdp) :: val
    
    class(hypre_vector_encap), pointer :: y_encap, f_encap
    
    y_encap => cast_as_hypre_vector(y)
    f_encap => cast_as_hypre_vector(f)

    !> Call Hypre to do the f_eval
    call HypreSolverFEval(this%c_hypre_solver_ptr, y_encap%c_hypre_vector_ptr, t, level_index, f_encap%c_hypre_vector_ptr, piece)
  end subroutine f_eval

  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f, piece)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    class(hypre_vector_encap), pointer :: y_encap, f_encap, rhs_encap
    real(pfdp) :: val

    y_encap => cast_as_hypre_vector(y)
    f_encap => cast_as_hypre_vector(f)
    rhs_encap => cast_as_hypre_vector(rhs)

    !> Call Hypre to do the f_comp
    call HypreSolverFComp(this%c_hypre_solver_ptr, &
                          y_encap%c_hypre_vector_ptr, &
                          t, &
                          dtq, &
                          rhs_encap%c_hypre_vector_ptr, &
                          level_index, &
                          f_encap%c_hypre_vector_ptr, &
                          piece)
  end subroutine f_comp

  function get_nrows(this, level_index) result(nrows)
     class(my_sweeper_t), intent(inout) :: this
     integer, intent(in) :: level_index
     integer :: nrows
     nrows = HypreSolverGetNumRowsLevel(this%c_hypre_solver_ptr, level_index)
  end function

  function get_extent(this, level_index, i) result(extent)
     class(my_sweeper_t), intent(inout) :: this
     integer, intent(in) :: level_index, i
     integer :: extent
     extent = HypreSolverGetExtentLevel(this%c_hypre_solver_ptr, level_index, i)
  end function

  !> Routine to set initial condition.
  subroutine initial(y_0)
    type(hypre_vector_encap), intent(inout) :: y_0
    call HypreVectorSetSinInitCond(y_0%c_hypre_vector_ptr)
  end subroutine initial

  function cast_as_my_sweeper_t(pf_sweeper_t_polymorph) result(my_sweeper_t_obj)
    class(pf_sweeper_t), intent(in), target :: pf_sweeper_t_polymorph
    type(my_sweeper_t), pointer :: my_sweeper_t_obj

    select type(pf_sweeper_t_polymorph)
    type is (my_sweeper_t)
       my_sweeper_t_obj => pf_sweeper_t_polymorph
    end select
  end function cast_as_my_sweeper_t
end module pf_my_sweeper
