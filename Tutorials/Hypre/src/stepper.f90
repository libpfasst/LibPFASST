module pf_my_stepper
  use encap
  use pf_mod_rkstepper
  use pfasst
  use pf_my_level
  use probin
  use pf_my_sweeper
  implicit none

  !>  extend the stepper type with stuff we need to compute rhs
  type, extends(pf_ark_stepper_t) :: my_stepper_t
   type(c_ptr) :: c_hypre_solver_ptr
   contains

     procedure :: f_eval => f_eval_stepper    !  Computes the explicit rhs terms
     procedure :: f_comp => f_comp_stepper    !  Does implicit solves 

     procedure :: initialize => initialize_stepper !  Overwrites stepper initialize_stepper
     procedure :: destroy => destroy_stepper    !  Overwrites stepper destroy_stepper

  end type my_stepper_t
contains


  !>  Routine to set up stepper variables and operators
  subroutine initialize_stepper(this, pf,level_index)
    class(my_stepper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    integer, intent(in) :: level_index
    integer :: nx, comm_color, space_dim, max_space_v_cycles
 
    !>  Call the stepper initialization
    call this%ark_initialize(pf,level_index)

    this%implicit=.TRUE.
    this%explicit=.FALSE.

    ! Space variables
    nx = pf%levels(level_index)%lev_shape(1)
    comm_color = pf%levels(level_index)%lev_shape(2)
    space_dim = pf%levels(level_index)%lev_shape(3)
    max_space_v_cycles = pf%levels(level_index)%lev_shape(4)

    !> Call the Hypre solver initialization
    call HypreSolverInit(this%c_hypre_solver_ptr, level_index, nx, comm_color, space_dim, max_space_v_cycles, pf%nlevels)
  end subroutine initialize_stepper

  !>  destroy_stepper the stepper type
  subroutine destroy_stepper(this, pf,level_index)
    class(my_stepper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    integer, intent(in) :: level_index

    !> Destroy Hypre solver variables
    call HypreSolverDestroy(this%c_hypre_solver_ptr, level_index)

    !>  Call the stepper destroy_stepper
    call this%ark_destroy(pf,level_index)
  end subroutine destroy_stepper

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval_stepper(this, y, t, level_index, f, piece)
    class(my_stepper_t), intent(inout) :: this
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
  end subroutine f_eval_stepper

  ! Solve for y and return f2 also.
  subroutine f_comp_stepper(this, y, t, dtq, rhs, level_index, f, piece)
    class(my_stepper_t), intent(inout) :: this
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
  end subroutine f_comp_stepper

  function cast_as_my_stepper_t(pf_stepper_t_polymorph) result(my_stepper_t_obj)
    class(pf_stepper_t), intent(in), target :: pf_stepper_t_polymorph
    type(my_stepper_t), pointer :: my_stepper_t_obj

    select type(pf_stepper_t_polymorph)
    type is (my_stepper_t)
       my_stepper_t_obj => pf_stepper_t_polymorph
    end select
  end function cast_as_my_stepper_t
end module pf_my_stepper
