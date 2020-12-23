module pf_my_level
  use encap
  use pf_mod_imex_sweeper

  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: my_level_t
   integer :: nlevels
   contains
     procedure :: restrict
     procedure :: interpolate
  end type my_level_t

  interface

     subroutine HypreRestrict(y_f, y_c, f_level, c_level) bind(c, name="HypreRestrict")
           use iso_c_binding
           type(c_ptr), value :: y_c, y_f
           integer, value :: f_level, c_level
     end subroutine HypreRestrict

     subroutine HypreProlong(y_f, y_c, f_level, c_level) bind(c, name="HypreProlong")
           use iso_c_binding
           type(c_ptr), value :: y_c, y_f
           integer, value :: f_level, c_level
     end subroutine HypreProlong

  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  These are the transfer functions that must be  provided for the level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>  Interpolate from coarse level to fine
  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(my_level_t), intent(inout) :: this
    real(pfdp),        intent(in   ) :: t
    class(pf_level_t), intent(inout) :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t), intent(inout) :: f_vec, c_vec  !  fine and coarse vectors
    integer, intent(in), optional :: flags
    integer :: f_level_index, c_level_index
    integer :: restrict_flag

    class(hypre_vector_encap), pointer :: y_f, y_c

    !>  Cast the abstract encap as my data type
    y_f => cast_as_hypre_vector(f_vec)
    y_c => cast_as_hypre_vector(c_vec)
 
    f_level_index = f_lev%index
    c_level_index = c_lev%index

    !> Here we use the identity map    
    !call y_f%copy(y_c)
    call HypreProlong(y_f%c_hypre_vector_ptr, &
                      y_c%c_hypre_vector_ptr, &
                      f_level_index, c_level_index) 
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),  intent(inout) :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp), intent(in) :: t             !<  time of solution
    integer, intent(in), optional :: flags
    integer :: f_level_index, c_level_index
    integer :: restrict_flag

    class(hypre_vector_encap), pointer :: y_f, y_c

    !>  Cast the abstract encap as my data type
    y_f => cast_as_hypre_vector(f_vec)
    y_c => cast_as_hypre_vector(c_vec)

    !f_level_index = abs(f_lev%index - this%nlevels)
    !c_level_index = abs(c_lev%index - this%nlevels)

    f_level_index = f_lev%index
    c_level_index = c_lev%index

    !> Here we use the identity map    
    !call y_c%copy(y_f)
    call HypreRestrict(y_f%c_hypre_vector_ptr, &
                       y_c%c_hypre_vector_ptr, &
                       f_level_index, c_level_index)
  end subroutine restrict

  function cast_as_my_level_t(level_polymorph) result(my_level_t_obj)
    class(pf_user_level_t), intent(in), target :: level_polymorph
    type(my_level_t), pointer :: my_level_t_obj

    select type(level_polymorph)
    type is (my_level_t)
       my_level_t_obj => level_polymorph
    end select
  end function cast_as_my_level_t

end module pf_my_level

