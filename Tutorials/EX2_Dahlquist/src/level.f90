!
! This file is part of LIBPFASST.
!
!> Level specification for Dahlquist example.
!>     u_t = lam1*u + lam2*u
module pf_my_level
  use encap
  use pf_mod_imex_sweeper

  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: my_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type my_level_t

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

    class(scalar_encap), pointer :: y_f, y_c

    !>  Cast the abstract encap as my data type
    y_f => cast_as_scalar(f_vec)
    y_c => cast_as_scalar(c_vec)

    !> Here we use the identity map    
    y_f%y=y_c%y
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(my_level_t), intent(inout)   :: this
    class(pf_level_t), intent(inout)   :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout) :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   )   :: t             !<  time of solution
    integer, intent(in), optional :: flags

    class(scalar_encap), pointer :: y_f, y_c

    !>  Cast the abstract encap as my data type
    y_f => cast_as_scalar(f_vec)
    y_c => cast_as_scalar(c_vec)

    !> Here we use the identity map    
    y_c%y = y_f%y
  end subroutine restrict

end module pf_my_level

