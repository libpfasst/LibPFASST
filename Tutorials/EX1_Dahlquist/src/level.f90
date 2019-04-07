!
! This file is part of LIBPFASST.
!
!
!> Level specification for Dahlquist example.
!>     u_t = lam1*u + lam2*u
module pf_my_level
  use pf_mod_dtype
  use pf_mod_ndarray
 

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
  !>  Interpolate from coarse  level to fine
  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags


    !> Here we use the identity map
    call f_vec%copy(c_vec)

  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags

    !> Here we use the identity map
    call c_vec%copy(f_vec)    
  end subroutine restrict


end module pf_my_level
