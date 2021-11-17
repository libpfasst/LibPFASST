!-------------------------------------------------------------------------------
! Copyright (c) 2017, Brandon Krull.  All rights reserved.
!--------------------------------------------------------------------------------
! MODULE: pf_my_level
! !> @author
!> Brandon Krull, Berkeley Lab
!
! Description:
!> This module extend the abstract level
module pf_my_level
  use pf_mod_dtype
  use pf_mod_imk
  use mod_zmkpair
  use utils

  implicit none

  external :: zgemm

  type, extends(pf_user_level_t) :: imk_context
   contains
     procedure :: restrict => restrict
     procedure :: interpolate => interpolate
  end type imk_context


contains


 subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
   class(imk_context), intent(inout) :: this
   class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
   class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
   real(pfdp),        intent(in   ) :: t
   integer, intent(in), optional :: flags

   class(zmkpair), pointer :: f, g
   f => cast_as_zmkpair(f_vec)
   g => cast_as_zmkpair(c_vec)

   g%array = f%array
   g%y = f%y
 end subroutine restrict

  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(imk_context), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
   real(pfdp),        intent(in   ) :: t
   integer, intent(in), optional :: flags

   class(zmkpair), pointer :: f, g
   f => cast_as_zmkpair(f_vec)
   g => cast_as_zmkpair(c_vec)

   f%array = g%array
   f%y = g%y
 end subroutine interpolate

end module pf_my_level
