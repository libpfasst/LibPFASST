!
! This file is part of LIBPFASST.
!
! Level specification for  zndarray encapsulations for rk steppers

module pf_my_level
  use pfasst
  use pf_mod_zndarray
  use pf_mod_fftpackage
  implicit none

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
    use pf_my_stepper, only: my_stepper_t, as_my_stepper
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    !  Local  variables
    type(pf_fft_t),     pointer :: fft

    class(my_stepper_t), pointer :: stepper_f, stepper_c  ! pointer to levels
    class(pf_zndarray_t), pointer :: f_z,c_z              ! pointer to encaps
    stepper_c => as_my_stepper(c_lev%ulevel%stepper)
    stepper_f => as_my_stepper(f_lev%ulevel%stepper)
    f_z => cast_as_zndarray(f_vec)
    c_z => cast_as_zndarray(c_vec)
    fft => stepper_f%fft_tool

    !  We use the stepper tmp pointers to avoid having to specify dimension
     call f_z%get_array(stepper_f%p_tmp) 
     call c_z%get_array(stepper_c%p_tmp)
    
     call fft%interp(stepper_c%p_tmp,stepper_f%p_tmp)
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use pf_my_stepper, only: my_stepper_t, as_my_stepper

    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags

    class(my_stepper_t), pointer :: stepper_f, stepper_c
    class(pf_zndarray_t), pointer :: f_z,c_z              ! pointer to encaps
    type(pf_fft_t),     pointer :: fft

    stepper_c => as_my_stepper(c_lev%ulevel%stepper)
    stepper_f => as_my_stepper(f_lev%ulevel%stepper)
    f_z => cast_as_zndarray(f_vec)
    c_z => cast_as_zndarray(c_vec)
    fft => stepper_f%fft_tool

    call f_z%get_array(stepper_f%p_tmp) 
    call c_z%get_array(stepper_c%p_tmp)

    call fft%restrict(stepper_f%p_tmp,stepper_c%p_tmp)
    
  end subroutine restrict


end module pf_my_level


