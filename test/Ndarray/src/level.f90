!
! This file is part of LIBPFASST.
!
! Level specification for  ndarray encapsulations

module pf_my_level
  use pf_mod_pfasst
  use pf_mod_ndarray
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
    use pf_my_sweeper, only: my_sweeper_t, as_my_sweeper
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    !  Local  variables
    type(pf_fft_t),     pointer :: fft_f
    type(pf_fft_t),     pointer :: fft_c

    class(my_sweeper_t), pointer :: sweeper_f, sweeper_c  ! pointer to levels
    class(pf_ndarray_t), pointer :: f_z,c_z              ! pointer to encaps
    sweeper_c => as_my_sweeper(c_lev%ulevel%sweeper)
    sweeper_f => as_my_sweeper(f_lev%ulevel%sweeper)
    f_z => cast_as_ndarray(f_vec)
    c_z => cast_as_ndarray(c_vec)
    fft_f => sweeper_f%fft_tool
    fft_c => sweeper_c%fft_tool

    !  We use the sweeper tmp pointers to avoid having to specify dimension
     call f_z%get_array(sweeper_f%p_tmp) 
     call c_z%get_array(sweeper_c%p_tmp)
    
     call fft_c%interp(sweeper_c%p_tmp,fft_f,sweeper_f%p_tmp)
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use pf_my_sweeper, only: my_sweeper_t, as_my_sweeper

    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags

    class(my_sweeper_t), pointer :: sweeper_f, sweeper_c
    class(pf_ndarray_t), pointer :: f_z,c_z              ! pointer to encaps
    type(pf_fft_t),     pointer :: fft

    sweeper_c => as_my_sweeper(c_lev%ulevel%sweeper)
    sweeper_f => as_my_sweeper(f_lev%ulevel%sweeper)
    f_z => cast_as_ndarray(f_vec)
    c_z => cast_as_ndarray(c_vec)
    fft => sweeper_f%fft_tool

    call f_z%get_array(sweeper_f%p_tmp) 
    call c_z%get_array(sweeper_c%p_tmp)

    call fft%restrict(sweeper_f%p_tmp,sweeper_c%p_tmp)
    
  end subroutine restrict


end module pf_my_level


