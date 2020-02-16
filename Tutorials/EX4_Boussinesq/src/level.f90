!
! This file is part of LIBPFASST.
!
! Level specification for 1-D advection/diffusion example: u_t + v*u_x = nu*u_xx

module pf_my_level
  use pf_mod_dtype
  use pf_mod_zndsysarray
  use pf_mod_fftpackage
  implicit none

  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  These are the transfer functions that must be  provided for the level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>  Interpolate from coarse  level to fine
  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use pf_my_sweeper, only: my_sweeper_t, as_my_sweeper
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    type(pf_fft_t),     pointer :: fft
    integer :: nx_f, nx_c
    complex(pfdp),         pointer :: yvec_f(:,:), yvec_c(:,:)

    class(my_sweeper_t), pointer :: sweeper_f, sweeper_c
    sweeper_c => as_my_sweeper(c_lev%ulevel%sweeper)
    sweeper_f => as_my_sweeper(f_lev%ulevel%sweeper)

    fft => sweeper_f%fft_tool
    
    yvec_f => get_array2d(f_vec,1) 
    yvec_c => get_array2d(c_vec,1)

    call fft%interp(yvec_c,yvec_f)
    yvec_f => get_array2d(f_vec,2) 
    yvec_c => get_array2d(c_vec,2)

    call fft%interp(yvec_c,yvec_f)

  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use pf_my_sweeper, only: my_sweeper_t, as_my_sweeper

    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags


    type(pf_fft_t),     pointer :: fft
    complex(pfdp), pointer :: yvec_f(:,:), yvec_c(:,:)  
    integer :: nx_f, nx_c


    class(my_sweeper_t), pointer :: sweeper_f, sweeper_c
    sweeper_c => as_my_sweeper(c_lev%ulevel%sweeper)
    sweeper_f => as_my_sweeper(f_lev%ulevel%sweeper)

    fft => sweeper_f%fft_tool

    yvec_f => get_array2d(f_vec,1) 
    yvec_c => get_array2d(c_vec,1)

    call fft%restrict(yvec_f,yvec_c)
    yvec_f => get_array2d(f_vec,2) 
    yvec_c => get_array2d(c_vec,2)

    call fft%restrict(yvec_f,yvec_c)
    
  end subroutine restrict



end module pf_my_level


