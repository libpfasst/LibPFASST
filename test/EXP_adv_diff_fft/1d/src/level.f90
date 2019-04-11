!
! This file is part of LIBPFASST.
!
! Level specification for 1-D advection/diffusion example: u_t + v*u_x = nu*u_xx

module pf_my_level
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_exp
  use phi_mod
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
    use pf_my_sweeper, only: ad_sweeper_t, as_ad_sweeper
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags


    integer :: nvarF, nvarG, xrat
    class(ad_sweeper_t), pointer :: sweeper_f, sweeper_c
    real(pfdp),         pointer :: yvec_f(:), yvec_c(:)
    complex(pfdp),         pointer ::  wk_f(:),wk_c(:)    
    type(pf_fft_t),     pointer :: fft_f,fft_c

    sweeper_c => as_ad_sweeper(c_lev%ulevel%sweeper)
    sweeper_f => as_ad_sweeper(f_lev%ulevel%sweeper)
    fft_c => sweeper_c%fft_tool
    fft_f => sweeper_f%fft_tool  

    yvec_f => get_array1d(f_vec) 
    yvec_c => get_array1d(c_vec)

    nvarF = size(yvec_f)
    nvarG = size(yvec_c)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       yvec_f = yvec_c
       return
    endif
    call fft_f%get_wk_ptr(wk_f)
    call fft_c%get_wk_ptr(wk_c)
    wk_c=yvec_c
    call fft_c%fftf()    
    wk_f = 0.0d0
    wk_f(1:nvarG/2) = wk_c(1:nvarG/2)
    wk_f(nvarF-nvarG/2+2:nvarF) = wk_c(nvarG/2+2:nvarG)

    wk_f=wk_f*2.0_pfdp

    call fft_f%fftb()
    yvec_f=real(wk_f)
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use pf_my_sweeper, only: ad_sweeper_t

    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags


    real(pfdp), pointer :: yvec_f(:), yvec_c(:)  

    integer :: irat

    yvec_f => get_array1d(f_vec)
    yvec_c => get_array1d(c_vec)

    irat  = size(yvec_f)/size(yvec_c)

    yvec_c = yvec_f(::irat)
  end subroutine restrict



end module pf_my_level


