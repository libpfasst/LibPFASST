!
! This file is part of LIBPFASST.
!
!> Level specification  for 2-D advection/diffusion example.
!>     u_t + a*u_x+b*u_y = nu*(u_xx+u_yy)
module pf_my_level
  use pf_mod_dtype
  use pf_mod_ndarray
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

  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use pf_my_sweeper, only: ad_sweeper_t, as_ad_sweeper    
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    integer      :: nx_f, nx_c, ny_f, ny_c
    integer      :: irat,jrat,i,j,ii,jj
    class(ad_sweeper_t), pointer :: sweeper_f, sweeper_c
    real(pfdp),         pointer :: yvec_f(:,:), yvec_c(:,:)
    complex(pfdp),         pointer ::  wk_f(:,:),wk_c(:,:)    
    type(pf_fft_t),     pointer :: fft_f,fft_c

    sweeper_c => as_ad_sweeper(c_lev%ulevel%sweeper)
    sweeper_f => as_ad_sweeper(f_lev%ulevel%sweeper)

    yvec_f => get_array2d(f_vec) 
    yvec_c => get_array2d(c_vec)

    nx_f=size(yvec_f,1)
    ny_f=size(yvec_f,2)
    nx_c=size(yvec_c,1)
    ny_c=size(yvec_c,2)
    irat  = nx_f/nx_c
    jrat  = ny_f/ny_c

    if (irat == 1 .and. jrat==1) then
       yvec_f = yvec_c
       return
    endif

    fft_c => sweeper_c%fft_tool
    fft_f => sweeper_f%fft_tool    

    call fft_c%interp(yvec_c,fft_f,yvec_f)
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags


    real(pfdp), pointer :: yvec_f(:,:), yvec_c(:,:)  
    integer      :: irat,jrat
    yvec_f => get_array2d(f_vec)
    yvec_c => get_array2d(c_vec)

    irat  = size(yvec_f,1)/size(yvec_c,1)
    jrat  = size(yvec_f,2)/size(yvec_c,2)

    if (irat == 1 .and. jrat==1) then
       yvec_c = yvec_f
    else
       yvec_c = yvec_f(::irat,::jrat)       
    endif


  end subroutine restrict

end module pf_my_level  
