!
! This file is part of LIBPFASST.
!
!> Level specification  for the K-S example.

module my_level
  use pf_mod_dtype
  use pf_mod_ndarray_oc
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

  !>  Interpolate from coarse  level to fine using FFT
  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use my_sweeper, only: my_sweeper_t, as_my_sweeper
    
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t                  !  Equation time
    integer, intent(in), optional :: flags                 !  Optional flags (not used here)

    !  Local variables
    integer :: nx_f, nx_c
    integer :: irat       !  Coarsening ratio
    class(my_sweeper_t), pointer :: sweeper_f, sweeper_c  !  fine and coarse sweepers
    real(pfdp),          pointer :: yvec_f(:), yvec_c(:)  !  fine and coarse solutions
    complex(pfdp),       pointer :: wk_f(:),wk_c(:)       !  fine and coarse FFT workspaces
    type(pf_fft_t),      pointer :: fft_f,fft_c           !  fine and coarse FFT packages
    integer :: which

    which = 1
    if(present(flags)) which = flags
    
    sweeper_c => as_my_sweeper(c_lev%ulevel%sweeper)
    sweeper_f => as_my_sweeper(f_lev%ulevel%sweeper)
    fft_c => sweeper_c%fft_tool
    fft_f => sweeper_f%fft_tool    

    yvec_f => get_array1d_oc(f_vec,which) 
    yvec_c => get_array1d_oc(c_vec,which)

    nx_f = size(yvec_f)
    nx_c = size(yvec_c)
    irat  = nx_f / nx_c

    !>  If 
    if (irat == 1) then !  Identity map
       yvec_f = yvec_c   
       return
    elseif (irat == 2) then  !  Use spectral space
       call fft_c%interp(yvec_c,fft_f,yvec_f)
    end if

  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags


    real(pfdp), pointer :: yvec_f(:), yvec_c(:)  

    integer :: irat, which
    
    which = 1
    if(present(flags)) which = flags

    !>  Grab the vectors from the encap
    yvec_f => get_array1d_oc(f_vec,which)
    yvec_c => get_array1d_oc(c_vec,which)
    
    irat  = size(yvec_f)/size(yvec_c)

    !>  Pointwise coarsening
    yvec_c = yvec_f(::irat)
  end subroutine restrict

end module my_level

