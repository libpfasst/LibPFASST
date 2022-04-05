!
! This file is part of LIBPFASST.
!
!> Level specification  for 1-D advection/diffusion example.
!>     u_t + v*u_x = nu*u_xx
module my_level
  use pf_mod_dtype
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

  !>  Interpolate from coarse  level to fine using FFT
  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use my_sweeper, only: my_sweeper_t, as_my_sweeper
    use my_stepper, only: my_stepper_t, as_my_stepper
    use probin
    
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t                  !  Equation time
    integer, intent(in), optional :: flags                 !  Optional flags (not used here)

    !  Local variables
    integer :: nx_f, nx_c
    integer :: irat       !  Coarsening ratio
    class(my_sweeper_t), pointer :: sweeper_f, sweeper_c  !  fine and coarse sweepers
    class(my_stepper_t), pointer :: stepper_f, stepper_c
    complex(pfdp),       pointer :: wk_f(:),wk_c(:)       !  fine and coarse FFT workspaces
    type(pf_fft_t),      pointer :: fft_f,fft_c           !  fine and coarse FFT packages

    real(pfdp),          pointer :: yvec_f_1d(:),     yvec_c_1d(:)  !  fine and coarse solutions
    real(pfdp),          pointer :: yvec_f_2d(:,:),   yvec_c_2d(:,:)
    real(pfdp),          pointer :: yvec_f_3d(:,:,:), yvec_c_3d(:,:,:)


    if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3)  .or. (solver_type .eq. 4)) then
       stepper_c => as_my_stepper(c_lev%ulevel%stepper)
       stepper_f => as_my_stepper(f_lev%ulevel%stepper)
       fft_c => stepper_c%fft_tool
       fft_f => stepper_f%fft_tool
    else
       sweeper_c => as_my_sweeper(c_lev%ulevel%sweeper)
       sweeper_f => as_my_sweeper(f_lev%ulevel%sweeper)
       fft_c => sweeper_c%fft_tool
       fft_f => sweeper_f%fft_tool
    end if

    if (ndim == 1) then
       yvec_f_1d => get_array1d(f_vec)
       yvec_c_1d => get_array1d(c_vec)

       nx_f = size(yvec_f_1d,1)
       nx_c = size(yvec_c_1d,1)
       irat  = nx_f / nx_c

       !>  If
       if (irat == 1) then !  Identity map
          yvec_f_1d = yvec_c_1d
       elseif (irat == 2) then  !  Use spectral space
          call fft_c%interp(yvec_c_1d,fft_f,yvec_f_1d)
       end if
    else if (ndim == 2) then
       yvec_f_2d => get_array2d(f_vec)
       yvec_c_2d => get_array2d(c_vec)

       nx_f = size(yvec_f_2d,1)
       nx_c = size(yvec_c_2d,1)
       irat  = nx_f / nx_c

       !>  If
       if (irat == 1) then !  Identity map
          yvec_f_2d = yvec_c_2d
       elseif (irat == 2) then  !  Use spectral space
          call fft_c%interp(yvec_c_2d,fft_f,yvec_f_2d)
       end if
    else if (ndim == 3) then
       yvec_f_3d => get_array3d(f_vec)
       yvec_c_3d => get_array3d(c_vec)

       nx_f = size(yvec_f_3d,1)
       nx_c = size(yvec_c_3d,1)
       irat  = nx_f / nx_c

       !>  If
       if (irat == 1) then !  Identity map
          yvec_f_3d = yvec_c_3d
       elseif (irat == 2) then  !  Use spectral space
          call fft_c%interp(yvec_c_3d,fft_f,yvec_f_3d)
       end if
    else

    end if

  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use probin, only:  ndim
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags


    real(pfdp), pointer :: yvec_f_1d(:),     yvec_c_1d(:)  
    real(pfdp), pointer :: yvec_f_2d(:,:),   yvec_c_2d(:,:)
    real(pfdp), pointer :: yvec_f_3d(:,:,:), yvec_c_3d(:,:,:)

    integer :: irat

    !>  Grab the vectors from the encap
    if (ndim == 1) then
       yvec_f_1d => get_array1d(f_vec)
       yvec_c_1d => get_array1d(c_vec)
       irat  = size(yvec_f_1d,1)/size(yvec_c_1d,1)
    else if (ndim == 2) then
       yvec_f_2d => get_array2d(f_vec)
       yvec_c_2d => get_array2d(c_vec)
       irat  = size(yvec_f_2d,1)/size(yvec_c_2d,1)
    else if (ndim == 3) then
       yvec_f_3d => get_array3d(f_vec)
       yvec_c_3d => get_array3d(c_vec)
       irat  = size(yvec_f_3d,1)/size(yvec_c_3d,1)
    else

    end if

    if (ndim == 1) then
       if (irat == 1) then !  Identity map
          yvec_c_1d = yvec_f_1d
       elseif (irat == 2) then !>  Pointwise coarsening
          yvec_c_1d = yvec_f_1d(::irat)
       end if
    else if (ndim == 2) then
       if (irat == 1) then !  Identity map
          yvec_c_2d = yvec_f_2d
       elseif (irat == 2) then !>  Pointwise coarsening
          yvec_c_2d = yvec_f_2d(::irat, ::irat)
       end if
    else if (ndim == 3) then
       if (irat == 1) then !  Identity map
          yvec_c_3d = yvec_f_3d
       elseif (irat == 2) then !>  Pointwise coarsening
          yvec_c_3d = yvec_f_3d(::irat, ::irat, ::irat)
       end if
    else

    end if
  end subroutine restrict

end module my_level

