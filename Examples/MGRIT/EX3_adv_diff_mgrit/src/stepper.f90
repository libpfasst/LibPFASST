!
! This file is part of LIBPFASST.
!
!
!> Stepper and RHS routines for advection/diffusion example.
!>     u_t + v*u_x = nu*u_xx
module my_stepper
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_rkstepper
  use pf_mod_fftpackage
  use pf_mod_solutions
  implicit none

  !>  extend the ark stepper type with stuff we need to compute rhs
  type, extends(pf_ark_stepper_t) :: my_stepper_t
     integer ::     nx   !  Grid size

     !>  FFT and Spectral derivatives
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: opE_1d(:), opE_2d(:,:), opE_3d(:,:,:) ! Explicit spectral operator
     complex(pfdp), allocatable :: opI_1d(:), opI_2d(:,:), opI_3d(:,:,:) ! Implicit spectral operator
     
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves
     procedure :: initialize
     procedure :: destroy

  end type my_stepper_t

contains

  !>  Helper function to return stepper pointer
  function as_my_stepper(stepper) result(r)
    class(pf_stepper_t), intent(inout), target :: stepper
    class(my_stepper_t), pointer :: r
    select type(stepper)
    type is (my_stepper_t)
       r => stepper
    class default
       stop
    end select
  end function as_my_stepper

  subroutine initialize(this, pf,level_index)
    use probin, only:  ark_stat, v, nu, ndim
    class(my_stepper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index

    complex(pfdp), allocatable :: lap_1d(:), lap_2d(:,:), lap_3d(:,:,:) ! Lapclacian operator
    complex(pfdp), allocatable :: ddx_1d(:), ddx_2d(:,:), ddx_3d(:,:,:) ! First derivative operator

    integer     :: nx, ny, nz

    !>  Set variables for explicit and implicit parts
    if (ark_stat .eq. 0 ) then
       this%explicit=.TRUE.
       this%implicit=.FALSE.
    elseif (ark_stat .eq. 1 ) then
       this%implicit=.TRUE.
       this%explicit=.FALSE.
    else
       this%implicit=.TRUE.
       this%explicit=.TRUE.
    end if

    !  Call the ark stepper initialize
    call this%ark_initialize(pf,level_index)

    nx=pf%levels(level_index)%lev_shape(1)  !  local convenient grid size
    ny=pf%levels(level_index)%lev_shape(2)
    nz=pf%levels(level_index)%lev_shape(3)

    !>  Set up the FFT
    allocate(this%fft_tool)

    !>  Choose the explicit and implicit operators depending on ark_stat
    if (ndim == 1) then
       call this%fft_tool%fft_setup([nx], ndim)

       allocate(lap_1d(nx))
       allocate(ddx_1d(nx))

       call this%fft_tool%make_lap(lap_1d)
       call this%fft_tool%make_deriv(ddx_1d)

       !> Allocate operators for implicit and explicit parts
       allocate(this%opE_1d(nx))
       allocate(this%opI_1d(nx))

       select case (ark_stat)
       case (0)  ! Fully Explicit
          this%opE_1d = -v(1)*ddx_1d + nu*lap_1d
          this%opI_1d = 0.0_pfdp
       case (1)  ! Fully Implicit
          this%opE_1d = 0.0_pfdp
          this%opI_1d = -v(1)*ddx_1d + nu * lap_1d
       case (2)  ! IMEX
          this%opE_1d = -v(1)*ddx_1d
          this%opI_1d =  nu*lap_1d
       case DEFAULT
          print *,'Bad case for ark_stat in f_eval ', ark_stat
          call exit(0)
       end select
       
       deallocate(lap_1d)
       deallocate(ddx_1d)
    else if (ndim == 2) then
       call this%fft_tool%fft_setup([nx,ny], ndim)

       allocate(lap_2d(nx,ny))
       allocate(ddx_2d(nx,ny))

       call this%fft_tool%make_lap(lap_2d)
       call this%fft_tool%make_deriv(ddx_2d,1)

       !> Allocate operators for implicit and explicit parts
       allocate(this%opE_2d(nx, ny))
       allocate(this%opI_2d(nx, ny))

       select case (ark_stat)
       case (0)  ! Fully Explicit
          this%opE_2d = -v(1)*ddx_2d + nu*lap_2d
          this%opI_2d = 0.0_pfdp
       case (1)  ! Fully Implicit
          this%opE_2d = 0.0_pfdp
          this%opI_2d = -v(1)*ddx_2d + nu * lap_2d
       case (2)  ! IMEX
          this%opE_2d = -v(1)*ddx_2d
          this%opI_2d =  nu*lap_2d
       case DEFAULT
          print *,'Bad case for ark_stat in f_eval ', ark_stat
          call exit(0)
       end select
 
       deallocate(lap_2d)
       deallocate(ddx_2d)
    else if (ndim == 3) then
       call this%fft_tool%fft_setup([nx,ny,nz], ndim)

       allocate(lap_3d(nx,ny,nz))
       allocate(ddx_3d(nx,ny,nz))

       call this%fft_tool%make_lap(lap_3d)
       call this%fft_tool%make_deriv(ddx_3d,1)

       !> Allocate operators for implicit and explicit parts
       allocate(this%opE_3d(nx, ny, nz))
       allocate(this%opI_3d(nx, ny, nz))

       select case (ark_stat)
       case (0)  ! Fully Explicit
          this%opE_3d = -v(1)*ddx_3d + nu*lap_3d
          this%opI_3d = 0.0_pfdp
       case (1)  ! Fully Implicit
          this%opE_3d = 0.0_pfdp
          this%opI_3d = -v(1)*ddx_3d + nu * lap_3d
       case (2)  ! IMEX
          this%opE_3d = -v(1)*ddx_3d
          this%opI_3d =  nu*lap_3d
       case DEFAULT
          print *,'Bad case for ark_stat in f_eval ', ark_stat
          call exit(0)
       end select

       deallocate(lap_3d)
       deallocate(ddx_3d)
    else

    end if

  end subroutine initialize

  subroutine destroy(this,pf,level_index)
    use probin, only: ndim
    class(my_stepper_t), intent(inout) :: this
    type(pf_pfasst_t),  target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    !>  Call base stepper destroy
    call this%ark_destroy(pf,level_index)

    !> Nuke the FFT operators
    if (ndim == 1) then
       deallocate(this%opE_1d)
       deallocate(this%opI_1d)
    else if (ndim == 2) then
       deallocate(this%opE_2d)
       deallocate(this%opI_2d)
    else if (ndim == 3) then
       deallocate(this%opE_3d)
       deallocate(this%opI_3d)
    else

    end if

    !>  Free up FFT stuff
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)

  end subroutine destroy
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    use probin, only:  ark_stat, nu, v, ndim
    class(my_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece  !  Which piece to solve for
    
    real(pfdp),      pointer :: yvec_1d(:), fvec_1d(:)
    real(pfdp),      pointer :: yvec_2d(:,:), fvec_2d(:,:)
    real(pfdp),      pointer :: yvec_3d(:,:,:), fvec_3d(:,:,:)
    type(pf_fft_t),     pointer :: fft


    if (ndim == 1) then    
       !  Grab the arrays from the encap
       yvec_1d  => get_array1d(y)
       fvec_1d => get_array1d(f)
       fft => this%fft_tool

       ! Apply spectral operators using the FFT convolution function
       select case (piece)
       case (1)  ! Explicit piece
          call fft%conv(yvec_1d,this%opE_1d,fvec_1d)
       case (2)  ! Implicit piece
          call fft%conv(yvec_1d,this%opI_1d,fvec_1d)
       case DEFAULT
          print *,'Bad case for piece in f_eval ', piece
          call exit(0)
       end select
    else if (ndim == 2) then
       !  Grab the arrays from the encap
       yvec_2d  => get_array2d(y)
       fvec_2d => get_array2d(f)
       fft => this%fft_tool

       ! Apply spectral operators using the FFT convolution function
       select case (piece)
       case (1)  ! Explicit piece
          call fft%conv(yvec_2d,this%opE_2d,fvec_2d)
       case (2)  ! Implicit piece
          call fft%conv(yvec_2d,this%opI_2d,fvec_2d)
       case DEFAULT
          print *,'Bad case for piece in f_eval ', piece
          call exit(0)
       end select
    else if (ndim == 3) then
       !  Grab the arrays from the encap
       yvec_3d  => get_array3d(y)
       fvec_3d => get_array3d(f)
       fft => this%fft_tool

       ! Apply spectral operators using the FFT convolution function
       select case (piece)
       case (1)  ! Explicit piece
          call fft%conv(yvec_3d,this%opE_3d,fvec_3d)
       case (2)  ! Implicit piece
          call fft%conv(yvec_3d,this%opI_3d,fvec_3d)
       case DEFAULT
          print *,'Bad case for piece in f_eval ', piece
          call exit(0)
       end select
    else

    end if

  end subroutine f_eval

  ! Solve for y and return f2 also
  !   y-dtq*f(y,t) = rhs
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  ark_stat, nu, v, ndim
    class(my_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y       !  The solution we seek
    real(pfdp),          intent(in   ) :: t       !  Equation time of implicit solve
    real(pfdp),          intent(in   ) :: dtq     !  The 
    class(pf_encap_t),   intent(in   ) :: rhs     !  The right hand side of the solve
    integer,             intent(in   ) :: level_index !  Which level this is
    class(pf_encap_t),   intent(inout) :: f       !  The function value
    integer,             intent(in   ) :: piece   !  Designates which piece to solve for (here implicit)

    real(pfdp),      pointer :: yvec_1d(:), rhsvec_1d(:), fvec_1d(:)
    real(pfdp),      pointer :: yvec_2d(:,:), rhsvec_2d(:,:), fvec_2d(:,:)
    real(pfdp),      pointer :: yvec_3d(:,:,:), rhsvec_3d(:,:,:), fvec_3d(:,:,:)
    type(pf_fft_t),     pointer :: fft

    if (ndim == 1) then
       !  Grab the arrays from the encaps
       yvec_1d  => get_array1d(y)
       rhsvec_1d => get_array1d(rhs)
       fvec_1d => get_array1d(f)

       if (ark_stat .eq. 0)  then
          print *,'We should not be calling fcomp for fully explicit'
          yvec_1d=rhsvec_1d
          fvec_1d=0.0_pfdp
          return
       endif

       ! Grab the fft workspace
       fft => this%fft_tool

       if (piece == 2) then
          ! Apply the inverse operator with the FFT convolution
          call fft%conv(rhsvec_1d,1.0_pfdp/(1.0_pfdp - dtq*this%opI_1d),yvec_1d)

          !  The function is easy to derive
          fvec_1d = (yvec_1d - rhsvec_1d) / dtq
       else
          print *,'Bad piece in f_comp ',piece
          call exit(0)
       end if
    else if (ndim == 2) then
       !  Grab the arrays from the encaps
       yvec_2d  => get_array2d(y)
       rhsvec_2d => get_array2d(rhs)
       fvec_2d => get_array2d(f)

       if (ark_stat .eq. 0)  then
          print *,'We should not be calling fcomp for fully explicit'
          yvec_2d=rhsvec_2d
          fvec_2d=0.0_pfdp
          return
       endif

       ! Grab the fft workspace
       fft => this%fft_tool

       if (piece == 2) then
          ! Apply the inverse operator with the FFT convolution
          call fft%conv(rhsvec_2d,1.0_pfdp/(1.0_pfdp - dtq*this%opI_2d),yvec_2d)

          !  The function is easy to derive
          fvec_2d = (yvec_2d - rhsvec_2d) / dtq
       else
          print *,'Bad piece in f_comp ',piece
          call exit(0)
       end if
    else if (ndim == 3) then
       !  Grab the arrays from the encaps
       yvec_3d  => get_array3d(y)
       rhsvec_3d => get_array3d(rhs)
       fvec_3d => get_array3d(f)

       if (ark_stat .eq. 0)  then
          print *,'We should not be calling fcomp for fully explicit'
          yvec_3d=rhsvec_3d
          fvec_3d=0.0_pfdp
          return
       endif

       ! Grab the fft workspace
       fft => this%fft_tool

       if (piece == 2) then
          ! Apply the inverse operator with the FFT convolution
          call fft%conv(rhsvec_3d,1.0_pfdp/(1.0_pfdp - dtq*this%opI_3d),yvec_3d)

          !  The function is easy to derive
          fvec_3d = (yvec_3d - rhsvec_3d) / dtq
       else
          print *,'Bad piece in f_comp ',piece
          call exit(0)
       end if
    else
    end if

  end subroutine f_comp

end module my_stepper
