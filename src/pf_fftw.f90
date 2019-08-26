!!  Module of FFT based routines using fftw
!!
! This file is part of LIBPFASST.
!
!>  Module for providing FFTs based on fftw
!!  To use this module, fftw must be installed.
!!  This can be done by the libpfasst make system with the comment
!!   > make fftw3
module pf_mod_fftpackage
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_fft_abs
  use, intrinsic :: iso_c_binding
  
  implicit none
  include 'fftw3.f03'
  
  !>  Variables and storage for FFTW
  type,extends(pf_fft_abs_t) :: pf_fft_t
     type(c_ptr) :: ffftw, ifftw                     !! fftw pointers
     real(pfdp) :: normfact                          !! normalization factor
   contains
     procedure :: fft_setup
     procedure :: fft_destroy
     procedure :: fftf
     procedure :: fftb
     !  Interpolate in spectral space
     procedure , private :: interp_1d,interp_2d, interp_3d,zinterp_1d,zinterp_2d, zinterp_3d
     generic   :: interp =>interp_1d,interp_2d,interp_3d,zinterp_1d,zinterp_2d,zinterp_3d          
  end type pf_fft_t 
  
contains

  !> Initialize the package
  subroutine fft_setup(this, grid_shape, dim, grid_size)
    class(pf_fft_t), intent(inout) :: this
    integer,             intent(in   ) :: dim
    integer,             intent(in   ) :: grid_shape(dim)    
    real(pfdp), optional, intent(in   ) :: grid_size(dim)    

    integer     :: nx,ny,nz
    real(pfdp)  :: kx,ky
    type(c_ptr) :: wk

    this%dim=dim
    nx=grid_shape(1)
    this%nx = nx

    ! Defaults for grid_size
    this%Lx = 1.0_pfdp
    this%Ly = 1.0_pfdp
    this%Lz = 1.0_pfdp
    

    select case (dim)
    case (1)            
       if(present(grid_size)) this%Lx = grid_size(1)
       this%normfact=real(nx,pfdp)
       wk = fftw_alloc_complex(int(nx, c_size_t))
       call c_f_pointer(wk, this%wk_1d, [nx])          
       
       this%ffftw = fftw_plan_dft_1d(nx, &
            this%wk_1d, this%wk_1d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifftw = fftw_plan_dft_1d(nx, &
            this%wk_1d, this%wk_1d, FFTW_BACKWARD, FFTW_ESTIMATE)
    case (2)  
       if(present(grid_size)) then
          this%Lx = grid_size(1)
          this%Ly = grid_size(2)
       end if
       
       ny=grid_shape(2)          
       this%ny = ny
       this%normfact=real(nx*ny,pfdp)
       ! create in-place, complex fft plans
       wk = fftw_alloc_complex(int(nx*ny, c_size_t))
       call c_f_pointer(wk, this%wk_2d, [nx,ny])

       this%ffftw = fftw_plan_dft_2d(nx,ny, &
            this%wk_2d, this%wk_2d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifftw = fftw_plan_dft_2d(nx,ny, &
            this%wk_2d, this%wk_2d, FFTW_BACKWARD, FFTW_ESTIMATE)
    case (3)
       if(present(grid_size))then
          this%Lx = grid_size(1)
          this%Ly = grid_size(2)
          this%Lz = grid_size(3)
       end if
    
       ny=grid_shape(2)          
       nz=grid_shape(3)          
       this%ny = ny
       this%nz = nz
       this%normfact=real(nx*ny*nz,pfdp)          
       wk = fftw_alloc_complex(int(nx*ny*nz, c_size_t))
       call c_f_pointer(wk, this%wk_3d, [nx,ny,nz])
       this%ffftw = fftw_plan_dft_3d(nx,ny,nz, &
            this%wk_3d, this%wk_3d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifftw = fftw_plan_dft_3d(nx,ny,nz, &
            this%wk_3d, this%wk_3d, FFTW_BACKWARD, FFTW_ESTIMATE)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)
    end select
  end subroutine fft_setup

  !>  Destroy the package
  subroutine fft_destroy(this)
    class(pf_fft_t), intent(inout) :: this
    call fftw_destroy_plan(this%ffftw)
    call fftw_destroy_plan(this%ifftw)
    select case (this%dim)
    case (1)            
       deallocate(this%wk_1d)
    case (2)            
       deallocate(this%wk_2d)
    case (3)            
       deallocate(this%wk_3d)
    case DEFAULT
         call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)
    end select
  end subroutine fft_destroy
  
  !>  Routine to take foreward FFT
  subroutine fftf(this)
    class(pf_fft_t), intent(inout) :: this
    select case (this%dim)       
    case (1)
       this%wk_1d=this%wk_1d/this%normfact                 
       call fftw_execute_dft(this%ffftw, this%wk_1d, this%wk_1d)
    case (2)            
       this%wk_2d=this%wk_2d/this%normfact          
       call fftw_execute_dft(this%ffftw, this%wk_2d, this%wk_2d)
    case (3)            
       this%wk_3d=this%wk_3d/this%normfact
       call fftw_execute_dft(this%ffftw, this%wk_3d, this%wk_3d)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)       
    end select
  end subroutine fftf

  !>  Routine to take inverse or backward FFT
  subroutine fftb(this)
    class(pf_fft_t), intent(inout) :: this
    !  Normalize the fft
    select case (this%dim)       
    case (1)
       call fftw_execute_dft(this%ifftw, this%wk_1d, this%wk_1d)
    case (2)
       call fftw_execute_dft(this%ifftw, this%wk_2d, this%wk_2d)
    case (3)            
       call fftw_execute_dft(this%ifftw, this%wk_3d, this%wk_3d)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)
    end select
  end subroutine fftb

      subroutine interp_1d(this, yvec_c, fft_f,yvec_f)
!      use pf_mod_fftpackage, only: pf_fft_t
!        class(pf_fft_abs_t), intent(inout) :: this
        class(pf_fft_t), intent(inout) :: this
    real(pfdp), intent(inout),  pointer :: yvec_f(:)
    real(pfdp), intent(inout),  pointer :: yvec_c(:)
    type(pf_fft_t),  pointer,intent(in) :: fft_f
    integer :: nx_f, nx_c

    complex(pfdp),         pointer :: wk_f(:), wk_c(:)

    call this%get_wk_ptr(wk_c)
    call fft_f%get_wk_ptr(wk_f)
    
    wk_c=yvec_c
    !  internal forward fft call    
    call this%fftf()    

    call this%zinterp_1d(wk_c, wk_f)
    wk_f=wk_f

    !  internal inverse fft call
    call fft_f%fftb()

    yvec_f=real(wk_f)
    
  end subroutine interp_1d
  subroutine interp_2d(this, yvec_c, fft_f,yvec_f)
        class(pf_fft_t), intent(inout) :: this
!    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout),  pointer :: yvec_f(:,:)
    real(pfdp), intent(inout),  pointer :: yvec_c(:,:)
    type(pf_fft_t),  pointer,intent(in) :: fft_f

    complex(pfdp),         pointer :: wk_f(:,:), wk_c(:,:)

    call this%get_wk_ptr(wk_c)
    call fft_f%get_wk_ptr(wk_f)
    
    wk_c=yvec_c
    !  internal forward fft call    
    call this%fftf()    

    call this%zinterp_2d(wk_c, wk_f)
    wk_f=wk_f

    !  internal inverse fft call
    call fft_f%fftb()

    yvec_f=real(wk_f)
    
  end subroutine interp_2d
  subroutine interp_3d(this, yvec_c, fft_f,yvec_f)
         class(pf_fft_t), intent(inout) :: this
!   class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout),  pointer :: yvec_f(:,:,:)
    real(pfdp), intent(inout),  pointer :: yvec_c(:,:,:)
    type(pf_fft_t),  pointer,intent(in) :: fft_f
    complex(pfdp),         pointer :: wk_f(:,:,:), wk_c(:,:,:)

    call this%get_wk_ptr(wk_c)
    call fft_f%get_wk_ptr(wk_f)
    
    wk_c=yvec_c
    !  internal forward fft call    
    call this%fftf()    

    call this%zinterp_3d(wk_c, wk_f)
    wk_f=wk_f

    !  internal inverse fft call
    call fft_f%fftb()

    yvec_f=real(wk_f)
    

  end subroutine interp_3d

  !>  Interpolate from coarse  level to fine
  subroutine zinterp_1d(this, yhat_c, yhat_f)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp),         pointer :: yhat_f(:), yhat_c(:)
    integer :: nx_f, nx_c

    nx_f = size(yhat_f)
    nx_c = size(yhat_c)

    yhat_f = 0.0_pfdp
    yhat_f(1:nx_c/2) = yhat_c(1:nx_c/2)
    yhat_f(nx_f-nx_c/2+2:nx_f) = yhat_c(nx_c/2+2:nx_c)
    
  end subroutine zinterp_1d
  !>  Interpolate from coarse  level to fine
  subroutine zinterp_2d(this, yhat_c, yhat_f)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp),         pointer :: yhat_f(:,:), yhat_c(:,:)

    integer :: nx_f(2), nx_c(2),nf1,nf2,nc1,nc2

    nx_f = shape(yhat_f)
    nx_c = shape(yhat_c)

    if (nx_f(1) .eq. nx_c(1) .and. nx_f(2) .eq. nx_c(2)) then
       yhat_f = yhat_c
       return
    end if
    
  
    nf1=nx_f(1)-nx_c(1)/2+2
    nf2=nx_f(2)-nx_c(2)/2+2
    nc1=nx_c(1)/2+2
    nc2=nx_c(2)/2+2

    yhat_f = 0.0_pfdp
    yhat_f(1:nx_c(1)/2,1:nx_c(2)/2) = yhat_c(1:nx_c(1)/2,1:nx_c(2)/2)
    yhat_f(nf1:nx_f(1),1:nx_c(2)/2) = yhat_c(nc1:nx_c(1),1:nx_c(2)/2)    
    yhat_f(1:nx_c(1)/2,nf2:nx_f(2)) = yhat_c(1:nx_c(1)/2,nc2:nx_c(2)) 
    yhat_f(nf1:nx_f(1),nf2:nx_f(2)) = yhat_c(nc1:nx_c(1),nc2:nx_c(2)) 
   
  end subroutine zinterp_2d

  !>  Interpolate from coarse  level to fine
  subroutine zinterp_3d(this, yhat_c, yhat_f)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp),         pointer :: yhat_f(:,:,:), yhat_c(:,:,:)

  integer :: nx_f(3), nx_c(3),nf1,nf2,nf3,nc1,nc2,nc3


    nx_f = size(yhat_f)
    nx_c = size(yhat_c)

    yhat_f = 0.0_pfdp
  
    nx_f = shape(yhat_f)
    nx_c = shape(yhat_c)
    
    nf1=nx_f(1)-nx_c(1)/2+2
    nf2=nx_f(2)-nx_c(2)/2+2
    nf3=nx_f(3)-nx_c(3)/2+2
    nc1=nx_c(1)/2+2
    nc2=nx_c(2)/2+2
    nc3=nx_c(3)/2+2

    yhat_f = 0.0_pfdp
    yhat_f(1:nx_c(1)/2,1:nx_c(2)/2,1:nx_c(3)/2) = yhat_c(1:nx_c(1)/2,1:nx_c(2)/2,1:nx_c(3)/2)
    yhat_f(nf1:nx_f(1),1:nx_c(2)/2,1:nx_c(3)/2) = yhat_c(nc1:nx_c(1),1:nx_c(2)/2,1:nx_c(3)/2)    
    yhat_f(1:nx_c(1)/2,nf2:nx_f(2),1:nx_c(3)/2) = yhat_c(1:nx_c(1)/2,nc2:nx_c(2),1:nx_c(3)/2) 
    yhat_f(nf1:nx_f(1),nf2:nx_f(2),1:nx_c(3)/2) = yhat_c(nc1:nx_c(1),nc2:nx_c(2),1:nx_c(3)/2) 

    yhat_f(1:nx_c(1)/2,1:nx_c(2)/2,nf3:nx_f(3)) = yhat_c(1:nx_c(1)/2,1:nx_c(2)/2,nc3:nx_c(3))
    yhat_f(nf1:nx_f(1),1:nx_c(2)/2,nf3:nx_f(3)) = yhat_c(nc1:nx_c(1),1:nx_c(2)/2,nc3:nx_c(3))    
    yhat_f(1:nx_c(1)/2,nf2:nx_f(2),nf3:nx_f(3)) = yhat_c(1:nx_c(1)/2,nc2:nx_c(2),nc3:nx_c(3)) 
    yhat_f(nf1:nx_f(1),nf2:nx_f(2),nf3:nx_f(3)) = yhat_c(nc1:nx_c(1),nc2:nx_c(2),nc3:nx_c(3)) 
    
  end subroutine zinterp_3d

end module pf_mod_fftpackage



