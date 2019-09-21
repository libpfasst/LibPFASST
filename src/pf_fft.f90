!!  Module of FFT based routines using fftpack
!
! This file is part of LIBPFASST.
!
!>  Module for using fftpack
module pf_mod_fft_abs
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  !>  Variables and storage for FFT
  type, abstract  :: pf_fft_abs_t
     integer ::   nx,ny,nz                         !! grid sizes
     integer ::   ndim                             !! number of spatial dimensions
     real(pfdp) :: Lx, Ly, Lz                      !! domain size
     complex(pfdp), pointer :: wk_1d(:)            ! work space
     complex(pfdp), pointer :: wk_2d(:,:)          ! work space
     complex(pfdp), pointer :: wk_3d(:,:,:)        ! work space                    
     real(pfdp), allocatable :: kx(:),ky(:),kz(:)              ! work space                    
   contains
     procedure(pf_fft_s_p),deferred :: fft_setup
     procedure(pf_fft_p),deferred :: fft_destroy
     procedure(pf_fft_p),deferred :: fftf   !  Forward FFT
     procedure(pf_fft_p),deferred :: fftb   !  Inverse (backward)  FFT
     !  FFT
     procedure, private  :: fft_1d, fft_2d, fft_3d, zfft_1d, zfft_2d, zfft_3d  
     generic :: fft => fft_1d, fft_2d, fft_3d, zfft_1d, zfft_2d, zfft_3d  
     !  Inverse FFT
     procedure, private  :: ifft_1d, ifft_2d, ifft_3d,izfft_1d, izfft_2d, izfft_3d  
     generic :: ifft => ifft_1d, ifft_2d, ifft_3d,izfft_1d, izfft_2d, izfft_3d    
     !  Convolution in spectral space     
     procedure, private  :: conv_1d, conv_2d, conv_3d  
     generic :: conv => conv_1d, conv_2d, conv_3d  
     !  Complex convolution in real space     
     procedure, private :: zconv_1d, zconv_2d, zconv_3d  
     generic :: zconv => zconv_1d, zconv_2d, zconv_3d  
     !  Convenience function to grab pointer to workspace
     procedure , private :: get_wk_ptr_1d, get_wk_ptr_2d, get_wk_ptr_3d  
     generic   :: get_wk_ptr =>get_wk_ptr_1d,get_wk_ptr_2d,get_wk_ptr_3d
     !  Construct spectral Laplacian
     procedure , private :: make_lap_1d,make_lap_2d, make_lap_3d
     generic   :: make_lap =>make_lap_1d,make_lap_2d,make_lap_3d
     !  Construct inverse spectral Laplacian
     procedure , private :: make_ilap_1d,make_ilap_2d, make_ilap_3d
     generic   :: make_ilap =>make_ilap_1d,make_ilap_2d,make_ilap_3d
     !  Construct spectral derivative
     procedure , private :: make_deriv_1d,make_deriv_2d, make_deriv_3d
     generic   :: make_deriv =>make_deriv_1d,make_deriv_2d,make_deriv_3d
     !  Restrict in spectral space
     procedure , private :: restrict_1d,restrict_2d, restrict_3d,zrestrict_1d,zrestrict_2d,zrestrict_3d
     generic   :: restrict =>restrict_1d,restrict_2d,restrict_3d,zrestrict_1d,zrestrict_2d,zrestrict_3d
     
  end type pf_fft_abs_t

  interface
     subroutine pf_fft_s_p(this, grid_shape, ndim, grid_size)
       import pf_fft_abs_t,pfdp
       class(pf_fft_abs_t), intent(inout) :: this
       integer,              intent(in   ) :: ndim
       integer,              intent(in   ) :: grid_shape(ndim)
       real(pfdp), optional, intent(in   ) :: grid_size(ndim)
     end subroutine pf_fft_s_p
     subroutine pf_fft_p(this)
       import pf_fft_abs_t,pfdp
       class(pf_fft_abs_t), intent(inout) :: this
     end subroutine pf_fft_p
  end interface
  
contains

 subroutine get_wk_ptr_1d(this,wk) 
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), pointer,intent(inout) :: wk(:)              ! work space
    wk=>this%wk_1d
  end subroutine get_wk_ptr_1d
  subroutine get_wk_ptr_2d(this,wk) 
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), pointer,intent(inout) :: wk(:,:)              ! work space
    wk=>this%wk_2d
  end subroutine get_wk_ptr_2d
  subroutine get_wk_ptr_3d(this,wk) 
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), pointer,intent(inout) :: wk(:,:,:)              ! work space
    wk=>this%wk_3d
  end subroutine get_wk_ptr_3d
  
  !++++++++++ Forward FFT real to complex  ++++++++++++++++
  subroutine fft_1d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:)
    complex(pfdp), intent(inout) :: ghat(:)

    this%wk_1d=g
    call this%fftf()
    ghat=this%wk_1d
  end subroutine fft_1d

  subroutine fft_2d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:)
    complex(pfdp), intent(inout) :: ghat(:,:)

    this%wk_2d=g
    call this%fftf()
    ghat=this%wk_2d
  end subroutine fft_2d

    subroutine fft_3d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:,:)
    complex(pfdp), intent(inout) :: ghat(:,:,:)

    this%wk_3d=g
    call this%fftf()
    ghat=this%wk_3d
  end subroutine fft_3d

  !++++++++++ Backward FFT complex to real   ++++++++++++++++
  subroutine ifft_1d(this,ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:)
    complex(pfdp), intent(in) :: ghat(:)

    this%wk_1d=ghat
    call this%fftb()
    g=REAL(this%wk_1d,pfdp)
  end subroutine ifft_1d

  subroutine ifft_2d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:,:)
    complex(pfdp), intent(in) :: ghat(:,:)

    this%wk_2d=ghat
    call this%fftb()
    g=REAL(this%wk_2d,pfdp)
  end subroutine ifft_2d

    subroutine ifft_3d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:,:,:)
    complex(pfdp), intent(in) :: ghat(:,:,:)

    this%wk_3d=ghat
    call this%fftb()
    g=REAL(this%wk_3d)
  end subroutine ifft_3d

  !++++++++++ Forward FFT complex to complex   ++++++++++++++++
  subroutine zfft_1d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: g(:)
    complex(pfdp), intent(inout) :: ghat(:)

    this%wk_1d=g
    call this%fftf()
    ghat=this%wk_1d
  end subroutine zfft_1d

    subroutine zfft_2d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: g(:,:)
    complex(pfdp), intent(inout) :: ghat(:,:)

    this%wk_2d=g
    call this%fftf()
    ghat=this%wk_2d
  end subroutine zfft_2d

  subroutine zfft_3d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: g(:,:,:)
    complex(pfdp), intent(inout) :: ghat(:,:,:)

    this%wk_3d=g
    call this%fftf()
    ghat=this%wk_3d
  end subroutine zfft_3d

  !++++++++++ Backward FFT complex to complex   ++++++++++++++++
  subroutine izfft_1d(this,ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: ghat(:)
    complex(pfdp), intent(inout) :: g(:)
    
    this%wk_1d=ghat
    call this%fftb()
    g=this%wk_1d
  end subroutine izfft_1d

  ! Take forward FFT
  subroutine izfft_2d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: ghat(:,:)
    complex(pfdp), intent(inout) :: g(:,:)

    this%wk_2d=ghat
    call this%fftb()
    g=this%wk_2d
  end subroutine izfft_2d

  subroutine izfft_3d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: ghat(:,:,:)
    complex(pfdp), intent(inout) :: g(:,:,:)

    this%wk_3d=ghat
    call this%fftb()
    g=this%wk_3d
  end subroutine izfft_3d

  ! Convolve g with spectral op and return in c
  subroutine conv_1d(this, g,op,c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:)
    complex(pfdp), intent(in) :: op(:)
    real(pfdp), intent(inout) :: c(:)

    this%wk_1d=g
    call this%fftf()
    this%wk_1d = this%wk_1d * op
    call this%fftb()
    c=REAL(this%wk_1d,pfdp)
  end subroutine conv_1d

  ! Convolve g with spectral op and return in c
  subroutine conv_2d(this, g,op,c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:)
    complex(pfdp), intent(in) :: op(:,:)
    real(pfdp), intent(inout) :: c(:,:)

    this%wk_2d=g
    call this%fftf()
    this%wk_2d = this%wk_2d * op
    call this%fftb()
    c=REAL(this%wk_2d,pfdp)        
  end subroutine conv_2d
  
  subroutine conv_3d(this, g,op,c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:,:)
    complex(pfdp), intent(in) :: op(:,:,:)
    real(pfdp), intent(inout) :: c(:,:,:)

    this%wk_3d=g
    call this%fftf()
    this%wk_3d = this%wk_3d * op
    call this%fftb()
    c=REAL(this%wk_3d,pfdp)            
  end subroutine conv_3d

  ! Convolution in real space 
  subroutine zconv_1d(this, ghat,op,chat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: ghat(:),op(:)
    complex(pfdp), intent(inout) :: chat(:)

    this%wk_1d=ghat
    call this%fftb()
    this%wk_1d = this%wk_1d * op
    call this%fftf()
    chat=this%wk_1d

  end subroutine zconv_1d

  subroutine zconv_2d(this, ghat,op,chat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: ghat(:,:),op(:,:)
    complex(pfdp), intent(inout) ::chat(:,:)
    
      
    this%wk_2d=ghat
    call this%fftb()
    this%wk_2d = this%wk_2d * op
    call this%fftf()
    chat=this%wk_2d
    
    
  end subroutine zconv_2d

  subroutine zconv_3d(this, ghat,op,chat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: ghat(:,:,:),op(:,:,:)
    complex(pfdp), intent(inout) ::chat(:,:,:)
    
    this%wk_3d=ghat
    call this%fftb()
    this%wk_3d = this%wk_3d * op
    call this%fftf()
    chat=this%wk_3d
    
  end subroutine zconv_3d

  !  Make the inverse of the Laplacian in spectral space
  subroutine make_ilap_1d(this, ilap)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: ilap(:)
    
    integer     :: i,nx
    
    nx=this%nx
    ilap(1) = 0.0_pfdp !    This sets DC component    
    do i = 2, nx
       ilap(i) = -1.0_pfdp/(this%kx(i)**2)
    end do
    
  end subroutine make_ilap_1d
  
  subroutine make_ilap_2d(this, ilap)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: ilap(:,:)
    
    integer     :: i,j,nx,ny
    
    nx=this%nx
    ny=this%ny
    
    ilap(1,1) = 0.0_pfdp !    This sets DC component
    do i = 2, nx
       ilap(i,1) = -1.0_pfdp/(this%kx(i)**2)       
    end do
    do j = 2, ny
       do i = 1, nx
          ilap(i,j) = -1.0_pfdp/(this%kx(i)**2+this%ky(j)**2)
       end do
    end do
  end subroutine make_ilap_2d
  
  subroutine make_ilap_3d(this, ilap)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: ilap(:,:,:)
    
    integer     :: i,j,k,nx,ny,nz
      
    nx=this%nx
    ny=this%ny
    nz=this%nz
    do k = 1,nz
       do j = 1, ny
          do i = 1, nx
             if (i .eq. 0  .and. j .eq. 0 ) then
                ilap(i,j,k) = 0.0_pfdp !    This sets DC component
             else
                ilap(i,j,k) = -1.0_pfdp/(this%kx(i)**2+this%ky(j)**2+this%kz(k)**2)
             end if
          end do
       end do
    end do
    
  end subroutine make_ilap_3d
  
    subroutine make_lap_1d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:)
      
      integer     :: i,nx

      nx=this%nx
      do i = 1, nx
         lap(i) = -(this%kx(i)**2)
      end do
    end subroutine make_lap_1d
    
    subroutine make_lap_2d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:,:)
      
      integer     :: i,j,nx,ny
      
      nx=this%nx
      ny=this%ny
      do j = 1, ny
         do i = 1, nx
            lap(i,j) = -(this%kx(i)**2+this%ky(j)**2)
         end do
      end do
    end subroutine make_lap_2d
    subroutine make_lap_3d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:,:,:)
      
      integer     :: i,j,k,nx,ny,nz
      
      nx=this%nx
      ny=this%ny
      nz=this%nz
      do k = 1,nz
         do j = 1, ny
            do i = 1, nx
               lap(i,j,k) = -(this%kx(i)**2+this%ky(j)**2+this%kz(k)**2)
            end do
         end do
      end do
      
    end subroutine make_lap_3d

    subroutine make_deriv_1d(this, ddx)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: ddx(:)
      
      integer     :: i,nx
      
      nx=this%nx
      do i = 1, nx
         ddx(i) = (0.0_pfdp, 1.0_pfdp) * this%kx(i)
      end do
    end subroutine make_deriv_1d
    subroutine make_deriv_2d(this, deriv,dir)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: deriv(:,:)
      integer, intent(in) :: dir
      
      integer     :: i,j,nx,ny
      
      nx=this%nx
      ny=this%ny
      
      if (dir .eq. 1) then
         do i = 1, nx
            deriv(i,:) = (0.0_pfdp,1.0_pfdp)*this%kx(i)
         end do
      else
         do j = 1, ny
            deriv(:,j) = (0.0_pfdp,1.0_pfdp)*this%ky(j)
         end do
      endif
    end subroutine make_deriv_2d

    subroutine make_deriv_3d(this, deriv,dir)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: deriv(:,:,:)
      integer, intent(in) :: dir
      
      integer     :: i,j,k,nx,ny,nz
      
      nx=this%nx
      ny=this%ny
      nz=this%nz       
      
      select case (dir)
      case (1)  
         do i = 1, nx
            deriv(i,:,:) = (0.0_pfdp,1.0_pfdp)*this%kx(i)
         end do
      case (2)
         do j = 1, ny
            deriv(:,j,:) = (0.0_pfdp,1.0_pfdp)*this%ky(j)
         end do
      case (3)
         do k = 1, nz
            deriv(:,:,k) = (0.0_pfdp,1.0_pfdp)*this%kz(k)
         end do
      case DEFAULT
         call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',dir)
      end select

    end subroutine make_deriv_3d

  !  Restrict routines that take a fine vector and produce a coarse version
  subroutine restrict_1d(this, yvec_f, yvec_c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp),         pointer :: yvec_f(:), yvec_c(:)
    integer :: nx_f, nx_c,irat
    nx_f = SIZE(yvec_f)
    nx_c = SIZE(yvec_c)
    irat  = nx_f/nx_c

    yvec_c = yvec_f(::irat)
  end subroutine restrict_1d
  subroutine restrict_2d(this, yvec_f, yvec_c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp),         pointer :: yvec_f(:,:), yvec_c(:,:)
    integer :: nx_f(2), nx_c(2), irat, jrat

    nx_f = SHAPE(yvec_f)
    nx_c = SHAPE(yvec_c)

    irat  = nx_f(1)/nx_c(1)
    jrat  = nx_f(2)/nx_c(2)
    
    yvec_c = yvec_f(::irat,::jrat)           
  end subroutine restrict_2d
  subroutine restrict_3d(this, yvec_f, yvec_c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp),         pointer :: yvec_f(:,:,:), yvec_c(:,:,:)
    integer :: nx_f(3), nx_c(3)
    integer :: irat, jrat,krat

    nx_f = SHAPE(yvec_f)
    nx_c = SHAPE(yvec_c)

    irat  = nx_f(1)/nx_c(1)
    jrat  = nx_f(2)/nx_c(2)
    krat  = nx_f(3)/nx_c(3)
    
    yvec_c = yvec_f(::irat,::jrat,::krat)           
  end subroutine restrict_3d
  
  subroutine zrestrict_1d(this, yhat_f, yhat_c)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp),  pointer :: yhat_f(:), yhat_c(:)
    integer :: nx_f, nx_c

    nx_f = SIZE(yhat_f)
    nx_c = SIZE(yhat_c)

    yhat_c=0.0_pfdp
    yhat_c(1:nx_c/2) = yhat_f(1:nx_c/2)
    yhat_c(nx_c/2+2:nx_c) = yhat_f(nx_f-nx_c/2+2:nx_f)
    
  end subroutine zrestrict_1d
  subroutine zrestrict_2d(this, yhat_f, yhat_c)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp),         pointer :: yhat_f(:,:), yhat_c(:,:)

    integer :: nx_f(2), nx_c(2),nf1,nf2,nc1,nc2
    nx_f = SHAPE(yhat_f)
    nx_c = SHAPE(yhat_c)
    
    nf1=nx_f(1)-nx_c(1)/2+2
    nf2=nx_f(2)-nx_c(2)/2+2
    nc1=nx_c(1)/2+2
    nc2=nx_c(2)/2+2

    yhat_c=0.0_pfdp
    yhat_c(1:nx_c(1)/2,1:nx_c(2)/2) = yhat_f(1:nx_c(1)/2,1:nx_c(2)/2)
    yhat_c(nc1:nx_c(1),1:nx_c(2)/2) = yhat_f(nf1:nx_f(1),1:nx_c(2)/2)
    yhat_c(1:nx_c(1)/2,nc2:nx_c(2)) = yhat_f(1:nx_c(1)/2,nf2:nx_f(2))
    yhat_c(nc1:nx_c(1),nc2:nx_c(2)) = yhat_f(nf1:nx_f(1),nf2:nx_f(2)) 
    
  end subroutine zrestrict_2d
  subroutine zrestrict_3d(this, yhat_f, yhat_c)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp),         pointer :: yhat_f(:,:,:), yhat_c(:,:,:)

    integer :: nx_f(3), nx_c(3),nf1,nf2,nf3,nc1,nc2,nc3
    
    
    
    yhat_c = 0.0_pfdp
    
    nx_f = SHAPE(yhat_f)
    nx_c = SHAPE(yhat_c)
    
    nf1=nx_f(1)-nx_c(1)/2+2
    nf2=nx_f(2)-nx_c(2)/2+2
    nf3=nx_f(3)-nx_c(3)/2+2
    nc1=nx_c(1)/2+2
    nc2=nx_c(2)/2+2
    nc3=nx_c(3)/2+2
    
    yhat_c(1:nx_c(1)/2,1:nx_c(2)/2,1:nx_c(3)/2) = yhat_f(1:nx_c(1)/2,1:nx_c(2)/2,1:nx_c(3)/2)
    yhat_c(nc1:nx_c(1),1:nx_c(2)/2,1:nx_c(3)/2) = yhat_f(nf1:nx_f(1),1:nx_c(2)/2,1:nx_c(3)/2)
    yhat_c(1:nx_c(1)/2,nc2:nx_c(2),1:nx_c(3)/2) = yhat_f(1:nx_c(1)/2,nf2:nx_f(2),1:nx_c(3)/2)
    yhat_c(nc1:nx_c(1),nc2:nx_c(2),1:nx_c(3)/2) = yhat_f(nf1:nx_f(1),nf2:nx_f(2),1:nx_c(3)/2) 
    
    yhat_c(1:nx_c(1)/2,1:nx_c(2)/2,nc3:nx_c(3)) = yhat_f(1:nx_c(1)/2,1:nx_c(2)/2,nf3:nx_f(3)) 
    yhat_c(nc1:nx_c(1),1:nx_c(2)/2,nc3:nx_c(3)) = yhat_f(nf1:nx_f(1),1:nx_c(2)/2,nf3:nx_f(3))
    yhat_c(1:nx_c(1)/2,nc2:nx_c(2),nc3:nx_c(3)) = yhat_f(1:nx_c(1)/2,nf2:nx_f(2),nf3:nx_f(3)) 
    yhat_c(nc1:nx_c(1),nc2:nx_c(2),nc3:nx_c(3)) = yhat_f(nf1:nx_f(1),nf2:nx_f(2),nf3:nx_f(3))
    
  end subroutine zrestrict_3d

  end module pf_mod_fft_abs
  
   



