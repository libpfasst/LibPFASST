!!  Module of FFT based routines using fftpack
!
! This file is part of LIBPFASST.
!
!>  Module for using fftpack
module pf_mod_fft_abs
  use pf_mod_dtype
  use pf_mod_utils
  implicit none
  
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp
  !>  Variables and storage for FFT
  type, abstract  :: pf_fft_abs_t
     integer ::   nx,ny,nz                         !! grid sizes
     integer ::   dim                              !! spatial dimension
     real(pfdp) :: Lx, Ly, Lz                      !! domain size
     complex(pfdp), pointer :: wk_1d(:)            ! work space
     complex(pfdp), pointer :: wk_2d(:,:)          ! work space
     complex(pfdp), pointer :: wk_3d(:,:,:)        ! work space                    
   contains
     procedure(pf_fft_s_p),deferred :: fft_setup
     procedure(pf_fft_p),deferred :: fft_destroy
     procedure(pf_fft_p),deferred :: fftf   !  Forward FFT
     procedure(pf_fft_p),deferred :: fftb   !  Inverse (backward)  FFT
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
     !  Construct spectral derivative
     procedure , private :: make_deriv_1d,make_deriv_2d, make_deriv_3d
     generic   :: make_deriv =>make_deriv_1d,make_deriv_2d,make_deriv_3d
  end type pf_fft_abs_t

  interface
     subroutine pf_fft_s_p(this, grid_shape, dim, grid_size)
       import pf_fft_abs_t,pfdp
       class(pf_fft_abs_t), intent(inout) :: this
       integer,              intent(in   ) :: dim
       integer,              intent(in   ) :: grid_shape(dim)
       real(pfdp), optional, intent(in   ) :: grid_size(dim)
     end subroutine pf_fft_s_p
     subroutine pf_fft_p(this)
       import pf_fft_abs_t,pfdp
       class(pf_fft_abs_t), intent(inout) :: this
     end subroutine pf_fft_p
  end interface
  
contains
  ! START private convolution procedures
  
  ! Convolve g with spectral op and return in c
  subroutine conv_1d(this, g,op,c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:)
    complex(pfdp), intent(in) :: op(:)
    real(pfdp), intent(inout) :: c(:)

    this%wk_1d=g
    call this%fftf()
    this%wk_1d = this%wk_1d * op
    call this%fftb()
    c=real(this%wk_1d,pfdp)
  end subroutine conv_1d

  ! Convolve g with spectral op and return in c
  subroutine conv_2d(this, g,op,c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:)
    complex(pfdp), intent(in) :: op(:,:)
    real(pfdp), intent(inout) :: c(:,:)
    this%wk_2d=g
    ! Compute Convolution
    call this%fftf()
    this%wk_2d = this%wk_2d * op
    call this%fftb()
    c=real(this%wk_2d,pfdp)        
  end subroutine conv_2d
  
  subroutine conv_3d(this, g,op,c)
        class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:,:)
    complex(pfdp), intent(in) :: op(:,:,:)
    real(pfdp), intent(inout) :: c(:,:,:)
    this%wk_3d=g
    ! Compute Convolution
    call this%fftf()
    this%wk_3d = this%wk_3d * op
    call this%fftb()
    c=real(this%wk_3d,pfdp)            
  end subroutine conv_3d

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
    
  subroutine zconv_1d(this, g)
    ! Variable Types
        class(pf_fft_abs_t), intent(inout) :: this
        real(pfdp), intent(in) :: g(:)
        ! Compute Convolution
        call this%fftb()
        this%wk_1d = this%wk_1d * g
        call this%fftf()
    end subroutine zconv_1d

    subroutine zconv_2d(this, g)
        ! Variable Types
        class(pf_fft_abs_t), intent(inout) :: this
        real(pfdp), intent(in) :: g(:,:)
        ! Compute Convolution
        call this%fftb()
        this%wk_2d = this%wk_2d * g
        call this%fftf()
    end subroutine zconv_2d

    subroutine zconv_3d(this, g)
        ! Variable Types
        class(pf_fft_abs_t), intent(inout) :: this
        real(pfdp), intent(in) :: g(:,:,:)
        ! Compute Convolution
        call this%fftb()
        this%wk_3d = this%wk_3d * g
        call this%fftf()
    end subroutine zconv_3d
    
    subroutine make_lap_1d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:)
      
      integer     :: i,nx
      real(pfdp)  :: kx, Lx
      
      nx=this%nx
      Lx=this%Lx
      do i = 1, nx
         if (i <= nx/2+1) then
            kx = two_pi / Lx * dble(i-1)
         else
            kx = two_pi / Lx * dble(-nx + i - 1)
         end if
         lap(i) = -kx**2
      end do
      
    end subroutine make_lap_1d
    subroutine make_deriv_1d(this, ddx)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: ddx(:)
      
      integer     :: i,nx
      real(pfdp)  :: kx, Lx
      
      nx=this%nx
      Lx=this%Lx
      do i = 1, nx
         if (i <= nx/2+1) then
            kx = two_pi / Lx * dble(i-1)
         else
            kx = two_pi / Lx * dble(-nx + i - 1)
         end if
         
         ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx
      end do
    end subroutine make_deriv_1d
     
    subroutine make_lap_2d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:,:)
      
      integer     :: i,j,nx,ny
      real(pfdp)  :: kx,ky,Lx,Ly
      
      nx=this%nx
      ny=this%ny
      Lx=this%Lx
      Ly=this%Ly
      
      do j = 1, ny
         if (j <= ny/2+1) then
            ky = two_pi / Ly * dble(j-1)
         else
            ky = two_pi / Ly * dble(-ny + j - 1)
         end if
         do i = 1, nx
            if (i <= nx/2+1) then
               kx = two_pi / Lx * dble(i-1)
            else
               kx = two_pi / Lx * dble(-nx + i - 1)
            end if
            
            lap(i,j) = -(kx**2+ky**2)
         end do
      end do
    end subroutine make_lap_2d

    subroutine make_deriv_2d(this, deriv,dir)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: deriv(:,:)
      integer, intent(in) :: dir
      
      integer     :: i,j,nx,ny
      real(pfdp)  :: kx,ky,Lx,Ly
      
      nx=this%nx
      ny=this%ny
      Lx=this%Lx
      Ly=this%Ly
      
      do j = 1, ny
         if (j <= ny/2+1) then
            ky = two_pi / Ly * dble(j-1)
         else
            ky = two_pi / Ly * dble(-ny + j - 1)
         end if
         do i = 1, nx
            if (i <= nx/2+1) then
               kx = two_pi / Lx * dble(i-1)
            else
               kx = two_pi / Lx * dble(-nx + i - 1)
            end if
            
            if (dir .eq. 1) then
               deriv(i,j) = (0.0_pfdp,1.0_pfdp)*kx
            else
               deriv(i,j) = (0.0_pfdp,1.0_pfdp)*ky
            endif
         end do
      end do
    end subroutine make_deriv_2d


    subroutine make_lap_3d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:,:,:)
      
      integer     :: i,j,k,nx,ny,nz
      real(pfdp)  :: kx,ky,kz,Lx,Ly,Lz
      
      nx=this%nx
      ny=this%ny
      nz=this%nz
      Lx=this%Lx
      Ly=this%Ly
      Lz=this%Lz
      do k = 1,nz
         if (k <= nz/2+1) then
            kz = two_pi / Lz * dble(k-1)
         else
            kz = two_pi / Ly * dble(-nz + k - 1)
         end if
         do j = 1, ny
            if (j <= ny/2+1) then
               ky = two_pi / Ly * dble(j-1)
            else
               ky = two_pi / Ly * dble(-ny + j - 1)
            end if
            do i = 1, nx
               if (i <= nx/2+1) then
                  kx = two_pi / Lx * dble(i-1)
               else
                  kx = two_pi / Lx * dble(-nx + i - 1)
               end if
               lap(i,j,k) = -(kx**2+ky**2+kz**2)
            end do
         end do
      end do
      
    end subroutine make_lap_3d

    subroutine make_deriv_3d(this, deriv,dir)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: deriv(:,:,:)
      integer, intent(in) :: dir
      
      integer     :: i,j,k,nx,ny,nz
      real(pfdp)  :: kx,ky,kz,Lx,Ly,Lz
      
      nx=this%nx
      ny=this%ny
      nz=this%nz       
      Lx=this%Lx
      Ly=this%Ly
      Lz=this%Lz
      
      select case (dir)
      case (1)  
         do i = 1, nx
            if (i <= nx/2+1) then
               kx = two_pi / Lx * dble(i-1)
            else
               kx = two_pi / Lx * dble(-nx + i - 1)
            end if
            deriv(i,:,:) = (0.0_pfdp,1.0_pfdp)*kx
         end do
      case (2)
         do j = 1, ny
            if (j <= ny/2+1) then
               ky = two_pi / Ly * dble(j-1)
            else
               ky = two_pi / Ly * dble(-ny + j - 1)
            end if
            deriv(:,j,:) = (0.0_pfdp,1.0_pfdp)*ky
         end do
      case (3)
         do k = 1, nz
            if (k <= nz/2+1) then
               kz = two_pi / Lz * dble(k-1)
            else
               kz = two_pi / Ly * dble(-nz + k - 1)
            end if
            deriv(:,:,k) = (0.0_pfdp,1.0_pfdp)*kz
         end do
      case DEFAULT
         call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',dir)
      end select

    end subroutine make_deriv_3d
    
  end module pf_mod_fft_abs
  
   



