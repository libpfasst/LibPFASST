!!  Module of FFT based routines using fftpack
!
! This file is part of LIBPFASST.
!
!>  Module for using fftpack
module pf_mod_fftpackage
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_fft_abs
  implicit none
  
  !>  Variables and storage for FFT
  type, extends(pf_fft_abs_t)  :: pf_fft_t
     integer ::   lensavx, lensavy, lensavz        !! workspace lengths
     real(pfdp) :: normfact                        !! normalization factor
     real(pfdp), allocatable :: wsavex(:)          ! work space
     real(pfdp), allocatable :: wsavey(:)          ! work space
     real(pfdp), allocatable :: wsavez(:)          ! work space     
     complex(pfdp), pointer :: workhatx(:)         ! work space
     complex(pfdp), pointer :: workhaty(:)         ! work space
     complex(pfdp), pointer :: workhatz(:)         ! work space
   contains
     procedure :: fft_setup 
     procedure :: fft_destroy
     procedure :: fftf  !  Forward FFT
     procedure :: fftb   !  Inverse (backward)  FFT
     !  Interpolate in spectral space
     procedure , private :: interp_1d,interp_2d, interp_3d,zinterp_1d,zinterp_2d, zinterp_3d
     generic   :: interp =>interp_1d,interp_2d,interp_3d,zinterp_1d,zinterp_2d,zinterp_3d     
  end type pf_fft_t
  
contains
  !>  Allocate and initialize FFT structure
  subroutine fft_setup(this, grid_shape, dim, grid_size)
    class(pf_fft_t), intent(inout) :: this
    integer,              intent(in   ) :: dim
    integer,              intent(in   ) :: grid_shape(dim)
    real(pfdp), optional, intent(in   ) :: grid_size(dim)    
      
    integer     :: nx,ny,nz
    this%dim=dim
    
    !  FFT Storage parameters
    nx=grid_shape(1)
    this%nx = nx
    this%lensavx = 4*nx + 15
    this%normfact=nx
    
    allocate(this%workhatx(nx))   !  complex transform
    allocate(this%wsavex(this%lensavx))
    this%Lx = 1.0_pfdp
    if(present(grid_size)) this%Lx = grid_size(1)
    !  Initialize FFT
    call ZFFTI( nx, this%wsavex )
    if (dim > 1) then
       !  FFT Storage
       ny=grid_shape(2)       
       this%ny = ny
       this%lensavy = 4*ny + 15
       this%normfact=nx*ny
       allocate(this%workhaty(ny))   !  complex transform
       allocate(this%wsavey(this%lensavy))
       this%Ly = 1.0_pfdp
       if(present(grid_size)) this%Ly = grid_size(2)
       !  Initialize FFT
       call ZFFTI( ny, this%wsavey)
       
       if (dim > 2) then
          !  FFT Storage
          nz=grid_shape(3)       
          this%nz = nz
          this%lensavz = 4*nz + 15
          this%normfact=nx*ny*nz             
          allocate(this%workhatz(nz))   !  complex transform
          allocate(this%wsavez(this%lensavz))
          this%Lz = 1.0_pfdp
          if(present(grid_size)) this%Lz = grid_size(3)
          !  Initialize FFT
          call ZFFTI( nz, this%wsavez)
       endif
    endif
    select case (this%dim)
    case (1)            
       allocate(this%wk_1d(nx))
    case (2)            
       allocate(this%wk_2d(nx,ny))
    case (3)            
       allocate(this%wk_3d(nx,ny,nz))
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)
    end select
  end subroutine fft_setup

  !>  Deallocate and destroy fft structures
  subroutine fft_destroy(this)
    class(pf_fft_t), intent(inout) :: this
    deallocate(this%workhatx)
    deallocate(this%wsavex)
    if (this%dim > 1) then
       deallocate(this%workhaty)
       deallocate(this%wsavey)
       if (this%dim > 2) then
          deallocate(this%workhatz)
          deallocate(this%wsavez)
       end if
    end if
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

  !>  Forward fft call
  subroutine fftf(this)
    class(pf_fft_t), intent(inout) :: this
    
    integer i,j,k
    
    select case (this%dim)       
    case (1)            
       call zfftf(this%nx, this%wk_1d, this%wsavex )
    case (2)
       do j = 1,this%ny
          this%workhatx =this%wk_2d(:,j)
          call zfftf(this%nx, this%workhatx, this%wsavex )
          this%wk_2d(:,j)=this%workhatx
       end do
       do i = 1,this%nx
          this%workhaty =this%wk_2d(i,:)
          call zfftf(this%ny, this%workhaty, this%wsavey )
          this%wk_2d(i,:)=this%workhaty
       end do
    case (3)            
       do k = 1,this%nz
          do j = 1,this%ny
             this%workhatx =this%wk_3d(:,j,k)
             call zfftf(this%nx, this%workhatx, this%wsavex )
             this%wk_3d(:,j,k)=this%workhatx
          end do
       end do
       do k = 1,this%nz
          do i = 1,this%nx
             this%workhaty =this%wk_3d(i,:,k)
             call zfftf(this%ny, this%workhaty, this%wsavey )
             this%wk_3d(i,:,k)=this%workhaty
          end do
       end do
       do j = 1,this%ny
          do i = 1,this%nx
             this%workhatz =this%wk_3d(i,j,:)
             call zfftf(this%nz, this%workhatz, this%wsavez)
             this%wk_3d(i,j,:)=this%workhatz
          end do
       end do
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)
    end select
  end subroutine fftf

  !  Backward FFT
  subroutine fftb(this)
    class(pf_fft_t), intent(inout) :: this
    
    integer i,j,k
    
    select case (this%dim)       
    case (1)
       this%wk_1d=this%wk_1d/this%normfact                                              
       call zfftb(this%nx, this%wk_1d, this%wsavex )
    case (2)
       this%wk_2d=this%wk_2d/this%normfact                                    
       do j = 1,this%ny
          this%workhatx =this%wk_2d(:,j)
          call zfftb(this%nx, this%workhatx, this%wsavex )
          this%wk_2d(:,j)=this%workhatx
       end do
       do i = 1,this%nx
          this%workhaty =this%wk_2d(i,:)
          call zfftb(this%ny, this%workhaty, this%wsavey )
          this%wk_2d(i,:)=this%workhaty
       end do
    case (3)            
       this%wk_3d=this%wk_3d/this%normfact                          
       do k = 1,this%nz
          do j = 1,this%ny
             this%workhatx =this%wk_3d(:,j,k)
             call zfftb(this%nx, this%workhatx, this%wsavex )
             this%wk_3d(:,j,k)=this%workhatx
          end do
       end do
       do k = 1,this%nz
          do i = 1,this%nx
             this%workhaty =this%wk_3d(i,:,k)
             call zfftb(this%ny, this%workhaty, this%wsavey )
             this%wk_3d(i,:,k)=this%workhaty
          end do
       end do
       do j = 1,this%ny
          do i = 1,this%nx
             this%workhatz =this%wk_3d(i,j,:)
             call zfftb(this%nz, this%workhatz, this%wsavez)
             this%wk_3d(i,j,:)=this%workhatz
          end do
       end do
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
    wk_f=wk_f*2.0_pfdp

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
    wk_f=wk_f*4.0_pfdp

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
    wk_f=wk_f*8.0_pfdp

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
   



