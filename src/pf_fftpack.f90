!!  Module of FFT based routines using fftpack
!
! This file is part of LIBPFASST.
!
!>  Module for using fftpack.  This is the default fft package when fftw is not used as
!!  indicated by setting the USE_FFTW=FALSE in Makefile.defauls
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
  subroutine fft_setup(this, grid_shape, ndim, grid_size)
    class(pf_fft_t), intent(inout) :: this
    integer,              intent(in   ) :: ndim
    integer,              intent(in   ) :: grid_shape(ndim)
    real(pfdp), optional, intent(in   ) :: grid_size(ndim)    
      
    integer     :: nx,ny,nz
    integer     :: i,j,k
    integer     :: ierr
    real(pfdp)  :: om
    this%ndim=ndim
    
    !  FFT Storage parameters
    nx=grid_SHAPE(1)
    this%nx = nx
    this%lensavx = 4*nx + 15
    this%normfact=REAL(nx,pfdp)
    
    allocate(this%workhatx(nx),stat=ierr)   !  complex transform
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)        
    allocate(this%wsavex(this%lensavx),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)        
    this%Lx = 1.0_pfdp
    if(present(grid_size)) this%Lx = grid_size(1)
    !  Initialize FFT
    call ZFFTI( nx, this%wsavex )
    
    if (ndim > 1) then
       !  FFT Storage
       ny=grid_SHAPE(2)       
       this%ny = ny
       this%lensavy = 4*ny + 15
       this%normfact=REAL(nx*ny,pfdp)
       allocate(this%workhaty(ny),stat=ierr)   !  complex transform
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
       allocate(this%wsavey(this%lensavy),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
       this%Ly = 1.0_pfdp
       if(present(grid_size)) this%Ly = grid_size(2)
       !  Initialize FFT
       call ZFFTI( ny, this%wsavey)
       
       if (ndim > 2) then
          !  FFT Storage
          nz=grid_SHAPE(3)       
          this%nz = nz
          this%lensavz = 4*nz + 15
          this%normfact=REAL(nx*ny*nz,pfdp)
          allocate(this%workhatz(nz),stat=ierr)   !  complex transform
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
          allocate(this%wsavez(this%lensavz),stat=ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
          this%Lz = 1.0_pfdp
          if(present(grid_size)) this%Lz = grid_size(3)
          !  Initialize FFT
          call ZFFTI( nz, this%wsavez)
       endif
    endif
    select case (this%ndim)
    case (1)            
       allocate(this%wk_1d(nx),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)           
    case (2)            
       allocate(this%wk_2d(nx,ny),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)           
    case (3)            
       allocate(this%wk_3d(nx,ny,nz),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%ndim)
    end select

    !  Assign wave numbers
    allocate(this%kx(nx),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)        
    om=two_pi/this%Lx                
    do i = 1, nx
       if (i <= nx/2+1) then
          this%kx(i) = om*REAL(i-1,pfdp)
       else
          this%kx(i) = om*REAL(-nx + i - 1,pfdp)
       end if
    end do

    if (ndim > 1) then
       allocate(this%ky(ny),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)           
       om=two_pi/this%Ly 
       do j = 1, ny
          if (j <= ny/2+1) then
             this%ky(j) = om*REAL(j-1,pfdp)
          else
             this%ky(j) = om*REAL(-ny + j - 1,pfdp)
          end if
       end do
    end if
    
    if (ndim > 2) then
       allocate(this%kz(nz),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)           
       om=two_pi / this%Lz 
       do k = 1,nz
          if (k <= nz/2+1) then
             this%kz(k) = om*REAL(k-1,pfdp)
          else
             this%kz(k) = om*REAL(-nz + k - 1,pfdp)
          end if
       end do
    end if

  end subroutine fft_setup

  !>  Deallocate and destroy fft structures
  subroutine fft_destroy(this)
    class(pf_fft_t), intent(inout) :: this
    deallocate(this%workhatx)
    deallocate(this%wsavex)
    deallocate(this%kx)   
    if (this%ndim > 1) then
       deallocate(this%workhaty)
       deallocate(this%wsavey)
       deallocate(this%ky)   
       if (this%ndim > 2) then
          deallocate(this%workhatz)
          deallocate(this%wsavez)
          deallocate(this%kz)   
       end if
    end if
    select case (this%ndim)
    case (1)            
       deallocate(this%wk_1d)
    case (2)            
       deallocate(this%wk_2d)
    case (3)            
       deallocate(this%wk_3d)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%ndim)
    end select

    !  Deallocate wave number arrays
    
  end subroutine fft_destroy

  !>  Forward fft call
  subroutine fftf(this)
    class(pf_fft_t), intent(inout) :: this
    
    integer i,j,k
    
    select case (this%ndim)       
    case (1)            
       call zfftf(this%nx, this%wk_1d, this%wsavex )
       this%wk_1d=this%wk_1d/this%normfact
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
       this%wk_2d=this%wk_2d/this%normfact
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
       this%wk_3d=this%wk_3d/this%normfact
       
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%ndim)
    end select
  end subroutine fftf

  !  Backward FFT
  subroutine fftb(this)
    class(pf_fft_t), intent(inout) :: this
    
    integer i,j,k
    
    select case (this%ndim)       
    case (1)
       this%wk_1d=this%wk_1d
       call zfftb(this%nx, this%wk_1d, this%wsavex )
    case (2)
       this%wk_2d=this%wk_2d
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
       this%wk_3d=this%wk_3d
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
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%ndim)
    end select
  end subroutine fftb
  

  subroutine interp_1d(this, yvec_c, fft_f,yvec_f,order)
    class(pf_fft_t), intent(inout) :: this
    real(pfdp), intent(inout), pointer :: yvec_f(:)
    real(pfdp), intent(in),    pointer :: yvec_c(:)
    type(pf_fft_t),intent(in), pointer :: fft_f
    integer, intent(in),optional :: order
    
    complex(pfdp),         pointer :: wk_f(:), wk_c(:)
    integer :: nx_f, nx_c,local_order
    real(pfdp) :: c1,c2,c3
    local_order=0
    if (present(order)) local_order = order

    nx_f = SIZE(yvec_f)
    nx_c = SIZE(yvec_c)
    if (nx_f .eq. nx_c) then
       yvec_f=yvec_c
       return
    end if
    
   select case (local_order)
   case (0)
      call this%get_wk_ptr(wk_c)
      call fft_f%get_wk_ptr(wk_f)
      
      wk_c=yvec_c
      
      call this%fftf()      !  internal forward fft call    
      call this%zinterp_1d(wk_c, wk_f)
      call fft_f%fftb()     !  internal inverse fft call
      
      yvec_f=REAL(wk_f,pfdp) !  grab the real part
   case (2)  !  This is for 2nd order Finite Difference in periodic domains
      yvec_f(1:nx_f-1:2)=yvec_c
      yvec_f(2:nx_f-2:2)=(yvec_c(1:nx_c-1)+yvec_c(2:nx_c))*0.5_pfdp
      yvec_f(nx_f)=(yvec_f(1) + yvec_f(nx_f-1))*0.5_pfdp
   case (4)
      if (nx_c .lt. 4) call pf_stop(__FILE__,__LINE__,'Coarse grid too small in interp_2d ',nx_c)
      yvec_f(1:nx_f-1:2)=yvec_c
      yvec_f(4:nx_f-4:2)=(-yvec_c(1:nx_c-3)+9.0_pfdp*(yvec_c(2:nx_c-2)+yvec_c(3:nx_c-1)) -yvec_c(4:nx_c) )*0.0625_pfdp
      yvec_f(2)=         (-yvec_f(nx_f-1)  +9.0_pfdp*(yvec_f(1)      + yvec_f(3)     )-yvec_f(5))*0.0625_pfdp         
      yvec_f(nx_f)=      (-yvec_f(3)       +9.0_pfdp*(yvec_f(1)      + yvec_f(nx_f-1))-yvec_f(nx_f-3))*0.0625_pfdp
      yvec_f(nx_f-2)=    (-yvec_f(1)       +9.0_pfdp*(yvec_f(nx_f-1) + yvec_f(nx_f-3))-yvec_f(nx_f-5))*0.0625_pfdp
   case (6)
      if (nx_c .lt. 6) call pf_stop(__FILE__,__LINE__,'Coarse grid too small in interp_1d ',nx_c)
      c1=  3.0_pfdp/256.0_pfdp
      c2=-25.0_pfdp/256.0_pfdp
      c3= 75.0_pfdp/128.0_pfdp
      yvec_f(1:nx_f-1:2)=yvec_c
      yvec_f(6:nx_f-6:2)=c1*(yvec_c(1:nx_c-5)+yvec_c(6:nx_c)) +c2*(yvec_c(2:nx_c-4)+yvec_c(5:nx_c-1)) +c3*(yvec_c(3:nx_c-3)+yvec_c(4:nx_c-2))
      yvec_f(2)=c1*(yvec_f(nx_f-3)+yvec_f(7)) +c2*(yvec_f(nx_f-1)+yvec_f(5)) +c3*(yvec_f(1)+yvec_f(3))         
      yvec_f(4)=c1*(yvec_f(nx_f-1)+yvec_f(9)) +c2*(yvec_f(1)+yvec_f(7)) +c3*(yvec_f(3)+yvec_f(5))         
      yvec_f(nx_f)=c1*(yvec_f(5)+yvec_f(nx_f-5)) +c2*(yvec_f(3)+yvec_f(nx_f-3)) +c3*(yvec_f(1)+yvec_f(nx_f-1))         
      yvec_f(nx_f-2)=c1*(yvec_f(3)+yvec_f(nx_f-7)) +c2*(yvec_f(1)+yvec_f(nx_f-5)) +c3*(yvec_f(nx_f-1)+yvec_f(nx_f-3))         
      yvec_f(nx_f-4)=c1*(yvec_f(1)+yvec_f(nx_f-9)) +c2*(yvec_f(nx_f-1)+yvec_f(nx_f-7)) +c3*(yvec_f(nx_f-3)+yvec_f(nx_f-5))         
   case DEFAULT
      call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',local_order)
   end select
    
  end subroutine interp_1d
  subroutine interp_2d(this, yvec_c, fft_f,yvec_f)
    class(pf_fft_t), intent(inout) :: this
    real(pfdp), intent(inout), pointer :: yvec_f(:,:)
    real(pfdp), intent(in),    pointer :: yvec_c(:,:)
    type(pf_fft_t),intent(in), pointer:: fft_f

    complex(pfdp),         pointer :: wk_f(:,:), wk_c(:,:)
    integer :: nx_f(2), nx_c(2)

    nx_f = SHAPE(yvec_f)
    nx_c = SHAPE(yvec_c)

    if (nx_f(1) .eq. nx_c(1) .and. nx_f(2) .eq. nx_c(2)) then
       yvec_f=yvec_c
       return
    end if

    call this%get_wk_ptr(wk_c)
    call fft_f%get_wk_ptr(wk_f)
    
    wk_c=yvec_c

    call this%fftf()       !  internal forward fft call     
    call this%zinterp_2d(wk_c, wk_f)
    call fft_f%fftb()     !  internal inverse fft call

    yvec_f=REAL(wk_f,pfdp)  !  grab the real part
    
  end subroutine interp_2d
  subroutine interp_3d(this, yvec_c, fft_f,yvec_f)
    class(pf_fft_t), intent(inout) :: this
    real(pfdp), intent(inout),  pointer :: yvec_f(:,:,:)
    real(pfdp), intent(inout),  pointer :: yvec_c(:,:,:)
    type(pf_fft_t),intent(in),  pointer :: fft_f
    complex(pfdp),         pointer :: wk_f(:,:,:), wk_c(:,:,:)
    integer :: nx_f(3), nx_c(3)

    nx_f = SHAPE(yvec_f)
    nx_c = SHAPE(yvec_c)

    if (nx_f(1) .eq. nx_c(1) .and. nx_f(2) .eq. nx_c(2) .and. nx_f(3) .eq. nx_c(3)) then
       yvec_f=yvec_c
       return
    end if

    call this%get_wk_ptr(wk_c)
    call fft_f%get_wk_ptr(wk_f)
    
    wk_c=yvec_c

    call this%fftf()      !  internal forward fft call    
    call this%zinterp_3d(wk_c, wk_f)
    call fft_f%fftb()     !  internal inverse fft call

    yvec_f=REAL(wk_f,pfdp)  !  grab the real part
  end subroutine interp_3d

  !>  Interpolate from coarse  level to fine in spectral space
  subroutine zinterp_1d(this, yhat_c, yhat_f)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp),   pointer,intent(inout) :: yhat_f(:) 
    complex(pfdp),   pointer,intent(in) :: yhat_c(:)

    integer :: nx_f, nx_c

    nx_f = SIZE(yhat_f)
    nx_c = SIZE(yhat_c)
    if (nx_f .eq. nx_c) then
       yhat_f=yhat_c
       return
    end if
    yhat_f = 0.0_pfdp
    yhat_f(1:nx_c/2) = yhat_c(1:nx_c/2)
    yhat_f(nx_f-nx_c/2+2:nx_f) = yhat_c(nx_c/2+2:nx_c)
    
  end subroutine zinterp_1d

  subroutine zinterp_2d(this, yhat_c, yhat_f)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp),   pointer,intent(inout) :: yhat_f(:,:) 
    complex(pfdp),   pointer,intent(in) :: yhat_c(:,:)

    integer :: nf1,nf2,nc1,nc2
    integer :: nx_f(2), nx_c(2)

    nx_f = SHAPE(yhat_f)
    nx_c = SHAPE(yhat_c)

    if (nx_f(1) .eq. nx_c(1) .and. nx_f(2) .eq. nx_c(2)) then
       yhat_f=yhat_c
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
    complex(pfdp),   pointer,intent(inout) :: yhat_f(:,:,:) 
    complex(pfdp),   pointer,intent(in) :: yhat_c(:,:,:)

    integer :: nf1,nf2,nf3,nc1,nc2,nc3
    integer :: nx_f(3), nx_c(3)

    nx_f = SHAPE(yhat_f)
    nx_c = SHAPE(yhat_c)
    if (nx_f(1) .eq. nx_c(1) .and. nx_f(2) .eq. nx_c(2) .and. nx_f(3) .eq. nx_c(3))  then
       yhat_f=yhat_c
       return
    end if
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
   



