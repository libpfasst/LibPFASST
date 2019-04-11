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
  
    
  end module pf_mod_fftpackage
   



