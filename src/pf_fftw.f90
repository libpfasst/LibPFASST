!!  Module of FFT based routines using fftw
!!
! This file is part of LIBPFASST.
!
!>  Module for providing FFTs based on fftw
!!  To use this module, fftw must be installed.
!!  This can be done by the libpfasst make system with the command
!!  make fftw3
!!  You also have to specify USE_FFTW=TRUE in Makefile.defaults and
!!  make clean is recommended before remake
module pf_mod_fftpackage
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_fft_abs
  use, intrinsic :: iso_c_binding  
  implicit none
  include 'fftw3.f03'
  
  !>  Variables and storage for FFTW
  type,extends(pf_fft_abs_t) :: pf_fft_t
     type(c_ptr) :: ffftw, ifftw         !! fftw pointers
     real(pfdp) :: normfact              !! normalization factor
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
  subroutine fft_setup(this, grid_shape, ndim, grid_size)
    class(pf_fft_t), intent(inout) :: this
    integer,             intent(in   ) :: ndim
    integer,             intent(in   ) :: grid_shape(ndim)    
    real(pfdp), optional, intent(in   ) :: grid_size(ndim)    

    integer     :: nx,ny,nz
    integer     :: i,j,k
    integer     :: ierr
    real(pfdp)  :: om
    type(c_ptr) :: wk

    this%ndim=ndim
    nx=grid_SHAPE(1)
    this%nx = nx

    ! Defaults for grid_size
    this%Lx = 1.0_pfdp
    this%Ly = 1.0_pfdp
    this%Lz = 1.0_pfdp

    select case (ndim)
    case (1)            
       if(present(grid_size)) this%Lx = grid_size(1)
       this%normfact=REAL(nx,pfdp)
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
       
       ny=grid_SHAPE(2)          
       this%ny = ny
       this%normfact=REAL(nx*ny,pfdp)
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
    
       ny=grid_SHAPE(2)          
       nz=grid_SHAPE(3)          
       this%ny = ny
       this%nz = nz
       this%normfact=REAL(nx*ny*nz,pfdp)          
       wk = fftw_alloc_complex(int(nx*ny*nz, c_size_t))
       call c_f_pointer(wk, this%wk_3d, [nx,ny,nz])
       this%ffftw = fftw_plan_dft_3d(nx,ny,nz, &
            this%wk_3d, this%wk_3d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifftw = fftw_plan_dft_3d(nx,ny,nz, &
            this%wk_3d, this%wk_3d, FFTW_BACKWARD, FFTW_ESTIMATE)
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

  !>  Destroy the package
  subroutine fft_destroy(this)
    class(pf_fft_t), intent(inout) :: this

    type(c_ptr) :: wk    
    call fftw_destroy_plan(this%ffftw)
    call fftw_destroy_plan(this%ifftw)
    select case (this%ndim)
    case (1)
       call fftw_free(c_loc(this%wk_1d))
       deallocate(this%kx)
    case (2)            
       call fftw_free(c_loc(this%wk_2d))
       !  Deallocate wave number arrays
       deallocate(this%kx)   
       deallocate(this%ky)   
    case (3)            
       call fftw_free(c_loc(this%wk_3d))

       !  Deallocate wave number arrays
       deallocate(this%kx)   
       deallocate(this%ky)   
       deallocate(this%kz)   
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%ndim)
    end select
      
  end subroutine fft_destroy
  
  !>  Routine to take foreward FFT
  subroutine fftf(this)
    class(pf_fft_t), intent(inout) :: this
    select case (this%ndim)       
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
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%ndim)       
    end select
  end subroutine fftf

  !>  Routine to take inverse or backward FFT
  subroutine fftb(this)
    class(pf_fft_t), intent(inout) :: this
    !  Normalize the fft
    select case (this%ndim)       
    case (1)
       call fftw_execute_dft(this%ifftw, this%wk_1d, this%wk_1d)
    case (2)
       call fftw_execute_dft(this%ifftw, this%wk_2d, this%wk_2d)
    case (3)            
       call fftw_execute_dft(this%ifftw, this%wk_3d, this%wk_3d)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%ndim)
    end select
  end subroutine fftb
  
  subroutine interp_1d(this, yvec_c, fft_f,yvec_f,order)
    class(pf_fft_t), intent(inout) :: this
    real(pfdp), intent(inout),  pointer :: yvec_f(:)
    real(pfdp), intent(in),     pointer :: yvec_c(:)
    type(pf_fft_t),intent(in),  pointer :: fft_f
    integer, intent(in),optional :: order

    complex(pfdp),         pointer :: wk_f(:), wk_c(:)
    integer :: nx_f, nx_c,local_order
    real :: c1,c2,c3
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
      
      call this%fftf()       !  internal forward fft call    
      call this%zinterp_1d(wk_c, wk_f)
      call fft_f%fftb()     !  internal inverse fft call
      
      yvec_f=REAL(wk_f,pfdp)     !  grab the real part
   case (2)  !  This is for 2nd order Finite Difference in periodic domains
      yvec_f(1:nx_f-1:2)=yvec_c
      yvec_f(2:nx_f-2:2)=(yvec_c(1:nx_c-1)+yvec_c(2:nx_c))*0.5_pfdp
      yvec_f(nx_f)=(yvec_f(1) + yvec_f(nx_f-1))*0.5_pfdp
   case (4)
      if (nx_c .lt. 4) call pf_stop(__FILE__,__LINE__,'Coarse grid too small in interp_1d ',nx_c)
      c2=-1.0_pfdp/16.0_pfdp
      c1= 9.0_pfdp/16.0_pfdp
      yvec_f(1:nx_f-1:2)=yvec_c
      yvec_f(2:nx_f:2)=c2*(cshift(yvec_c,-1)+cshift(yvec_c,2))+c1*(yvec_c+cshift(yvec_c,1))      
   case (6)
      if (nx_c .lt. 6) call pf_stop(__FILE__,__LINE__,'Coarse grid too small in interp_1d ',nx_c)
      c3=  3.0_pfdp/256.0_pfdp
      c2=-25.0_pfdp/256.0_pfdp
      c1= 75.0_pfdp/128.0_pfdp
      yvec_f(1:nx_f-1:2)=yvec_c
      yvec_f(2:nx_f:2)=c3*(cshift(yvec_c,-2)+cshift(yvec_c,3))+c2*(cshift(yvec_c,-1)+cshift(yvec_c,2))+c1*(yvec_c+cshift(yvec_c,1))
   case DEFAULT
      call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',local_order)
   end select
    
  end subroutine interp_1d
  subroutine interp_2d(this, yvec_c, fft_f,yvec_f)
    class(pf_fft_t), intent(inout) :: this
    real(pfdp), intent(inout), pointer :: yvec_f(:,:)
    real(pfdp), intent(in), pointer    :: yvec_c(:,:)
    type(pf_fft_t),intent(in), pointer :: fft_f

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
    call this%fftf()                  !  internal forward fft call    
    call this%zinterp_2d(wk_c, wk_f)  !  interpolate in spectral space
    call fft_f%fftb()                 !  internal inverse fft call
    yvec_f=REAL(wk_f,pfdp)                 !  grab the real part
    
  end subroutine interp_2d
  subroutine interp_3d(this, yvec_c, fft_f,yvec_f)
    class(pf_fft_t), intent(inout) :: this
    real(pfdp), intent(inout),  pointer :: yvec_f(:,:,:)
    real(pfdp), intent(in),  pointer    :: yvec_c(:,:,:)
    type(pf_fft_t), intent(in), pointer :: fft_f

    complex(pfdp),  pointer :: wk_f(:,:,:), wk_c(:,:,:)
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
    call this%fftf()                 !  internal forward fft call    
    call this%zinterp_3d(wk_c, wk_f) ! interpolate in spectral space
    call fft_f%fftb()                ! internal inverse fft call
    yvec_f=REAL(wk_f,pfdp)                ! grab the real part

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
    yhat_f(nx_c/2+1) = yhat_c(nx_c/2+1)*0.5_pfdp
    yhat_f(nx_f-nx_c/2+1) = yhat_c(nx_c/2+1)*0.5_pfdp
    
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

    yhat_f = 0.0_pfdp
  
    nx_f = SHAPE(yhat_f)
    nx_c = SHAPE(yhat_c)
    
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



