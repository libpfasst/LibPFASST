!
! This file is part of LIBPFASST.
!
!>  Module for providing FFTs based on fftw
!!  To use this module, fftw must be installed.
!!  This can be done by the libpfasst make system with the comment
!!   > make fftw3
module pf_mod_fftpackage
  use pf_mod_dtype
  use pf_mod_utils
  
  implicit none
  
  include 'fftw3.f03'
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp
  
  !>  Define the fft package
  type :: pf_fft_t
     type(c_ptr) :: ffft, ifft                       !! fftw pointers
     integer ::  dim                                 !! spatial number of dimensions
     integer ::  nx,ny,nz                            !! grid sizes
     real(pfdp) :: Lx, Ly, Lz                        !! domain size
     real(pfdp) :: normfact                          !! normalization factor
     complex(pfdp), pointer :: wk_1d(:)              !! work space
     complex(pfdp), pointer :: wk_2d(:,:)            !! work space
     complex(pfdp), pointer :: wk_3d(:,:,:)          !! work space               
   contains
     procedure :: fft_setup
     procedure :: fft_destroy
     procedure :: fftf
     procedure :: fftb
     procedure :: get_wk_ptr_1d
     procedure :: get_wk_ptr_2d
     procedure :: get_wk_ptr_3d
     procedure :: make_lap_1d
     procedure :: make_lap_2d
     procedure :: make_lap_3d
     procedure :: make_deriv_1d
     procedure :: make_deriv_2d
     procedure :: make_deriv_3d                    
  end type pf_fft_t 
  
contains

  !> Routines to return the pointer to the work variable
  function get_wk_ptr_1d(this) result(wk)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), pointer :: wk(:)              ! work space
    wk=>this%wk_1d
  end function get_wk_ptr_1d
  function get_wk_ptr_2d(this) result(wk)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), pointer :: wk(:,:)              ! work space
    wk=>this%wk_2d
  end function get_wk_ptr_2d
  function get_wk_ptr_3d(this) result(wk)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), pointer :: wk(:,:,:)              ! work space
    wk=>this%wk_3d
  end function get_wk_ptr_3d
    
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
    
    if(present(grid_size)) this%Lx = grid_size(1)
    if(present(grid_size)) this%Ly = grid_size(2)
    if(present(grid_size)) this%Lz = grid_size(3)

    select case (dim)
    case (1)            
       this%normfact=real(nx,pfdp)
       wk = fftw_alloc_complex(int(nx, c_size_t))
       call c_f_pointer(wk, this%wk_1d, [nx])          
       
       this%ffft = fftw_plan_dft_1d(nx, &
            this%wk_1d, this%wk_1d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifft = fftw_plan_dft_1d(nx, &
            this%wk_1d, this%wk_1d, FFTW_BACKWARD, FFTW_ESTIMATE)
    case (2)  
       ny=grid_shape(2)          
       this%ny = ny
       this%normfact=real(nx*ny,pfdp)
       ! create in-place, complex fft plans
       wk = fftw_alloc_complex(int(nx*ny, c_size_t))
       call c_f_pointer(wk, this%wk_2d, [nx,ny])

       this%ffft = fftw_plan_dft_2d(nx,ny, &
            this%wk_2d, this%wk_2d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifft = fftw_plan_dft_2d(nx,ny, &
            this%wk_2d, this%wk_2d, FFTW_BACKWARD, FFTW_ESTIMATE)
    case (3)
       ny=grid_shape(2)          
       nz=grid_shape(3)          
       this%ny = ny
       this%nz = nz
       this%normfact=real(nx*ny*nz,pfdp)          
       wk = fftw_alloc_complex(int(nx*ny*nz, c_size_t))
       call c_f_pointer(wk, this%wk_3d, [nx,ny,nz])
       this%ffft = fftw_plan_dft_3d(nx,ny,nz, &
            this%wk_3d, this%wk_3d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifft = fftw_plan_dft_3d(nx,ny,nz, &
            this%wk_3d, this%wk_3d, FFTW_BACKWARD, FFTW_ESTIMATE)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)
    end select
  end subroutine fft_setup

  !>  Destroy the package
  subroutine fft_destroy(this)
    class(pf_fft_t), intent(inout) :: this
    call fftw_destroy_plan(this%ffft)
    call fftw_destroy_plan(this%ifft)
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
       call fftw_execute_dft(this%ffft, this%wk_1d, this%wk_1d)
    case (2)            
       call fftw_execute_dft(this%ffft, this%wk_2d, this%wk_2d)
    case (3)            
       call fftw_execute_dft(this%ffft, this%wk_3d, this%wk_3d)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)       
    end select
  end subroutine fftf

  !>  Routine to take inverser or backward FFT
  subroutine fftb(this)
    class(pf_fft_t), intent(inout) :: this
    !  Normalize the fft
    select case (this%dim)       
    case (1)
       this%wk_1d=this%wk_1d/this%normfact          
       call fftw_execute_dft(this%ifft, this%wk_1d, this%wk_1d)
    case (2)
       this%wk_2d=this%wk_2d/this%normfact          
       call fftw_execute_dft(this%ifft, this%wk_2d, this%wk_2d)
    case (3)            
       this%wk_3d=this%wk_3d/this%normfact
       call fftw_execute_dft(this%ifft, this%wk_3d, this%wk_3d)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',this%dim)
    end select
  end subroutine fftb

  !>  Routines to construct spectral differential operators
  subroutine make_lap_1d(this, lap)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: lap(:)
    
    integer     :: i,nx
    real(pfdp)  :: kx
    
    nx=this%nx
    do i = 1, nx
       if (i <= nx/2+1) then
          kx = two_pi / this%Lx * dble(i-1)
       else
          kx = two_pi / this%Lx * dble(-nx + i - 1)
       end if
       lap(i) = -kx**2
    end do
  end subroutine make_lap_1d
  subroutine make_deriv_1d(this, ddx)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: ddx(:)
    
    integer     :: i,nx
    real(pfdp)  :: kx
    
    nx=this%nx
    do i = 1, nx
       if (i <= nx/2+1) then
          kx = two_pi /this%Lx * dble(i-1)
       else
          kx = two_pi / this%Lx * dble(-nx + i - 1)
       end if
       
       ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx
    end do
  end subroutine make_deriv_1d
     
  subroutine make_lap_2d(this, lap)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: lap(:,:)
    
    integer     :: i,j,nx,ny
    real(pfdp)  :: kx,ky


    nx=this%nx
    ny=this%ny

    do j = 1, ny
       if (j <= ny/2+1) then
          ky = two_pi / this%Ly * dble(j-1)
       else
          ky = two_pi / this%Ly * dble(-ny + j - 1)
       end if
       do i = 1, nx
          if (i <= nx/2+1) then
             kx = two_pi / this%Lx * dble(i-1)
          else
             kx = two_pi / this%Lx * dble(-nx + i - 1)
          end if
          
          lap(i,j) = -(kx**2+ky**2)
       end do
    end do
  end subroutine make_lap_2d

  subroutine make_deriv_2d(this, deriv,dir)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: deriv(:,:)
    integer, intent(in) :: dir
    
    integer     :: i,j,nx,ny
    real(pfdp)  :: kx,ky
    
    nx=this%nx
    ny=this%ny
    
    do j = 1, ny
       if (j <= ny/2+1) then
          ky = two_pi / this%Ly * dble(j-1)
       else
          ky = two_pi / this%Ly * dble(-ny + j - 1)
       end if
       do i = 1, nx
          if (i <= nx/2+1) then
             kx = two_pi / this%Lx * dble(i-1)
          else
             kx = two_pi / this%Lx * dble(-nx + i - 1)
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
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: lap(:,:,:)
    
    integer     :: i,j,k,nx,ny,nz
    real(pfdp)  :: kx,ky,kz
    
    nx=this%nx
    ny=this%ny
    nz=this%nz
    do k = 1,nz
       if (k <= nz/2+1) then
          kz = two_pi / this%Lz * dble(k-1)
       else
          kz = two_pi / this%Lz * dble(-nz + k - 1)
       end if
       do j = 1, ny
          if (j <= ny/2+1) then
             ky = two_pi / this%Ly * dble(j-1)
          else
             ky = two_pi / this%Ly * dble(-ny + j - 1)
          end if
          do i = 1, nx
             if (i <= nx/2+1) then
                kx = two_pi / this%Lx * dble(i-1)
             else
                kx = two_pi / this%Lx * dble(-nx + i - 1)
             end if
             lap(i,j,k) = -(kx**2+ky**2+kz**2)
          end do
       end do
    end do
  end subroutine make_lap_3d

  subroutine make_deriv_3d(this, deriv,dir)
    class(pf_fft_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: deriv(:,:,:)
    integer, intent(in) :: dir
    
    integer     :: i,j,k,nx,ny,nz
    real(pfdp)  :: kx,ky,kz

    nx=this%nx
    ny=this%ny
    nz=this%nz       
    
    select case (dir)
    case (1)  
       do i = 1, nx
          if (i <= nx/2+1) then
             kx = two_pi / this%Lx * dble(i-1)
          else
             kx = two_pi / this%Lx * dble(-nx + i - 1)
          end if
          deriv(i,:,:) = (0.0_pfdp,1.0_pfdp)*kx
       end do
       
    case (2)
       do j = 1, ny
          if (j <= ny/2+1) then
             ky = two_pi / this%Ly * dble(j-1)
          else
             ky = two_pi / this%Ly * dble(-ny + j - 1)
          end if
          deriv(:,j,:) = (0.0_pfdp,1.0_pfdp)*ky
       end do
          
    case (3)
       do k = 1, nz
          if (k <= nz/2+1) then
             kz = two_pi / this%Lz * dble(k-1)
          else
             kz = two_pi / this%Lz * dble(-nz + k - 1)
          end if
          deriv(:,:,k) = (0.0_pfdp,1.0_pfdp)*kz
       end do
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',dir)
    end select
  end subroutine make_deriv_3d

end module pf_mod_fftpackage



