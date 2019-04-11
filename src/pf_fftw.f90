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
  
  implicit none
  
  include 'fftw3.f03'
  
  !>  Variables and storage for FFTW
  type,extends(pf_fft_abs_t) :: pf_fft_t
     type(c_ptr) :: ffft, ifft                       !! fftw pointers
     real(pfdp) :: normfact                          !! normalization factor
   contains
     procedure :: fft_setup
     procedure :: fft_destroy
     procedure :: fftf
     procedure :: fftb
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
       
       this%ffft = fftw_plan_dft_1d(nx, &
            this%wk_1d, this%wk_1d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifft = fftw_plan_dft_1d(nx, &
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

       this%ffft = fftw_plan_dft_2d(nx,ny, &
            this%wk_2d, this%wk_2d, FFTW_FORWARD, FFTW_ESTIMATE)
       this%ifft = fftw_plan_dft_2d(nx,ny, &
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

  !>  Routine to take inverse or backward FFT
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

end module pf_mod_fftpackage



