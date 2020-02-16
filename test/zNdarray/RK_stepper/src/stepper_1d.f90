!>  RK stepper
module pf_my_stepper
  use pf_mod_dtype
  use pf_mod_zndarray
  use pf_mod_rkstepper
  use pf_mod_fftpackage
  use pf_mod_solutions
  use pf_mod_zutils
  use pf_mod_fftops
  implicit none
  
  
  ! Define the derived stepper type
  type, extends(pf_ark_stepper_t) :: my_stepper_t
     integer ::     nx
     ! fft object and differentiaton matrices
     type(pf_fft_t), pointer :: fft_tool
     type(pf_fft_ops_t), pointer :: fft_ops
     complex(pfdp), allocatable :: tmp(:)                ! Temp space for feval
     complex(pfdp), pointer :: p_tmp(:)   ! a useful pointer
     
     !  Useful for making routines  generic by dimension
     class(pf_zndarray_t), pointer :: f_encap,rhs_encap,y_encap 
     complex(pfdp),        pointer :: yvec(:), rhsvec(:), fvec(:)

   contains

     procedure :: f_eval       !  Computes the advection and diffusion terms
     procedure :: f_comp       !  Computes the advection and diffusion terms
     procedure :: initialize
     procedure :: initialize_tmp
     procedure :: destroy
     
  end type my_stepper_t
  
contains
  subroutine initialize_tmp(this)
    class(my_stepper_t), intent(inout) :: this
    
    integer :: istat
    allocate(this%tmp(this%nx),STAT=istat)
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
    
  end subroutine initialize_tmp
  
  ! The rest of the stuff is dimension independent
  include 'stepper_include.f90'
  
end module pf_my_stepper
