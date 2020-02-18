! MODULE: pf_my_sweeper
! IMEX Sweeper for 3-D  example in real space:
module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imex_sweeper
  use pf_mod_fftpackage
  use pf_mod_solutions
  use pf_mod_rutils
  use pf_mod_fftops
  implicit none

  !>  extend the IMEX sweeper
  type, extends(pf_imex_sweeper_t) :: my_sweeper_t
     integer ::     nx
     ! fft object and differentiaton matrices
     type(pf_fft_t), pointer :: fft_tool
     type(pf_fft_ops_t), pointer :: fft_ops
     real(pfdp), allocatable :: tmp(:,:,:) ! temp space for feval

     real(pfdp), pointer :: p_tmp(:,:,:) ! Useful pointer

     !  Useful for making f_eval and f_comp generic
     class(pf_ndarray_t), pointer :: f_encap,rhs_encap,y_encap 
     real(pfdp),      pointer :: yvec(:,:,:), rhsvec(:,:,:), fvec(:,:,:)     

   contains

     procedure :: f_eval       !  Computes the advection and diffusion terms
     procedure :: f_comp       !  Computes the advection and diffusion terms
     procedure :: initialize
     procedure :: initialize_tmp
     procedure :: destroy
  end type my_sweeper_t

contains

  !>  Routine to initialize tmp array for local storage
  subroutine initialize_tmp(this)
    class(my_sweeper_t), intent(inout) :: this

    integer :: istat
    allocate(this%tmp(this%nx,this%nx,this%nx),STAT=istat)
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)    
  end subroutine initialize_tmp

  
  ! The rest of the stuff is dimension independent
  include 'sweeper_include.f90'

  
  
end module pf_my_sweeper


