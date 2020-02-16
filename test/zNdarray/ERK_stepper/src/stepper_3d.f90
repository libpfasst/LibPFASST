!> Exponential RK stepper module
module pf_my_stepper
  use pf_mod_dtype
  use pf_mod_zndarray
  use pf_mod_erkstepper
  use pf_mod_fftpackage
  use pf_mod_solutions
  use phi_mod
  use pf_mod_zutils
  use pf_mod_fftops
  implicit none
  
  ! Define the derived stepper type
  type, extends(pf_erk_stepper_t) :: my_stepper_t
     integer ::     nx
     ! phi storage space
     complex(pfdp),  allocatable :: A_phi(:,:,:,:)          ! coefficients for A vector
     complex(pfdp),  allocatable :: b_phi(:,:,:,:)          ! coefficients for b vector
     complex(pfdp),  allocatable :: d_phi(:,:,:,:)          ! coefficients for d vector
     real(pfdp)                  :: dt_phi              ! dt used to store coefficients
     ! fft object and differentiaton matrices
     type(pf_fft_t), pointer :: fft_tool
     type(pf_fft_ops_t), pointer :: fft_ops     
     complex(pfdp), allocatable :: tmp(:,:,:) ! temp space for feval
     complex(pfdp), pointer :: p_tmp(:,:,:)   ! a useful pointer
     !  Useful for making f_eval and f_comp generic
     class(pf_zndarray_t), pointer :: f_encap,rhs_encap,y_encap
     complex(pfdp),      pointer :: yvec(:,:,:), rhsvec(:,:,:), fvec(:,:,:)          
   contains
     
     procedure :: f_eval       !  Computes the advection and diffusion terms
     procedure :: compA        !  applies phi functions for tableau matrix A(i,j)
     procedure :: compB        !  applied phi functions for tableau vector B(i)
     procedure :: compD        !  applies phi functions for tableau vector D(i) 
     procedure :: initialize
     procedure :: destroy
     procedure, private :: initRKCoeff
     procedure, private :: mult_dphi
     procedure, private :: mult_bphi
     procedure, private :: mult_Aphi
     procedure, private :: initialize_tmp
  end type my_stepper_t
  
contains
  subroutine mult_Aphi(this,i,y,f)
    class(my_stepper_t), intent(inout) :: this
    integer, intent(in) :: i
    class(pf_encap_t), intent(in) :: y
    class(pf_encap_t), intent(in) :: f

    complex(pfdp),        pointer :: yvec(:,:,:), fvec(:,:,:)
    class(pf_zndarray_t), pointer :: y_encap,f_encap
    y_encap => cast_as_zndarray(y)
    call y_encap%get_array(yvec)  
    f_encap => cast_as_zndarray(f)
    call f_encap%get_array(fvec)  

    yvec=this%A_phi(:,:,:, i)*fvec
  end subroutine mult_Aphi

  subroutine mult_dphi(this,i,y,f)
    class(my_stepper_t), intent(inout) :: this
    integer, intent(in) :: i
    class(pf_encap_t), intent(in) :: y
    class(pf_encap_t), intent(in) :: f

    complex(pfdp),        pointer :: yvec(:,:,:), fvec(:,:,:)
    class(pf_zndarray_t), pointer :: y_encap,f_encap
    y_encap => cast_as_zndarray(y)
    call y_encap%get_array(yvec)  
    f_encap => cast_as_zndarray(f)
    call f_encap%get_array(fvec)  

    yvec=this%d_phi(:,:,:, i)*fvec
  end subroutine mult_dphi

  subroutine mult_bphi(this,i,y,f)
    class(my_stepper_t), intent(inout) :: this
    integer, intent(in) :: i
    class(pf_encap_t), intent(in) :: y
    class(pf_encap_t), intent(in) :: f

    complex(pfdp),        pointer :: yvec(:,:,:), fvec(:,:,:)
    class(pf_zndarray_t), pointer :: y_encap,f_encap
    y_encap => cast_as_zndarray(y)
    call y_encap%get_array(yvec)  
    f_encap => cast_as_zndarray(f)
    call f_encap%get_array(fvec)  

    yvec=this%b_phi(:,:,:, i)*fvec
  end subroutine mult_bphi
  
  subroutine initialize_tmp(this)
    class(my_stepper_t), intent(inout) :: this

    integer :: istat,nx
    nx=this%nx
    allocate(this%tmp(nx,nx,nx),STAT=istat)
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)

    ! allocate space for phi functions
    allocate(this%A_phi(nx,nx,nx,this%nnz_A),STAT=istat)
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
    allocate(this%b_phi(nx,nx,nx, this%nstages),STAT=istat)
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
    allocate(this%d_phi(nx,nx,nx, this%nstages),STAT=istat)
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
    
  end subroutine initialize_tmp

  subroutine initRKCoeff(this, opL, dt)
    
    ! Arguments
    class(my_stepper_t), intent(inout) :: this
    complex(pfdp), intent(in) :: opL(:,:,:)
    real(pfdp), intent(in) :: dt
    
    ! Variables
    integer :: i, j, k, num_phi, A_ij
    complex(pfdp), allocatable :: phi(:,:,:,:)
    real(pfdp) :: c1, c2
    num_phi = size(this%A,2)
    allocate(phi(num_phi, this%nx, this%nx, this%nx))
    
    ! init A Tablaeu coefficients
    this%A_phi = (0.0_pfdp, 0.0_pfdp)
    do i = 1, this%nstages - 1
       do j = 1, i
          A_ij = this%AF(j,i)
          if(A_ij .ne. 0) then ! A(i,j) is nonzero
             do k = 1, num_phi
                c1 = this%A(1, k, j, i)
                c2 = this%A(2, k, j, i)
                if(c1 .ne. 0.0_pfdp) then ! A(i,j) contains \varphi_{k-1}(c2 opL)
                   call phi_zmatrix3d(dt*c2*opL, num_phi-1, phi)
                   this%A_phi(:,:,:, A_ij) = this%A_phi(:,:,:, A_ij)+dt*c1*phi(k,:, :,:)
                endif
             enddo
          endif
       enddo
    enddo
         
    ! init b vector coefficients
    this%b_phi = (0.0_pfdp, 0.0_pfdp)
    do i = 1, this%nstages
       do k = 1, num_phi
          c1 = this%b(1, k, i)
          c2 = this%b(2, k, i)
          if(c1 .ne. 0.0_pfdp) then ! b(i) contains \varphi_{k-1}(c2 opL)
             call phi_zmatrix3d(dt*c2*opL, num_phi - 1, phi)
             this%b_phi(:,:,:, i) = this%b_phi(:,:,:, i) + dt*c1*phi(k,:, :,:)
          endif
       enddo
    enddo

    ! init d vector coefficients        
    do i = 1, this%nstages
       c1 = this%d(1, i)
       c2 = this%d(2, i)
       if(c1 .ne. 0.0_pfdp) then ! d(i) contains \exp(c2 opL)
          call phi_zmatrix3d(dt*c2*opL, num_phi-1, phi)
          this%d_phi(:,:,:, i) = c1*phi(1,:,:,:)
       endif
    enddo
    deallocate(phi)
  end subroutine initRKCoeff

  ! The rest of the stuff is dimension independent
  include 'stepper_include.f90'

  
end module pf_my_stepper
