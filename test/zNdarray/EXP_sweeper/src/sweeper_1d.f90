! MODULE: pf_my_sweeper
! Exponential Sweeper for 1-D  examples
module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_zndarray
  use pf_mod_fexp_sweeper
  use phi_mod
  use pf_mod_fftpackage
  use pf_mod_solutions
  use pf_mod_zutils
  use pf_mod_fftops
  implicit none

  !>  extend the exponential sweeper
  type, extends(pf_fexp_sweeper_t) :: my_sweeper_t
     integer ::     nx
     integer ::     nnodes
     ! fft object and differentiaton matrices
     type(pf_fft_t), pointer :: fft_tool
     type(pf_fft_ops_t), pointer :: fft_ops
     complex(pfdp), allocatable :: tmp(:) ! temp space for feval
     
     complex(pfdp), pointer :: p_tmp(:) ! Useful pointer
     
     !  Useful for making f_eval and f_comp generic
     class(pf_zndarray_t), pointer :: f_encap,rhs_encap,y_encap 
     complex(pfdp),        pointer :: yvec(:), rhsvec(:), fvec(:)
     ! phi storage and scratch space
     complex(pfdp), allocatable :: P_sweep(:,:,:)  ! phi_0 and phi_1 for sweeper
     complex(pfdp), allocatable :: W_sweep(:,:,:)	! W for residual
     real(pfdp) :: dt_sweep    = real(0.0, pfdp)   ! stepsize used for sweep
     real(pfdp) :: dt_residual = real(0.0, pfdp)   ! stepsize used for residual
   contains
     
     procedure :: f_eval      	      !  Computes the advection and diffusion terms
     procedure :: expResidualSubstep   
     procedure :: expSweepSubstep      
     procedure :: initWSweep          !  Computes W_sweep and P_sweep
     procedure :: initialize
     procedure :: initialize_tmp
     procedure :: destroy
     procedure,private :: P_axpy
     procedure,private :: P_daxpy
     procedure,private :: W_axpy
  end type my_sweeper_t

contains
  subroutine initialize_tmp(this)
    class(my_sweeper_t), intent(inout) :: this
    integer :: istat
    allocate(this%tmp(this%nx),STAT=istat)
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)

    ! allocate space for phi functions
    allocate(this%P_sweep(this%nx, this%nnodes - 1, 2))
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
    allocate(this%W_sweep(this%nx, this%nnodes, this%nnodes - 1))
    if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
    
  end subroutine initialize_tmp

  subroutine P_axpy(this,y,f,j,i)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(in) :: y
    class(pf_encap_t), intent(in) :: f
    integer, intent(in) :: j,i
    complex(pfdp),        pointer :: yvec(:), fvec(:)
    class(pf_zndarray_t), pointer :: f_encap,y_encap     
    y_encap => cast_as_zndarray(y)
    call y_encap%get_array(yvec)  
    f_encap => cast_as_zndarray(f)
    call f_encap%get_array(fvec)  
    
    yvec = yvec + this%P_sweep(:, j, i) * fvec
  end subroutine P_axpy
  subroutine P_daxpy(this,y,f1,f2,j,i)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(in) :: y
    class(pf_encap_t), intent(in) :: f1
    class(pf_encap_t), intent(in) :: f2
    integer, intent(in) :: j,i
    complex(pfdp),        pointer :: yvec(:), f1vec(:),f2vec(:)
    class(pf_zndarray_t), pointer :: y_encap,f1_encap, f2_encap
    y_encap => cast_as_zndarray(y)
    call y_encap%get_array(yvec)  
    f1_encap => cast_as_zndarray(f1)
    call f1_encap%get_array(f1vec)  
    f2_encap => cast_as_zndarray(f2)
    call f2_encap%get_array(f2vec)  
    
    yvec = yvec + this%P_sweep(:, j, i) * (f1vec-f2vec)
  end subroutine P_daxpy

  subroutine W_axpy(this,y,f,i,j)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(in) :: y
    class(pf_encap_t), intent(in) :: f
    integer, intent(in) :: j,i

    complex(pfdp),        pointer :: yvec(:), fvec(:)
    class(pf_zndarray_t), pointer :: y_encap,f_encap
    y_encap => cast_as_zndarray(y)
    call y_encap%get_array(yvec)  
    f_encap => cast_as_zndarray(f)
    call f_encap%get_array(fvec)  
    yvec = yvec + this%W_sweep(:, i, j) * fvec
    
  end subroutine W_axpy
  
  ! The rest of the stuff is dimension independent
  include 'sweeper_include.f90'



  ! =======================================================================
  ! INITW   Initializes ETDSDC W functions for vector opL using weights
  !         function by Fornberg.
  !
  ! Arguments
  !
  !   dt   (input) DOUBLE
  !       timestep
  !
  !  This sets the sweeper quantities
  !   W_sweep COMPLEX*16 array, dimension(nnodes-1, nodes, size(opL))
  !       contains the ETD Coefficients so that W(:,:,j) is the ETD integration
  !       matrix cooresponding to lambda=L(j)
  !
  !   P_sweep  COMPLEX*16 array, dimension(2,nnodes-1,size(opL))
  !       contains phi_0 and phi_1 functions needed for ETD Euler method, where
  !       P_sweep(:,:,1,k) = Exp(dt_k L) and P_sweep(:,:,2,k) = dt*dt_i*\phi_{1}(h_k L)
  ! =======================================================================
  subroutine initWSweep(this,dt)

    ! Arguments
    class(my_sweeper_t),  intent(inout) :: this    
    real(pfdp),    intent(in)  :: dt

    ! Local Variables
    real(pfdp) :: eta(this%nnodes-1)
    real(pfdp) :: q(this%nnodes)
    real(pfdp) :: A(this%nnodes,this%nnodes)
    complex(pfdp), allocatable :: p(:,:)
    integer :: i,j,m,nnodes

    ! Initialize variables
    nnodes     = this%nnodes
    eta   = this%nodes(2:nnodes) - this%nodes(1:nnodes-1)
    allocate(p(nnodes+1,this%nx))

    ! Calculate W functions
    do m=1,nnodes-1
      q = (this%nodes - this%nodes(m))/eta(m)                   ! scaled quadrature ppints
      call weights(0.d0,q,nnodes-1,A)                  ! finite difference matrix
      call phi_zvector(eta(m)*dt*this%fft_ops%opL,nnodes,p)            ! phi functions 0 to n
      ! store ith row of Integration matrix
      !      this%W_sweep(:,:,m)   = transpose(dt*eta(m)*matmul(a,p(2:nnodes+1,:)))
      do i=1,nnodes
         do j=1,nnodes
            this%W_sweep(:,i,m)   = this%W_sweep(:,i,m)+ dt*eta(m)*A(i,j)*p(1+j,:)
         end do
      end do
      
      this%P_sweep(:,m,1) = p(1,:)                    ! store exp(dt_i L)
      this%P_sweep(:,m,2) = dt*eta(m)*p(2,:)         ! store dt*eta(i)*phi_1(dt_i L)
    enddo
    deallocate(p)

  end subroutine initWSweep
    
end module pf_my_sweeper


