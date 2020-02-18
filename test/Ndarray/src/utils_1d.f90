module pf_mod_rutils
  use pf_mod_solutions
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_fftpackage  
  implicit none
  
contains

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Routine to return the exact solution
  subroutine exact(fft,t, y_exact)  

    type(pf_fft_t), pointer, intent(in) :: fft
    real(pfdp), intent(in)  :: t
    type(pf_ndarray_t), intent(inout) :: y_exact
    
    real(pfdp), pointer :: yex(:)

    call y_exact%get_array(yex)    
    call exact_realspace(t,yex)

  end subroutine exact

  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, yex)
    use probin, only: eq_type,lam1,lam2,nu, a, t00, kfreq,Lx,beta,ic_type
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    yex=real(0.0,pfdp)

    select case (eq_type)
    case (0)
       yex=exp((lam1+lam2)*t)
    case (1)  !  Advection diffusion
       if (ic_type .eq. 1) then
          call exact_ad_cos(t,yex,nu,a,kfreq,Lx)
       else
          call exact_ad_exp(t,yex,nu,a,Lx)
       end if
    case (2)  !  Burgers
       call exact_burg_sin(t,yex,nu,Lx)
    case (4)  !  KdV
       call exact_kdv(t,yex,beta,Lx)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
    
  end subroutine exact_realspace

  !> Routine to return set the linear and nonlinear operators
  subroutine set_ops(opL,opNL,ddx,lap)
    use probin, only: eq_type,lam1,lam2,nu, a,beta,splitting
    complex(pfdp), intent(inout) :: opL(:)
    complex(pfdp), intent(inout) :: opNL(:)
    complex(pfdp), intent(in) :: ddx(:),lap(:)
    select case (eq_type)
    case (0)  !  Dahlquist
       select case (splitting)
       case (1)  !  Exponential on diffusion
          opL = lam1          
       case (2)  !  Exponential on both
          opL = lam1+lam2          
       case (3)  !  Exponential on advection
          opL = lam2          
       end select
    case (1)  !  Advection diffusion
       select case (splitting)
       case (1)  !  Exponential on diffusion
          opL = nu*lap          
       case (2)  !  Exponential on both
          opL = nu*lap -a*ddx          
       case (3)  !  Exponential on advection
          opL = -a*ddx
       end select
    case (2)  !  Burgers
       opL = nu * lap
    case (4)  !  KdV
       opL = -ddx*ddx*ddx/(4.0_pfdp*beta*beta)
    case (5)  !  KS
       opL =-lap*lap
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select

    select case (eq_type)
    case (0)  !  Dahlquist
       select case (splitting)
       case (1)  !  Exponential on diffusion
          opNL= lam2
       case (2)  !  Exponential on both
          opNL= 0.0_pfdp
       case (3)  !  Exponential on advection
          opNL= lam1
       end select
    case (1)  !  Advection diffusion
       select case (splitting)
       case (1)  !  Exponential on diffusion
          opNL= -a*ddx
       case (2)  !  Exponential on both
          opNL= 0.0_pfdp
       case (3)  !  Exponential on advection
          opNL= nu*lap
       end select
    case (2,5)  !  Burgers, KS
       opNL=ddx
    case (4)  !  KDV
       opNL=ddx
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
    
  end subroutine set_ops

  !> Routine to return set the linear and nonlinear operators
  subroutine f_NL(yvec,fvec,opNL,tmp,fft)
    use probin, only: eq_type
    real(pfdp), intent(in) :: yvec(:)
    real(pfdp), intent(inout) :: fvec(:)
    complex(pfdp), intent(in) :: opNL(:)
    real(pfdp), intent(inout) :: tmp(:)
    type(pf_fft_t), intent(in),    pointer :: fft

!    this%fvec=this%yvec
!    if (nx>3) then
!       this%fvec(K:nx-K+1)=0.0
    !    end if
    select case (eq_type)
    case (0)  !  Dahlquist
       fvec=yvec*opNL
    case (1)  !  advection
       call fft%conv(yvec,opNL,fvec)
    case (2,5)  !  Burgers and KS       
       tmp=-0.5_pfdp*yvec*yvec
       call fft%conv(tmp,opNL,fvec)
    case (4)  !  KdV
       tmp=-1.5_pfdp*yvec*yvec
       call fft%conv(tmp,opNL,fvec)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
  end subroutine f_NL
  
  
end module pf_mod_rutils

module pf_mod_fftops
  use pf_mod_dtype
  use pf_mod_fftpackage  
  use pf_mod_rutils
  implicit none


  type, public  :: pf_fft_ops_t
     complex(pfdp), allocatable :: lap(:) ! Laplacian operators
     complex(pfdp), allocatable :: ddx(:) ! first derivative operator
     complex(pfdp), allocatable :: opL(:) ! implcit operator
     complex(pfdp), allocatable :: opNL(:) ! explicit operator
   contains
        procedure :: init  =>  fftops_init
        procedure :: destroy  =>  fftops_destroy
  end type pf_fft_ops_t

  contains

    subroutine fftops_init(this,fft,nx)
      class(pf_fft_ops_t), intent(inout)    :: this
      type(pf_fft_t), pointer, intent(in) :: fft
      integer, intent(in) :: nx
     
      integer :: istat
      allocate(this%lap(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddx(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opL(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opNL(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      
      call fft%make_deriv(this%ddx) !  First derivative
      call fft%make_lap(this%lap)  !  Second derivative
      
      ! initialize  operators
      call set_ops(this%opL,this%opNL,this%ddx,this%lap)
      deallocate(this%lap)
      deallocate(this%ddx)
    end subroutine fftops_init

    subroutine fftops_destroy(this)
      class(pf_fft_ops_t), intent(inout)    :: this

      deallocate(this%opL)
      deallocate(this%opNL)
    end subroutine fftops_destroy
    
    
end module pf_mod_fftops
