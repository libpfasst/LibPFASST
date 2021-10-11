!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Some helpful routines which are problem/DIM dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module pf_mod_zutils
  use pfasst
  use pf_mod_solutions
  use pf_mod_zndarray
  use pf_mod_fftpackage
  use fnpy
  implicit none
contains  

  !> Routine to return the exact solution
  subroutine exact(fft,t, y_exact)  
    type(pf_fft_t), pointer, intent(in) :: fft
    real(pfdp), intent(in)  :: t
    type(pf_zndarray_t), intent(inout) :: y_exact
    
    complex(pfdp), pointer :: yex(:)
    complex(pfdp), pointer :: yreal(:)  !  Real space exact solution

    allocate(yreal(y_exact%arr_shape(1)))
    call y_exact%get_array(yex)    
    call exact_realspace(t,yreal)

    call fft%fft(yreal,yex)

    deallocate(yreal)
    
  end subroutine exact
  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, yex)
    use probin, only: eq_type,lam1,lam2,nu, a, t00, kfreq,dom_size,beta,ic_type
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(out) :: yex(:)

    yex=cmplx(0.0,0.0,pfdp)

    select case (eq_type)
    case (0) 
       yex=exp(zi*(lam1+lam2)*t)
    case (1)  !  Advection diffusion
       if (ic_type .eq. 1) then
          call exact_ad_cos(t,yex,nu,a,kfreq,dom_size(1))
       else
          call exact_ad_exp(t,yex,nu,a,dom_size(1))
       end if
    case (2)  !  Burgers
       call exact_burg_sin(t,yex,nu,dom_size(1))
    case (3)  !  NLS
       select case (ic_type)
       case (1)
          call exact_nls(t,yex,dom_size(1))
       case (2)
          call exact_nls_sg(t,yex,dom_size(1))
       case (3)
          call ic_nls_paper(t,yex,dom_size(1))
       end select
    case (4)  !  KdV
       call exact_kdv(t,yex,beta,dom_size(1))
    case (5)  !  KS
       call exact_ks_1d(t,yex,dom_size(1))
    case(6)
       call exact_rd(t,yex,dom_size(1))
    case (7)
       call exact_nls_pert(t,yex,dom_size(1))
       
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
    
  end subroutine exact_realspace

  !  Initial condition for K-S in 1d
  subroutine exact_ks_1d(t,uex,Lx)
    real(pfdp), intent(in) :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
!       uex(i) = 0.25+cos(x/32.0_pfdp)*(1.0_pfdp+sin(x/16.0_pfdp))
       uex(i) = cos(x/16.0_pfdp)*(1.0_pfdp+sin(x/16.0_pfdp))
    end do

  end subroutine exact_ks_1d
  !  Initial condition for a-d-r in 1d  
  subroutine exact_rd(t,uex,Lx)
    real(pfdp), intent(in) :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x,sigma
    sigma=0.2_pfdp
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) -Lx/2.0_pfdp
       uex(i) = 0.45_pfdp + 0.55_pfdp*exp(-x*x*x*x/sigma)
    end do

  end subroutine exact_rd
  subroutine ic_nls_paper(t,uex,Lx)
    real(pfdp), intent(in) :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x,eps

    nx = SIZE(uex)
    eps = 0.01_pfdp
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) -Lx/2.0_pfdp
       uex(i) = 1.0_pfdp + eps*exp(zi*x*0.25_pfdp)
    end do

  end subroutine ic_nls_paper

  !> Routine to return set the linear and nonlinear operators
  subroutine set_ops(opL,opNL,ddx,lap)
    use probin, only: eq_type,lam1,lam2,nu, a,beta,splitting
    complex(pfdp), intent(inout) :: opL(:)
    complex(pfdp), intent(inout) :: opNL(:)
    complex(pfdp), intent(in) :: ddx(:),lap(:)
    real(pfdp) :: cst
    select case (eq_type)
    case (0)  !  Dahlquist
       select case (splitting)
       case (1)  !  Exponential on lam1
          opL = zi*lam1          
       case (2)  !  Exponential on both
          opL = zi*(lam1+lam2)          
       case (3)  !  Exponential on lam2
          opL = zi*lam2          
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
    case (2)  !  Burgers (exponential on advection)
       opL = nu * lap
    case (3)  !  NLS  
       cst = 0.0_pfdp  !  Hyperdiffusion coeff
       opL = zi*lap -cst*lap*lap  !  exponential on i*Lap  and hyperdiffusion
    case (4)  !  KdV  
       opL = -ddx*ddx*ddx/(4.0_pfdp*beta*beta)  ! exponential on third derivative
    case (5)  !  KS  
       opL =-lap*lap-lap  ! exponential on third Lap and Lap^2
    case (6)  !  advection-diffusion-reaction
       opL =nu*lap-a*ddx  ! exponential on third derivative
    case (7)  !  Zero dispersion
       opL = ddx*ddx*ddx
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select

    !  In some cases the "nonlinear term is also done in spectral space"
    select case (eq_type)
    case (0)  !  Dahlquist
       select case (splitting)
       case (1)  !  Exponential on lam1
          opNL= zi*lam2
       case (2)  !  Exponential on both
          opNL= 0.0_pfdp
       case (3)  !  Exponential on lam2
          opNL= zi*lam1
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
    case (2,4,5,6)  !  Burgers, KS
       opNL=ddx
    case (3,7)  !  NLS, zero dispersion
       opNL=0.0_pfdp
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select

  end subroutine set_ops

  !> Routine to compute the nonlinear operators
  subroutine f_NL(yvec,fvec,opNL,tmp,fft)
    use probin, only: eq_type,gamma,beta
    complex(pfdp), intent(in) :: yvec(:)
    complex(pfdp), intent(inout) :: fvec(:)
    complex(pfdp), intent(in) :: opNL(:)
    complex(pfdp), intent(inout) :: tmp(:)
    type(pf_fft_t), intent(in),    pointer :: fft

    logical :: do_dealias=.FALSE.

    !  Set this flag to do dealias
    !do_dealias=.TRUE.

    !  Cases where nonlinear term is actually linear
    if (eq_type .eq. 0 .or. eq_type .eq. 1) then
       fvec=opNL*yvec
       return
    end if

    !  More general cases    
    fvec=yvec
    tmp=yvec
    select case (eq_type)
    case (2,5)  !  Burgers and KS
       if (do_dealias) call fft%dealias(tmp,2)
       call fft%ifft(tmp,tmp)
       tmp=-0.5_pfdp*tmp*tmp
       call fft%fft(tmp,fvec)
       fvec=fvec*opNL
    case (3,7)  !  NLS and zero dispersion
       if (do_dealias) call fft%dealias(tmp,3)
       call fft%ifft(tmp,tmp)
       fvec=conjg(tmp)*tmp*tmp
       call fft%fft(fvec,fvec)
       fvec=2.0_pfdp*zi*fvec
    case (4)  !  KdV
       if (do_dealias)  call fft%dealias(tmp,2)
       call fft%ifft(tmp,tmp)
       tmp=-1.5_pfdp*tmp*tmp
       call fft%fft(tmp,fvec)
       fvec=fvec*opNL
    case (6)  !  advection-reaction diffusion
       if (do_dealias) call fft%dealias(tmp,3)
       call fft%ifft(tmp,tmp)
       fvec=beta*tmp*(1.0_pfdp-tmp)*(0.0_pfdp-tmp)
       call fft%fft(fvec,fvec)
       tmp=gamma*tmp*tmp
       call fft%fft(tmp,tmp)
       fvec=fvec+tmp*opNL
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select

  end subroutine f_NL
  
end module pf_mod_zutils

!  Module to hold the spectral operators
module pf_mod_fftops
  use pf_mod_dtype
  use pf_mod_fftpackage  
  use pf_mod_zutils
  implicit none

  type, public  :: pf_fft_ops_t
     complex(pfdp), allocatable :: lap(:) ! Laplacian operators
     complex(pfdp), allocatable :: ddx(:) ! first derivative operator
     complex(pfdp), allocatable :: opL(:) ! implcit operator
     complex(pfdp), allocatable :: opNL(:) ! explicit operator
     complex(pfdp), allocatable :: opDamp(:) ! Damping operator
   contains
        procedure :: init  =>  fftops_init
        procedure :: destroy  =>  fftops_destroy
  end type pf_fft_ops_t

  contains

    subroutine fftops_init(this,fft,nx)
      use probin, only: d0,d1,r0,r1,split_damp,split_rho
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

      ! Create damping op
      if (split_damp) then
         allocate(this%opDamp(nx),STAT=istat)
         if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
         !      this%opDamp= d0*this%ddx**r0 + d1*abs(this%ddx)**r1
         this%opDamp= abs(this%ddx)**2/tan(two_pi/4.0_pfdp + split_rho)
      
         !  add damping to linear operator
         this%opL=this%opL+this%opDamp
      end if
      deallocate(this%lap)
      deallocate(this%ddx)
    end subroutine fftops_init

    subroutine fftops_destroy(this)
      use probin, only: split_damp
      class(pf_fft_ops_t), intent(inout)    :: this

      deallocate(this%opL)
      deallocate(this%opNL)
      if (split_damp) then
         deallocate(this%opDamp)
      end if
      
    end subroutine fftops_destroy
    
  !> Routine to return out put the solution to numpy (dimension dependent)
  subroutine numpy_dump(fft,t, y_out,fname, do_complex_in)
    type(pf_fft_t), pointer, intent(in) :: fft
    real(pfdp), intent(in)  :: t
    type(pf_zndarray_t), intent(inout) :: y_out
    character(len=*),  intent(in   ) :: fname
    logical,           intent(in   ), optional :: do_complex_in

    complex(pfdp), pointer :: y(:)
    complex(pfdp), pointer :: yreal(:)  !  Real space exact solution
    logical :: do_complex
    
    do_complex=.false.
    if (present(do_complex_in)) do_complex=do_complex_in

    call y_out%get_array(y)  !  Grab the solution from encapsulationi
    
    if (do_complex) then
       call save_complex_double( fname, shape(y), y)
    else  ! output real solution
       allocate(yreal(y_out%arr_shape(1)))
       call fft%ifft(y,yreal)  !  compute the solution in real space
       call save_complex_double( fname, shape(yreal), yreal)
       deallocate(yreal)
    end if
    
  end subroutine numpy_dump
  
    
end module pf_mod_fftops
