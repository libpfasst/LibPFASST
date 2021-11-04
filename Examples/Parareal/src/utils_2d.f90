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
    
    complex(pfdp), pointer :: yex(:,:)
    complex(pfdp), pointer :: yreal(:,:)  !  Real space exact solution

    allocate(yreal(y_exact%arr_shape(1),y_exact%arr_shape(2)))

    call y_exact%get_array(yex)    
    call exact_realspace(t,yreal)

    call fft%fft(yreal,yex)

    deallocate(yreal)
    
  end subroutine exact
  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, yex)
    use probin, only: eq_type,lam1,lam2,nu, a,b, t00, kfreq,dom_size,beta,ic_type
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(out) :: yex(:,:)

    yex=cmplx(0.0,0.0,pfdp)

    select case (eq_type)
    case (0)
       yex=exp((lam1+lam2)*t)
    case (1)  !  Advection diffusion
       if (ic_type .eq. 1) then
          call exact_ad_cos(t,yex,nu,[a,b],kfreq,dom_size)
       else
          call exact_ad_exp(t,yex,nu,[a,b],dom_size)
       end if
    case (2)  !  Burgers
       call exact_burg_sin(t,yex,nu,dom_size)
    case (3)  !  NLS
       select case (ic_type)
       case (1)
          call exact_nls(t,yex,dom_size)
       case (2)
          call exact_nls_sg(t,yex,dom_size)
       case (3)
          call ic_nls_paper(t,yex,dom_size)
       end select
    case (4)  !  KdV
       !       call exact_kdv(t,yex,beta,dom_size)
    case (7)  !  NLS
!       call init_nls_twowavez(t,yex,dom_size)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case  for eq_type ',eq_type)
    end select
    
  end subroutine exact_realspace
  subroutine ic_nls_paper(t, uex,Lx)
    use probin, only: beta
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,eps
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    eps=0.02_pfdp
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = 1.0_pfdp + eps*exp(zi*(x+y)*0.25_pfdp) + eps/10.0*exp(zi*(x-y)*0.5_pfdp)
       end do
    end do
       
  end subroutine ic_nls_paper
  
  !> Routine to return set the linear and nonlinear operators
  subroutine set_ops(opL,opNL,ddx,ddy,lap)
    use probin, only: eq_type,lam1,lam2,nu, a,b,beta,splitting
    complex(pfdp), intent(inout) :: opL(:,:)
    complex(pfdp), intent(inout) :: opNL(:,:)
    complex(pfdp), intent(in) :: ddx(:,:),ddy(:,:),lap(:,:)
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
          opL = nu*lap -a*ddx-b*ddy
       case (3)  !  Exponential on advection
          opL = -a*ddx-b*ddy                    
       end select
    case (2)  !  Burgers
       opL = nu * lap
    case (3)  !  NLS
       opL = zi*lap
    case (4)  !  KdV
       opL = -ddx*ddx*ddx/(4.0_pfdp*beta*beta)
    case (5)  !  KS
       opL =-lap*lap
    case (7)  !  Zero dispersion
       opL=ddx*ddx*ddx+ddy*ddy*ddy
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select

    !  In some cases the "nonlinear term is also done in spectral space"
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
          opNL= -a*ddx-b*ddy                    
       case (2)  !  Exponential on both
          opNL= 0.0_pfdp
       case (3)  !  Exponential on advection
          opNL= nu*lap
       end select
    case (2,5)  !  Burgers, KS
       opNL=ddx+ddy
    case (3)  !  NLS
       opNL=0.0_pfdp
    case (4)  !  kdv
       opNL=ddx
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
    
  end subroutine set_ops

  !> Routine to compute the nonlinear operators
  subroutine f_NL(yvec,fvec,opNL,tmp,fft)
    use probin, only: eq_type,beta,gamma
    complex(pfdp), intent(in) :: yvec(:,:)
    complex(pfdp), intent(inout) :: fvec(:,:)
    complex(pfdp), intent(in) :: opNL(:,:)
    complex(pfdp), intent(inout) :: tmp(:,:)
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
       call fft%ifft(yvec,tmp)
       tmp=-0.5_pfdp*tmp*tmp
       call fft%fft(tmp,fvec)
       fvec=fvec*opNL
    case (3,7)  !  NLS
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
     complex(pfdp), allocatable :: lap(:,:) ! Laplacian operators
     complex(pfdp), allocatable :: ddx(:,:) ! first derivative operator
     complex(pfdp), allocatable :: ddy(:,:) ! first derivative operator
     complex(pfdp), allocatable :: opL(:,:) ! implcit operator
     complex(pfdp), allocatable :: opNL(:,:) ! explicit operator
     complex(pfdp), allocatable :: opDamp(:,:) ! explicit operator
   contains
        procedure :: init  =>  fftops_init
        procedure :: destroy  =>  fftops_destroy
  end type pf_fft_ops_t

  contains

    subroutine fftops_init(this,fft,nx)
      use probin, only: split_damp,split_rho
      class(pf_fft_ops_t), intent(inout)    :: this
      type(pf_fft_t), pointer, intent(in) :: fft
      integer, intent(in) :: nx

      integer :: istat
      allocate(this%lap(nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddx(nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddy(nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opL(nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opNL(nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      
      call fft%make_deriv(this%ddx,1) !  First derivative
      call fft%make_deriv(this%ddy,2) !  First derivative
      call fft%make_lap(this%lap)  !  Second derivative
      
      ! initialize  operators
      call set_ops(this%opL,this%opNL,this%ddx,this%ddy,this%lap)

      ! Create damping op
      if (split_damp) then
         allocate(this%opDamp(nx,nx),STAT=istat)
         if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
         !      this%opDamp= d0*this%ddx**r0 + d1*abs(this%ddx)**r1
         this%opDamp= abs(this%lap)/tan(two_pi/4.0_pfdp + split_rho)
      
         !  add damping to linear operator
         this%opL=this%opL+this%opDamp
      end if

      
      deallocate(this%lap)
      deallocate(this%ddx)
      deallocate(this%ddy)
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

    complex(pfdp), pointer :: y(:,:)
    complex(pfdp), pointer :: yreal(:,:)  !  Real space exact solution
    logical :: do_complex
    
    do_complex=.false.
    if (present(do_complex_in)) do_complex=do_complex_in

    call y_out%get_array(y)  !  Grab the solution from encapsulationi
    
    if (do_complex) then
       call save_complex_double( fname, shape(y), y)
    else  ! output real solution
       allocate(yreal(y_out%arr_shape(1),y_out%arr_shape(2)))
       call fft%ifft(y,yreal)  !  compute the solution in real space
       call save_complex_double( fname, shape(yreal), yreal)
       deallocate(yreal)
    end if
    
  end subroutine numpy_dump
  
    
end module pf_mod_fftops
