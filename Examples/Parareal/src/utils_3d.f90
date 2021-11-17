module pf_mod_zutils
  use pfasst
  use pf_mod_solutions
  use pf_mod_zndarray
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
    type(pf_zndarray_t), intent(inout) :: y_exact
    
    complex(pfdp), pointer :: yex(:,:,:)
    complex(pfdp), pointer :: yreal(:,:,:)  !  Real space exact solution

    allocate(yreal(y_exact%arr_shape(1),y_exact%arr_shape(2),y_exact%arr_shape(3)))

    call y_exact%get_array(yex)    
    call exact_realspace(t,yreal)

    call fft%fft(yreal,yex)

    deallocate(yreal)
    
  end subroutine exact

  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, yex)
    use probin, only: eq_type,lam1,lam2,nu, a,b,c, t00, kfreq,beta,ic_type,dom_size
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(out) :: yex(:,:,:)

    yex=cmplx(0.0,0.0,pfdp)

    select case (eq_type)
    case (0)
       yex=exp((lam1+lam2)*t)
    case (1)  !  Advection diffusion
       if (ic_type .eq. 1) then 
          call exact_ad_cos(t,yex,nu,[a,b,c],kfreq,dom_size)
       else
          call exact_ad_exp(t,yex,nu,[a,b,c],dom_size)
       end if
    case (2)  !  Burgers
       call exact_burg_sin(t,yex,nu,dom_size)
    case (3)  !  NLS
       call exact_nls(t,yex,dom_size)
    case (4)  !  KdV
!       call exact_kdv(t,yex,beta,dom_size)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case  for eq_type ',eq_type)
    end select
    
  end subroutine exact_realspace

  !> Routine to return set the linear and nonlinear operators
  subroutine set_ops(opL,opNL,ddx,ddy,ddz,lap)
    use probin, only: eq_type,lam1,lam2,nu, a,b,c,beta,splitting
    complex(pfdp), intent(inout) :: opL(:,:,:)
    complex(pfdp), intent(inout) :: opNL(:,:,:)
    complex(pfdp), intent(in) :: ddx(:,:,:),ddy(:,:,:),ddz(:,:,:),lap(:,:,:)
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
          opL = nu*lap -a*ddx-b*ddy-c*ddz                                        
       case (3)  !  Exponential on advection
          opL = -a*ddx-b*ddy-c*ddz                                                            
       end select
    case (2)  !  Burgers
       opL = nu * lap
    case (3)  !  NLS
       opL = zi*lap
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
          opNL= -a*ddx-b*ddy-c*ddz
       case (2)  !  Exponential on both
          opNL= 0.0_pfdp
       case (3)  !  Exponential on advection
          opNL= nu*lap
       end select
    case (2,5)  !  Burgers, KS
       opNL=ddx+ddy+ddz
    case (3)  !  NSLS
       opNL=0.0_pfdp
    case (4)  !  kdv
       opNL=ddx
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
  end subroutine set_ops
  
  !> Routine to return set the linear and nonlinear operators
  subroutine f_NL(yvec,fvec,opNL,tmp,fft)
    use probin, only: eq_type
    complex(pfdp), intent(in) :: yvec(:,:,:)
    complex(pfdp), intent(inout) :: fvec(:,:,:)
    complex(pfdp), intent(in) :: opNL(:,:,:)
    complex(pfdp), intent(inout) :: tmp(:,:,:)
    type(pf_fft_t), intent(in),    pointer :: fft


    if (eq_type .eq. 0 .or. eq_type .eq. 1) then
       fvec=opNL*yvec
       return
    end if

    call fft%ifft(yvec,tmp)
    select case (eq_type)
    case (2,5)  !  Burgers and KS
       tmp=-0.5_pfdp*tmp*tmp
       call fft%fft(tmp,fvec)
       fvec=fvec*opNL
    case (3)  !  NLS
       tmp=6.0_pfdp*zi*abs(tmp)*abs(tmp)*tmp
       call fft%fft(tmp,fvec)
    case (4)  !  KdV
       tmp=-1.5_pfdp*tmp*tmp
       call fft%fft(tmp,fvec)
       fvec=fvec*opNL
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
  end subroutine f_NL


end module pf_mod_zutils
module pf_mod_fftops
  use pf_mod_dtype
  use pf_mod_fftpackage  
  use pf_mod_zutils
  use fnpy
  implicit none


  type, public  :: pf_fft_ops_t
     complex(pfdp), allocatable :: lap(:,:,:) ! Laplacian operators
     complex(pfdp), allocatable :: ddx(:,:,:) ! first derivative operator
     complex(pfdp), allocatable :: ddy(:,:,:) ! first derivative operator
     complex(pfdp), allocatable :: ddz(:,:,:) ! first derivative operator
     complex(pfdp), allocatable :: opL(:,:,:) ! implcit operator
     complex(pfdp), allocatable :: opNL(:,:,:) ! explicit operator
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
      allocate(this%lap(nx,nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddx(nx,nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddy(nx,nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddz(nx,nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opL(nx,nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opNL(nx,nx,nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      
      call fft%make_deriv(this%ddx,1) !  First derivative
      call fft%make_deriv(this%ddy,2) !  First derivative
      call fft%make_deriv(this%ddz,3) !  First derivative
      call fft%make_lap(this%lap)  !  Second derivative
      
      ! initialize  operators
      call set_ops(this%opL,this%opNL,this%ddx,this%ddy,this%ddz,this%lap)
      deallocate(this%lap)
      deallocate(this%ddx)
      deallocate(this%ddy)
      deallocate(this%ddz)
    end subroutine fftops_init

    subroutine fftops_destroy(this)
      class(pf_fft_ops_t), intent(inout)    :: this

      deallocate(this%opL)
      deallocate(this%opNL)
    end subroutine fftops_destroy
    
  !> Routine to return out put the solution to numpy (dimension dependent)
  subroutine numpy_dump(fft,t, y_out,fname, do_complex_in)
    type(pf_fft_t), pointer, intent(in) :: fft
    real(pfdp), intent(in)  :: t
    type(pf_zndarray_t), intent(inout) :: y_out
    character(len=*),  intent(in   ) :: fname
    logical,           intent(in   ), optional :: do_complex_in

    complex(pfdp), pointer :: y(:,:,:)
    complex(pfdp), pointer :: yreal(:,:,:)  !  Real space exact solution
    logical :: do_complex
    
    do_complex=.false.
    if (present(do_complex_in)) do_complex=do_complex_in

    call y_out%get_array(y)  !  Grab the solution from encapsulationi
    
    if (do_complex) then
       call save_complex_double( fname, shape(y), y)
    else  ! output real solution
       allocate(yreal(y_out%arr_shape(1),y_out%arr_shape(2),y_out%arr_shape(3)))
       call fft%ifft(y,yreal)  !  compute the solution in real space
       call save_complex_double( fname, shape(yreal), yreal)
       deallocate(yreal)
    end if
    
  end subroutine numpy_dump

    
end module pf_mod_fftops
