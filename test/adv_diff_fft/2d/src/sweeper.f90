!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 2-D advection/diffusion example.
!>     u_t + a*u_x+b*u_y = nu*(u_xx+u_yy)
module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imex_sweeper
  use pf_mod_fftpackage
  use pf_mod_solutions
  implicit none


  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: ad_sweeper_t
     integer ::     nx,ny
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: opE(:,:) ! Explicit part operator
     complex(pfdp), allocatable :: opI(:,:) ! Implicit part operator
     
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves 
     procedure :: destroy     !  Will be called from base sweeper destroy
     procedure :: initialize  !  Will be called from base sweeper initialize

  end type ad_sweeper_t

contains

  !>  Helper function to return sweeper pointer
  function as_ad_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(ad_sweeper_t), pointer :: r
    select type(sweeper)
    type is (ad_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_ad_sweeper

  !>  Routine to initialize variables (called from base sweeper initialize)
  subroutine initialize(this, pf,level_index)
    use probin, only:  imex_stat,nu,a,b
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index


    integer     :: i,nx,ny
    complex(pfdp), allocatable :: lap(:,:)         ! Lapclacian operators
    complex(pfdp), allocatable :: ddx(:,:) ! First derivative operators
    complex(pfdp), allocatable :: ddy(:,:) ! First derivative operators
    
    nx=pf%levels(level_index)%lev_shape(1)
    ny=pf%levels(level_index)%lev_shape(2)
    
    !  Call the imex sweeper initialize
    call this%imex_initialize(pf,level_index)
    !>  Set variables for explicit and implicit parts
    if (imex_stat .eq. 0 ) then
       this%explicit=.TRUE.
       this%implicit=.FALSE.
    elseif (imex_stat .eq. 1 ) then
       this%implicit=.TRUE.
       this%explicit=.FALSE.
    else
       this%implicit=.TRUE.
       this%explicit=.TRUE.
    end if
    

    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nx,ny],2)

    allocate(lap(nx,ny))
    allocate(ddx(nx,ny))
    allocate(ddy(nx,ny))
    allocate(this%opE(nx,ny))
    allocate(this%opI(nx,ny))

    call this%fft_tool%make_lap(lap)
    call this%fft_tool%make_deriv(ddx,1)
    call this%fft_tool%make_deriv(ddy,2)
    select case (imex_stat)
    case (0)  ! Fully Explicit        
       this%opE = -a*ddx-b*ddy + nu * lap
       this%opI = 0.0_pfdp          
    case (1)  ! Fully Implicit
       this%opE = 0.0_pfdp
       this%opI = -a*ddx-b*ddy + nu * lap 
    case (2)  ! IMEX
       this%opE = -a*ddx-b*ddy
       this%opI =  nu * lap           
    case DEFAULT
       print *,'Bad case for imex_stat in f_eval ', imex_stat
       call exit(0)
    end select

    deallocate(lap)
    deallocate(ddx)
    deallocate(ddy)    

  end subroutine initialize

  !>  Destroy sweeper (called from base sweeper destroy)
  subroutine destroy(this,pf,level_index)
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),  target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    
    call this%imex_destroy(pf,level_index)

    deallocate(this%opE)
    deallocate(this%opI)
    
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)

  end subroutine destroy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    use probin, only:  imex_stat ,nu, a,b
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:,:), fvec(:,:)
    type(pf_fft_t),     pointer :: fft

    yvec  => get_array2d(y)
    fvec => get_array2d(f)
    fft => this%fft_tool
    

    ! Apply spectral operators using the FFT convolution function
    select case (piece)
    case (1)  ! Explicit piece
       call fft%conv(yvec,this%opE,fvec)            
    case (2)  ! Implicit piece
       call fft%conv(yvec,this%opI,fvec)            
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select


  end subroutine f_eval

  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  imex_stat ,nu,a,b
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),      pointer :: yvec(:,:), rhsvec(:,:), fvec(:,:)
    type(pf_fft_t),     pointer :: fft

    fft => this%fft_tool

    if (piece == 2) then
       yvec  => get_array2d(y)
       rhsvec => get_array2d(rhs)
       fvec => get_array2d(f)

       ! Apply the inverse operator with the FFT convolution
       call fft%conv(rhsvec,1.0_pfdp/(1.0_pfdp - dtq*this%opI),yvec)

       !  The function is easy to derive
       fvec = (yvec - rhsvec) / dtq
    else
       print *,'Bad piece in f_comp ',piece
       call exit(0)
    end if
  end subroutine f_comp



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial(y_0)
    type(ndarray), intent(inout) :: y_0
    real(pfdp), pointer :: yvec(:,:)
    yvec => get_array2d(y_0)    
    call exact(0.0_pfdp, yvec)
  end subroutine initial

  !> Routine to return the exact solution
  subroutine exact(t, yex)
    use probin, only: nprob,nu, a,b, kfreq, t00
    real(pfdp), intent(in)  :: t
    
    real(pfdp), intent(out) :: yex(:,:)

    integer    :: nx,ny, i, j
    real(pfdp) :: tol, x,y, r2, omega

    nx = size(yex,1)
    ny = size(yex,2)

    !  Using sin wave initial condition
    if (nprob .eq. 1) then
!!$       omega = two_pi*kfreq
!!$       do j = 1, ny
!!$          y = dble(j-1)/dble(ny)-0.5_pfdp - t*b 
!!$          do i = 1, nx
!!$             x = dble(i-1)/dble(nx)-0.5_pfdp - t*a 
!!$             yex(i,j) = sin(omega*x)*sin(omega*y)*exp(-2.0_pfdp*omega*omega*nu*t)
!!$          end do
!!$       end do
       call exact_ad_cos(t,yex,nu,[a,b],[kfreq,kfreq],[1.0_pfdp,1.0_pfdp])
    else
       do j = 1, ny
          y = dble(j-1)/dble(ny)-0.5_pfdp - t*b 
          do i = 1, nx
             x = dble(i-1)/dble(nx)-0.5_pfdp - t*a
             r2=x*x+y*y
             yex(i,j) = t00/(t00+t)*exp(-r2/(4.0*nu*(t00+t)))             
          end do
       end do
    endif
  end subroutine exact
end module pf_my_sweeper

