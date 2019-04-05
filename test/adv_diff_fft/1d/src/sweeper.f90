!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 1-D advection/diffusion example.
!>     u_t + v*u_x = nu*u_xx
module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imex_sweeper
  use pf_mod_fftpackage
  implicit none

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: ad_sweeper_t
     integer ::     nx   !  Grid size
     !>  Spectral derivative operators
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: lap(:) ! Lapclacian operators
     complex(pfdp), allocatable :: ddx(:) ! First derivative operators
     complex(pfdp), allocatable :: opE(:)         ! Lapclacian operators
     complex(pfdp), allocatable :: opI(:)         ! First derivative operators
     
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves
     procedure :: initialize  !  Bypasses base sweeper initialize
     procedure :: destroy     !  Bypasses base sweeper destroy

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


  !>  Routine to initialize variables (bypasses imex sweeper initialize)
  subroutine initialize(this, pf,level_index)
    use probin, only:  imex_stat, v ,nu
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index


    integer     :: nx

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

    nx=pf%levels(level_index)%shape(1)  !  local convenient grid size

    !>  Set up the FFT based operators
    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nx],1)

    allocate(this%lap(nx))
    allocate(this%ddx(nx))
    allocate(this%opE(nx))
    allocate(this%opI(nx))
    call this%fft_tool%make_lap_1d(this%lap)
    call this%fft_tool%make_deriv_1d(this%ddx)    

    select case (imex_stat)
    case (0)  ! Fully Explicit        
       this%opE = -v*this%ddx + nu * this%lap
       this%opI = 0.0_pfdp          
    case (1)  ! Fully Implicit
       this%opE = 0.0_pfdp
       this%opI = -v*this%ddx + nu * this%lap 
    case (2)  ! IMEX
       this%opE = -v*this%ddx
       this%opI =  nu * this%lap           
    case DEFAULT
       print *,'Bad case for imex_stat in f_eval ', imex_stat
       call exit(0)
    end select
    deallocate(this%lap)
    deallocate(this%ddx)
    
  end subroutine initialize
  

  !>  Destroy sweeper (bypasses base sweeper destroy)
  subroutine destroy(this,pf,level_index)
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),  target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    !>  Call base sweeper destroy
    call this%imex_destroy(pf,level_index)

    !> Nuke the local stuff
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
    use probin, only:  imex_stat ,nu, v
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:), fvec(:)
    type(pf_fft_t),     pointer :: fft
    complex(pfdp),      pointer :: wk(:)

    yvec  => get_array1d(y)
    fvec => get_array1d(f)
    fft => this%fft_tool
    wk => fft%get_wk_ptr_1d()
    
    ! Load the solution into the FFT
    wk=yvec

    ! Apply spectral operators using the FFT convolution function
    select case (piece)
    case (1)  ! Explicit piece
       call fft%conv(this%opE)            
    case (2)  ! Implicit piece
       call fft%conv(this%opI)            
    case DEFAULT
       print *,'Bad case for piece in f_eval ', piece
       call exit(0)
    end select

    fvec=real(wk)
  end subroutine f_eval

  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  imex_stat ,nu,v
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),      pointer :: yvec(:), rhsvec(:), fvec(:)
    complex(pfdp),      pointer :: wk(:)
    type(pf_fft_t),     pointer :: fft

    fft => this%fft_tool
    wk => fft%get_wk_ptr_1d()
    if (piece == 2) then
       yvec  => get_array1d(y)
       rhsvec => get_array1d(rhs)
       fvec => get_array1d(f)

       ! Apply the inverse opeator with the FFT convolution
       wk=rhsvec
       call fft%conv(1.0_pfdp/(1.0_pfdp - dtq*this%opI))

       yvec=real(wk)
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
    call exact(0.0_pfdp, y_0%flatarray)
  end subroutine initial

  !> Routine to return the exact solution
  subroutine exact(t, yex)
    use probin, only: nprob,nu, v, t00, kfreq
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    integer    :: nx, i, ii, k,nbox
    real(pfdp) :: tol, t0,Dx, omega
    real(pfdp), allocatable ::  x

    nx = size(yex)
    Dx = 1.0d0/dble(nx)

    if (nprob .eq. 1) then
       !  Using sin wave initial condition
       omega = two_pi*kfreq
       do i = 1, nx
          x = dble(i-1-nx/2)/dble(nx) - t*v 
          yex(i) = cos(omega*x)*exp(-omega*omega*nu*t)
       end do
    else  !  Use periodic image of Gaussians
       yex=0
       if (nu .gt. 0.0) then
          nbox = ceiling(sqrt(4.0*nu*(t00+t)*37.0d0))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = dble(i-1-nx/2)/dble(nx) - t*v  + dble(k)
                yex(i) = yex(i) + sqrt(t00)/sqrt(t00+t)*exp(-x*x/(4.0*nu*(t00+t)))
             end do
          end do
       else
          nbox = ceiling(sqrt(37.0d0))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = dble(i-1-nx/2)/dble(nx) - t*v  + dble(k)
                yex(i) = yex(i) + exp(-x*x)
             end do
          end do
       end if  ! nbox

    end if

  end subroutine exact


end module feval
