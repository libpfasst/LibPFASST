!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 1-D advection/diffusion example.
!>     u_t + v*u_x = nu*u_xx
module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imex_sweeper
  use pf_mod_fftpackage
  use pf_mod_solutions
  implicit none

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: ad_sweeper_t
     integer ::     nx   !  Grid size

     !>  FFT and Spectral derivatives
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: opE(:) ! Explicit spectral operator
     complex(pfdp), allocatable :: opI(:) ! Implicit spectral operator
     
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


  !>  Routine to initialize sweeper (bypasses imex sweeper initialize)
  subroutine initialize(this, pf,level_index)
    use probin, only:  imex_stat, v ,nu
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index

    complex(pfdp), allocatable :: lap(:) ! Lapclacian operator
    complex(pfdp), allocatable :: ddx(:) ! First derivative operator

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

    nx=pf%levels(level_index)%lev_shape(1)  !  local convenient grid size

    !>  Set up the FFT 
    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nx],1)

    !>  Define spectral derivatitive operators
    allocate(lap(nx))
    allocate(ddx(nx))
    call this%fft_tool%make_lap(lap)
    call this%fft_tool%make_deriv(ddx)
    
    !> Allocate operators for implicit and explicit parts
    allocate(this%opE(nx))
    allocate(this%opI(nx))
    
    !>  Choose the explicit and implicit operators depending on imex_stat
    select case (imex_stat)
    case (0)  ! Fully Explicit        
       this%opE = -v*ddx + nu*lap
       this%opI = 0.0_pfdp          
    case (1)  ! Fully Implicit
       this%opE = 0.0_pfdp
       this%opI = -v*ddx + nu * lap 
    case (2)  ! IMEX
       this%opE = -v*ddx
       this%opI =  nu*lap           
    case DEFAULT
       print *,'Bad case for imex_stat in f_eval ', imex_stat
       call exit(0)
    end select

    !>  Clean up locals
    deallocate(lap)
    deallocate(ddx)
    
  end subroutine initialize
  

  !>  Destroy sweeper (bypasses base sweeper destroy)
  subroutine destroy(this,pf,level_index)
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),  target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    !>  Call base sweeper destroy
    call this%imex_destroy(pf,level_index)

    !> Nuke the FFT operators 
    deallocate(this%opE)
    deallocate(this%opI)

    !>  Free up FFT stuff
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
    integer,             intent(in   ) :: piece  !  Which piece to solve for
    
    real(pfdp),      pointer :: yvec(:), fvec(:)
    type(pf_fft_t),     pointer :: fft

    !  Grab the arrays from the encap
    yvec  => get_array1d(y)
    fvec => get_array1d(f)
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

  ! Solve for y and return f2 also
  !   y-dtq*f(y,t) = rhs
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  imex_stat ,nu,v
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y       !  The solution we seek
    real(pfdp),          intent(in   ) :: t       !  Equation time of implicit solve
    real(pfdp),          intent(in   ) :: dtq     !  The 
    class(pf_encap_t),   intent(in   ) :: rhs     !  The right hand side of the solve
    integer,             intent(in   ) :: level_index !  Which level this is
    class(pf_encap_t),   intent(inout) :: f       !  The function value
    integer,             intent(in   ) :: piece   !  Designates which piece to solve for (here implicit)

    real(pfdp),      pointer :: yvec(:), rhsvec(:), fvec(:)
    type(pf_fft_t),     pointer :: fft

    !  Grab the arrays from the encaps
    yvec  => get_array1d(y)
    rhsvec => get_array1d(rhs)
    fvec => get_array1d(f)

    if (imex_stat .eq. 0)  then
       print *,'We should not be calling fcomp for fully explicit'
       yvec=rhsvec
       fvec=0.0_pfdp
       return
    endif

    ! Grab the fft workspace
    fft => this%fft_tool

    if (piece == 2) then
       ! Apply the inverse opeator with the FFT convolution
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
  subroutine initial_sol(y_0)
    type(pf_ndarray_t), intent(inout) :: y_0
    call exact(0.0_pfdp, y_0%flatarray)
  end subroutine initial_sol

  !> Routine to return the exact solution
  subroutine exact(t, yex)
    use probin, only: nu, v,kfreq,Lx
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)


    !  Call exact solution from Libpfasst for ad problem
    call exact_ad_cos_1d(t,yex,nu,v,kfreq,Lx)
  end subroutine exact


end module pf_my_sweeper
