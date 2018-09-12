!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 1-D advection/diffusion example.
!>     u_t + v*u_x = nu*u_xx
module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imexQ


  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imexQ_t) :: ad_sweeper_t
     real(pfdp), allocatable :: wsave(:)          ! work space
     complex(pfdp), allocatable :: workhat(:)     ! work space
     integer ::   lensav,  nx
     complex(pfdp), allocatable :: ddx(:), lap(:) ! spectral operators
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves 

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

  !>  Routine to set up sweeper variables and operators
  subroutine sweeper_setup(sweeper, grid_shape)
    use probin, only:  imex_stat 
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: grid_shape(1)

    class(ad_sweeper_t), pointer :: this
    integer     :: i,nx
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    nx=grid_shape(1)
    this => as_ad_sweeper(sweeper)

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

    !  FFT Storage
    this%nx = nx 
    this%lensav = 4*nx + 15

    ! create complex fft plans
    allocate(this%workhat(nx))   !  complex transform
    allocate(this%wsave(this%lensav))
    !  Initialize FFT
    call ZFFTI( nx, this%wsave )

    ! create spectral operators
    allocate(this%ddx(nx))
    allocate(this%lap(nx))
    do i = 1, nx
       if (i <= nx/2+1) then
          kx = two_pi * dble(i-1)
       else
          kx = two_pi * dble(-nx + i - 1)
       end if

       this%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          this%lap(i) = 0.0_pfdp
       else
          this%lap(i) = -kx**2
       end if
    end do
  end subroutine sweeper_setup

  !>  destroy the sweeper type
  subroutine destroy(this, lev)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev

    deallocate(this%workhat)
    deallocate(this%wsave)
    deallocate(this%ddx)
    deallocate(this%lap)
    
    call this%imexQ_destroy(lev)

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

    yvec  => get_array1d(y)
    fvec => get_array1d(f)


    ! Take the fft of y
    this%workhat = yvec
    call zfftf(this%nx, this%workhat, this%wsave )    

    ! Apply spectral operators
    select case (piece)
    case (1)  ! Explicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          this%workhat = (-v * this%ddx + nu * this%lap) *  this%workhat
       case (1)  ! Fully Implicit
          print *,'Should not be in this case in feval'
       case (2)  ! IMEX
          this%workhat = -v * this%ddx *  this%workhat
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case (2)  ! Implicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          print *,'Should not be in this case in feval'
       case (1)  ! Fully Implicit
          this%workhat = (-v * this%ddx + nu * this%lap) *  this%workhat
       case (2)  ! IMEX
          this%workhat = (nu * this%lap) *  this%workhat
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select

    !  Normalize the fft
    this%workhat=this%workhat/dble(this%nx)

    !  Do the inverse fft
    call zfftb(this%nx, this%workhat, this%wsave)

    !  The function is the real part
    fvec = real(this%workhat)

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
    
    if (piece == 2) then
       yvec  => get_array1d(y)
       rhsvec => get_array1d(rhs)
       fvec => get_array1d(f)
       this%workhat = rhsvec

       !  Take the FFT
       call zfftf(this%nx, this%workhat, this%wsave )

       !  Apply spectral inverse operators
       if (imex_stat .eq. 2) then
          this%workhat =  this%workhat/ (1.0_pfdp - nu*dtq*this%lap)
       else  ! fully implicit
          this%workhat =  this%workhat/ (1.0_pfdp - dtq*(-v * this%ddx +nu*this%lap))
       end if

       !  Normalize the FFT
       this%workhat=this%workhat/dble(this%nx)

       !  Take the inverse FFT
       call zfftb(this%nx, this%workhat,this%wsave)

       !  The solution is the real part
       yvec  = real(this%workhat)

       !  The function is easy to derive
       fvec = (yvec - rhsvec) / dtq
    else
       print *,'Bad piece in f_comp ',piece
       call exit(0)
    end if
  end subroutine f_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  These are the transfer functions that must be  provided for the level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF,qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags


    integer :: nvarF, nvarG, xrat
    class(ad_sweeper_t), pointer :: sweeper_f, sweeper_c
    real(pfdp),         pointer :: yvec_f(:), yvec_c(:)

    sweeper_c => as_ad_sweeper(levelG%ulevel%sweeper)
    sweeper_f => as_ad_sweeper(levelf%ulevel%sweeper)

    yvec_f => get_array1d(qF); 
    yvec_c => get_array1d(qG)

    nvarF = size(yvec_f)
    nvarG = size(yvec_c)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       yvec_f = yvec_c
       return
    endif

    sweeper_c%workhat=yvec_c
    call zfftf(sweeper_c%nx, sweeper_c%workhat, sweeper_c%wsave )    

    sweeper_f%workhat = 0.0d0
    sweeper_f%workhat(1:nvarG/2) = sweeper_c%workhat(1:nvarG/2)
    sweeper_f%workhat(nvarF-nvarG/2+2:nvarF) = sweeper_c%workhat(nvarG/2+2:nvarG)

    sweeper_f%workhat=sweeper_f%workhat/dble(sweeper_c%nx)
    call zfftb(sweeper_f%nx, sweeper_f%workhat, sweeper_f%wsave)    

    yvec_f = real(sweeper_f%workhat)

  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, levelf, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelf  !<  fine level
    class(pf_level_t), intent(inout) :: levelG  !<  coarse level
    class(pf_encap_t), intent(inout) :: qF    !<  fine solution
    class(pf_encap_t), intent(inout) :: qG    !<  coarse solution
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags


    real(pfdp), pointer :: yvec_f(:), yvec_c(:)  

    integer :: irat

    yvec_f => get_array1d(qF)
    yvec_c => get_array1d(qG)

    irat  = size(yvec_f)/size(yvec_c)

    yvec_c = yvec_f(::irat)
  end subroutine restrict


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
    real(pfdp) :: tol, x, t0,Dx, omega

    nx = size(yex)
    Dx = 1.0d0/dble(nx)

    if (nprob .eq. 1) then
       !  Using sin wave initial condition
       omega = two_pi*kfreq
       do i = 1, nx
          x = dble(i-1-nx/2)/dble(nx) - t*v 
          yex(i) = sin(omega*x)*exp(-omega*omega*nu*t)
       end do
    else  !  Use periodic image of Gaussians
       yex=0
       if (nu .gt. 0) then
          nbox = ceiling(sqrt(4.0*nu*(t00+t)*37.0d0))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = (i-1)*Dx-0.5d0 - t*v + dble(k)
                yex(i) = yex(i) + sqrt(t00)/sqrt(t00+t)*exp(-x*x/(4.0*nu*(t00+t)))
             end do
          end do
       else
          nbox = ceiling(sqrt(37.0d0))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = i*Dx-0.5d0 - t*v + dble(k)
                yex(i) = yex(i) + exp(-x*x)
             end do
          end do
       end if  ! nbox

    end if
  end subroutine exact


end module feval
