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
  implicit none

  real(pfdp), parameter ::  Lx     = 1.0_pfdp    ! domain size

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
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
     real(pfdp), allocatable :: work(:)           ! work space
     complex(pfdp), allocatable :: workhat(:)     ! work space
     integer ::  lenwrk, lensav, ierror, nx
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
    integer     :: i,ierror,nx
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
    this%lenwrk = 2*nx 
    this%lensav = 2*nx + int(log(real(nx,kind=8))/log(2.0d+00))+4

    ! create complex fft plans
    allocate(this%workhat(nx))   !  complex transform
    allocate(this%work(this%lenwrk))
    allocate(this%wsave(this%lensav))
    call cfft1i ( nx, this%wsave, this%lensav, ierror )
    if (ierror .ne. 0) then
       stop "error  initializing fft"
    end if

    ! create spectral operators
    allocate(this%ddx(nx))
    allocate(this%lap(nx))
    do i = 1, nx
       if (i <= nx/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-nx + i - 1)
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

    deallocate(this%work)
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
    integer ::  ierror

    yvec  => get_array1d(y)
    fvec => get_array1d(f)

    this%workhat = yvec

    call cfft1f (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in f_eval"
    end if

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

    call cfft1b (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in f_eval"
    end if

    fvec = real(this%workhat)

  end subroutine f_eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
    integer ::  ierror
    
    if (piece == 2) then
       yvec  => get_array1d(y)
       rhsvec => get_array1d(rhs)
       fvec => get_array1d(f)
       this%workhat = rhsvec
       
       call cfft1f (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
       if (ierror .ne. 0) &
          stop "error  calling cfft1f in f_comp"
       if (imex_stat .eq. 2) then
          this%workhat =  this%workhat/ (1.0_pfdp - nu*dtq*this%lap)
       else  ! fully implicit
          this%workhat =  this%workhat/ (1.0_pfdp - dtq*(-v * this%ddx +nu*this%lap))
       end if
       call cfft1b (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )      
       if (ierror .ne. 0) &
                      stop "error  calling cfft1f in f_comp"
      yvec  = real(this%workhat)
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
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),         pointer :: f(:), g(:)

    integer ::  ierror
    adG => as_ad_sweeper(levelG%ulevel%sweeper)
    adF => as_ad_sweeper(levelf%ulevel%sweeper)

    f => get_array1d(qF); 
    g => get_array1d(qG)

    nvarF = size(f)
    nvarG = size(g)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       f = g
       return
    endif

    adG%workhat=g
    call cfft1f (adG%nx, 1, adG%workhat, adG%nx, adG%wsave, adG%lensav, adG%work, adG%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in interpolate"
    end if

    adF%workhat = 0.0d0
    adF%workhat(1:nvarG/2) = adG%workhat(1:nvarG/2)
    adF%workhat(nvarF-nvarG/2+2:nvarF) = adG%workhat(nvarG/2+2:nvarG)


    call cfft1b (adF%nx, 1, adF%workhat, adF%nx, adF%wsave, adF%lensav, adF%work, adF%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in interpolate"
    end if

    f = real(adF%workhat)

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


    real(pfdp), pointer :: f(:), c(:)  

    integer :: irat

    f => get_array1d(qF)
    c => get_array1d(qG)

    irat  = size(f)/size(c)

    c = f(::irat)
  end subroutine restrict


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0
    call exact(0.0_pfdp, q0%flatarray)
  end subroutine initial

  !> Routine to set initial condition.
  subroutine exact(t, yex)
    use probin, only: nprob,nu, v, t00, kfreq
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    integer    :: nx, i, ii, k,nbox
    real(pfdp) :: tol, x, t0,Dx, omega

    nx = size(yex)
    Dx = Lx/dble(nx)

    if (nprob .eq. 1) then
       !  Using sin wave initial condition
       omega = 2*pi*kfreq
       do i = 1, nx
          x = Lx*dble(i-1-nx/2)/dble(nx) - t*v 
          yex(i) = dsin(omega*x)*dexp(-omega*omega*nu*t)
       end do
    else  !  Use periodic image of Gaussians
       yex=0
       if (nu .gt. 0) then
          nbox = ceiling(sqrt(4.0*nu*(t00+t)*37.0d0/(Lx*Lx)))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = (i-1)*Dx-0.5d0 - t*v + dble(k)*Lx
                yex(i) = yex(i) + dsqrt(t00)/dsqrt(t00+t)*dexp(-x*x/(4.0*nu*(t00+t)))
             end do
          end do
       else
          nbox = ceiling(sqrt(37.0d0/(Lx*Lx)))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = i*Dx-0.5d0 - t*v + dble(k)*Lx
                yex(i) = yex(i) + dexp(-x*x)
             end do
          end do
       end if  ! nbox

    end if
  end subroutine exact

end module feval
