!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for 1-D advection/diffusion example.
!     u_t + v*u_x = nu*u_xx
module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imex
  use pf_mod_imexQ
  implicit none
!  include 'fftw3.f03'

  real(pfdp), parameter :: &
       Lx     = 1.0_pfdp, &    ! domain size
       kfreq  = 8.0_pfdp, &    ! Frequency of initial conditions
       v      = 1.0_pfdp, &    ! advection velocity
       nu     = 0.02_pfdp, &   ! viscosity
       t00    = 0.15_pfdp      ! initial time for exact solution starting from delta function (Gaussian)

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_imexQ_t) :: ad_sweeper_t
!     type(c_ptr) :: ffft, ifft
!     complex(pfdp), pointer :: wk(:)              ! work space
     real(pfdp), allocatable :: wsave(:)           ! work space
     real(pfdp), allocatable :: work(:)             ! work space
     complex(pfdp), allocatable :: workhat(:)      ! work space
     integer ::  lenwrk, lensav, ierror,nvars
     complex(pfdp), allocatable :: ddx(:), lap(:) ! operators
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves 

  end type ad_sweeper_t

contains

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

  subroutine setup(sweeper, nvars)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: nvars

    class(ad_sweeper_t), pointer :: this
    integer     :: i,ierror
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    this => as_ad_sweeper(sweeper)

    ! create in-place, complex fft plans
    
!    wk = fftw_alloc_complex(int(nvars, c_size_t))
!    call c_f_pointer(wk, this%wk, [nvars])
!    this%ffft = fftw_plan_dft_1d(nvars, &
!         this%wk, this%wk, FFTW_FORWARD, FFTW_ESTIMATE)
!    this%ifft = fftw_plan_dft_1d(nvars, &
!         this%wk, this%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    this%nvars = nvars 
    this%lenwrk = 2*nvars 
    this%lensav = 2*nvars + int(log(real(nvars,kind=8))/log(2.0d+00))+4  
    allocate(this%workhat(nvars))   !  complex transform
    allocate(this%work(this%lenwrk))
    allocate(this%wsave(this%lensav))
    call cfft1i ( nvars, this%wsave, this%lensav, ierror )
    if (ierror .ne. 0) then
       stop "error  initializing fft"
    end if

    ! create operators
    allocate(this%ddx(nvars))
    allocate(this%lap(nvars))
    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-nvars + i - 1)
       end if

       this%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          this%lap(i) = 0.0_pfdp
       else
          this%lap(i) = -kx**2
       end if
    end do
  end subroutine setup

  subroutine destroy(this, lev)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev

    deallocate(this%work)
    deallocate(this%workhat)
    deallocate(this%wsave)
    deallocate(this%ddx)
    deallocate(this%lap)
!    call fftw_destroy_plan(this%ffft)
!    call fftw_destroy_plan(this%ifft)
    
    call this%imexQ_destroy(lev)

  end subroutine destroy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0
    call exact(0.0_pfdp, q0%flatarray)
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, yex)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    integer    :: nvars, i, ii, nbox
    real(pfdp) :: tol, x

    nvars = size(yex)
    yex   = 0

    if (nu .gt. 0) then
       do i = 1, nvars
          x = Lx*dble(i-nvars/2-1)/dble(nvars) - t*v
          yex(i) = yex(i) + dsin(2.0_pfdp*pi*x)*dexp(-4.0_pfdp*pi*pi*nu*t)
       end do
    else

       ! decide how many images so that contribution is neglible
       tol  = 1e-16
       !nbox = 1 + ceiling( sqrt( -(4.0*t00)*log((4.0*pi*(t00))**(0.5)*tol) ))
       nbox = 2

       do ii = -nbox, nbox
          do i = 1, nvars
             x = Lx*dble(i-nvars/2-1)/dble(nvars) + ii*Lx - t*v
             yex(i) = yex(i) + 1.0/(4.0*pi*t00)**(0.5)*dexp(-x**2/(4.0*t00))
             x = Lx*dble(i)/dble(nvars)  - t*v
             yex(i) =  dsin(two_pi*kfreq*x)
          end do
       end do
    end if
  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level, f, piece)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:), fvec(:)
    integer ::  ierror

    yvec  => array1(y)
    fvec => array1(f)

    this%workhat = yvec
    !call fftw_execute_dft(this%ffft, wk, wk)

    call cfft1f (this%nvars, 1, this%workhat, this%nvars, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in f_eval"
    end if

    select case (piece)
    case (1)  ! Explicit piece
!      this%workhat = -v * this%ddx *  this%workhat/ real(this%nvars,kind=8)
      this%workhat = -v * this%ddx *  this%workhat
    case (2)  ! Implicit piece
!      this%workhat = nu * this%lap *  this%workhat/ real(this%nvars,kind=8)
      this%workhat = nu * this%lap *  this%workhat
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select

    call cfft1b (this%nvars, 1, this%workhat, this%nvars, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in f_eval"
    end if

    fvec = real(this%workhat)

  end subroutine f_eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dt, rhs, level, f,piece)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dt
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),      pointer :: yvec(:), rhsvec(:), fvec(:)
    integer ::  ierror
    
    if (piece == 2) then
      yvec  => array1(y)
      rhsvec => array1(rhs)
      fvec => array1(f)
      this%workhat = rhsvec
      
      call cfft1f (this%nvars, 1, this%workhat, this%nvars, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in f_comp"
    end if

!      call fftw_execute_dft(this%ffft, wk, wk)
!      this%workhat =  this%workhat/ (1.0_pfdp - nu*dt*this%lap) / real(this%nvars, kind=8)
      this%workhat =  this%workhat/ (1.0_pfdp - nu*dt*this%lap) 
!      call fftw_execute_dft(this%ifft, wk, wk)
      call cfft1b (this%nvars, 1, this%workhat, this%nvars, this%wsave, this%lensav, this%work, this%lenwrk, ierror )      
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in f_comp"
    end if

      yvec  = real(this%workhat)
      fvec = (yvec - rhsvec) / dt
    else
      print *,'Bad piece in f_comp ',piece
      call exit(0)
    end if
  end subroutine f_comp

  subroutine interpolate(this, levelF, levelG, qF, qG, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t


    integer :: nvarF, nvarG, xrat
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),         pointer :: f(:), g(:)
!    complex(kind=8),    pointer :: wkF(:), wkG(:)
    integer ::  ierror
    adG => as_ad_sweeper(levelG%ulevel%sweeper)
    adF => as_ad_sweeper(levelF%ulevel%sweeper)

    f => array1(qF); 
    g => array1(qG)

    nvarF = size(f)
    nvarG = size(g)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       f = g
       return
    endif

    adG%workhat=g
    call cfft1f (adG%nvars, 1, adG%workhat, adG%nvars, adG%wsave, adG%lensav, adG%work, adG%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in interpolate"
    end if

!    call fftw_execute_dft(adG%ffft, wkG, wkG)
!    adG%workhat = adG%workhat / real(nvarG, kind=8)


    adF%workhat = 0.0d0
    adF%workhat(1:nvarG/2) = adG%workhat(1:nvarG/2)
    adF%workhat(nvarF-nvarG/2+2:nvarF) = adG%workhat(nvarG/2+2:nvarG)


!    call fftw_execute_dft(adF%ifft, adF%wk, adF%wk)
    call cfft1b (adF%nvars, 1, adF%workhat, adF%nvars, adF%wsave, adF%lensav, adF%work, adF%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in interpolate"
    end if

    f = real(adF%workhat)

  end subroutine interpolate

  subroutine restrict(this, levelF, levelG, qF, qG, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t

    real(pfdp), pointer :: f(:), g(:)

    integer :: nvarF, nvarG, xrat

    f => array1(qF)
    g => array1(qG)

    nvarF = size(f)
    nvarG = size(g)
    xrat  = nvarF / nvarG

    g = f(::xrat)
  end subroutine restrict

end module feval
