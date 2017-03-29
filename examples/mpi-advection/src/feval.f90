!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imex
  use pf_mod_imexQ
  implicit none
  include 'fftw3.f03'

  real(pfdp), parameter :: &
       Lx     = 1.0_pfdp, &  ! domain size
       kfreq  = 32.0_pfdp, & ! domain size
       v      = 1.0_pfdp, &  ! velocity
       nu     = 0.02_pfdp, &  ! viscosity
       t00    = 0.15_pfdp    ! initial time for exact solution

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_imexQ_t) :: ad_sweeper_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), allocatable :: ddx(:), lap(:) ! operators
   contains
     procedure :: f1eval
     procedure :: f2eval
     procedure :: f2comp
!     final :: destroy0, destroy1
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
    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    this => as_ad_sweeper(sweeper)

    ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nvars, c_size_t))
    call c_f_pointer(wk, this%wk, [nvars])

    this%ffft = fftw_plan_dft_1d(nvars, &
         this%wk, this%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    this%ifft = fftw_plan_dft_1d(nvars, &
         this%wk, this%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

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

  subroutine destroy0(this)
    type(ad_sweeper_t), intent(inout) :: this
    deallocate(this%wk)
    deallocate(this%ddx)
    deallocate(this%lap)
    call fftw_destroy_plan(this%ffft)
    call fftw_destroy_plan(this%ifft)
  end subroutine destroy0

  subroutine destroy1(this)
    type(ad_sweeper_t), intent(inout) :: this(:)
    integer :: i
    do i = 1, size(this)
       call destroy0(this(i))
    end do
  end subroutine destroy1

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
          yex(i) = yex(i) + dcos(2.0_pfdp*pi*x)*dexp(-4.0_pfdp*pi*pi*nu*t)
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
  subroutine f1eval(this, y, t, level, f)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level

    real(pfdp),      pointer :: yvec(:), f1vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    f1vec => array1(f)
    wk    => this%wk

    wk = yvec

    call fftw_execute_dft(this%ffft, wk, wk)
    wk = -v * this%ddx * wk / size(wk)
    call fftw_execute_dft(this%ifft, wk, wk)

    f1vec = real(wk)
  end subroutine f1eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine f2eval(this, y, t, level, f)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level

    real(pfdp),      pointer :: yvec(:), f2vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    f2vec => array1(f)
    wk => this%wk

    wk = yvec
    call fftw_execute_dft(this%ffft, wk, wk)
    wk = nu * this%lap * wk / size(wk)
    call fftw_execute_dft(this%ifft, wk, wk)

    f2vec = real(wk)

  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine f2comp(this, y, t, dt, rhs, level, f)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: rhs
    class(pf_encap_t),   intent(inout) :: y, f
    real(pfdp),          intent(in   ) :: t, dt
    integer,             intent(in   ) :: level

    real(pfdp),      pointer :: yvec(:), rhsvec(:), f2vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    rhsvec => array1(rhs)
    f2vec => array1(f)
    wk => this%wk

    wk = rhsvec
    call fftw_execute_dft(this%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*this%lap) / size(wk)
    call fftw_execute_dft(this%ifft, wk, wk)

    yvec  = real(wk)
    f2vec = (yvec - rhsvec) / dt

  end subroutine f2comp

  subroutine interpolate(this, levelF, levelG, qF, qG, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t


    integer :: nvarF, nvarG, xrat
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),         pointer :: f(:), g(:)
    complex(kind=8),    pointer :: wkF(:), wkG(:)

    adG => as_ad_sweeper(levelG%ulevel%sweeper)
    adF => as_ad_sweeper(levelF%ulevel%sweeper)

    f => array1(qF); g => array1(qG)

    nvarF = size(f)
    nvarG = size(g)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       f = g
       return
    endif

    wkF => adF%wk; wkG => adG%wk

    wkG = g
    call fftw_execute_dft(adG%ffft, wkG, wkG)
    wkG = wkG / nvarG

    adF%wk = 0.0d0
    adF%wk(1:nvarG/2) = wkG(1:nvarG/2)
    adF%wk(nvarF-nvarG/2+2:nvarF) = wkG(nvarG/2+2:nvarG)

    call fftw_execute_dft(adF%ifft, adF%wk, adF%wk)

    f = real(adF%wk)

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
