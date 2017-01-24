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
       nu     = 0.0_pfdp, &  ! viscosity
       t00    = 0.15_pfdp    ! initial time for exact solution

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type, extends(pf_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_imexQ_t) :: ad_sweeper_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), allocatable :: ddx(:), lap(:)     ! operators
   contains
     procedure :: setup
     procedure :: f1eval
     procedure :: f2eval
     procedure :: f2comp
     ! procedure :: restrict
     ! procedure :: interpolate
     final :: destroy0, destroy1
  end type ad_sweeper_t

contains

  subroutine setup(this, nvars)
    class(ad_sweeper_t), intent(inout) :: this
    integer,             intent(in   ) :: nvars

    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx

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
  subroutine f1eval(this, y, t, level, f1)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f1
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level

    real(pfdp),      pointer :: yvec(:), f1vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    f1vec => array1(f1)
    wk    => this%wk

    wk = yvec

    call fftw_execute_dft(this%ffft, wk, wk)
    wk = -v * this%ddx * wk / size(wk)
    call fftw_execute_dft(this%ifft, wk, wk)

    f1vec = real(wk)
  end subroutine f1eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine f2eval(this, y, t, level, f2)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f2
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level

    real(pfdp),      pointer :: yvec(:), f2vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    f2vec => array1(f2)
    wk => this%wk

    wk = yvec
    call fftw_execute_dft(this%ffft, wk, wk)
    wk = nu * this%lap * wk / size(wk)
    call fftw_execute_dft(this%ifft, wk, wk)

    f2vec = real(wk)

  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine f2comp(this, y, t, dt, rhs, level, f2)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: rhs
    class(pf_encap_t),   intent(inout) :: y, f2
    real(pfdp),          intent(in   ) :: t, dt
    integer,             intent(in   ) :: level

    real(pfdp),      pointer :: yvec(:), rhsvec(:), f2vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    rhsvec => array1(rhs)
    f2vec => array1(f2)
    wk => this%wk

    wk = rhsvec
    call fftw_execute_dft(this%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*this%lap) / size(wk)
    call fftw_execute_dft(this%ifft, wk, wk)

    yvec  = real(wk)
    f2vec = (yvec - rhsvec) / dt

  end subroutine f2comp

  subroutine interpolate(levelF, levelG, qFp, qGp, t)
    class(ad_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qFp, qGp
    real(pfdp),        intent(in   ) :: t

    real(pfdp),      pointer :: qF(:), qG(:)
    complex(kind=8), pointer :: wkF(:), wkG(:)

    integer :: nvarF, nvarG, xrat

    type(ad_sweeper_t), pointer :: adF, adG

    select type(swp => levelG%sweeper)
    type is (ad_sweeper_t)
       adG => swp
    class default
       stop
    end select

    select type(swp => levelF%sweeper)
    type is (ad_sweeper_t)
       adF => swp
    class default
       stop
    end select

    qF => array1(qFp)
    qG => array1(qGp)

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       qF = qG
       return
    endif

    wkF => adF%wk
    wkG => adG%wk

    wkG = qG
    call fftw_execute_dft(adG%ffft, wkG, wkG)
    wkG = wkG / nvarG

    wkF = 0.0d0
    wkF(1:nvarG/2) = wkG(1:nvarG/2)
    wkF(nvarF-nvarG/2+2:nvarF) = wkG(nvarG/2+2:nvarG)

    call fftw_execute_dft(adF%ifft, wkF, wkF)

    qF = real(wkF)
  end subroutine interpolate

  subroutine restrict(levelF, levelG, qFp, qGp, t)
    class(ad_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qFp, qGp
    real(pfdp),        intent(in   ) :: t

    real(pfdp), pointer :: qF(:), qG(:)

    integer :: nvarF, nvarG, xrat

    qF => array1(qFp)
    qG => array1(qGp)

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    qG = qF(::xrat)
  end subroutine restrict

end module feval
