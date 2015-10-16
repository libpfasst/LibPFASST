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

  type :: ad_work_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), pointer :: ddx(:), lap(:)     ! operators
  end type ad_work_t

  type, extends(pf_imexQ_t) :: ad_sweeper_t
     type(ad_work_t) :: work
   contains
     procedure :: f1eval      => eval_f1
     procedure :: f2eval      => eval_f2
     procedure :: f2comp      => comp_f2
     procedure :: restrict    => restrict
     procedure :: interpolate => interpolate
  end type ad_sweeper_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine feval_create_workspace(work, nvars)
    type(ad_work_t), intent(out) :: work
    integer,         intent(in)  :: nvars

    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nvars, c_size_t))
    call c_f_pointer(wk, work%wk, [nvars])

    work%ffft = fftw_plan_dft_1d(nvars, &
         work%wk, work%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    work%ifft = fftw_plan_dft_1d(nvars, &
         work%wk, work%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
    allocate(work%ddx(nvars))
    allocate(work%lap(nvars))
    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-nvars + i - 1)
       end if

       work%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          work%lap(i) = 0.0_pfdp
       else
          work%lap(i) = -kx**2
       end if
    end do
  end subroutine feval_create_workspace

  ! subroutine feval_destroy_workspace(levelctx)
  !   type(c_ptr), intent(in) :: levelctx
  !   type(ad_work_t), pointer :: work

  !   call c_f_pointer(levelctx, work)

  !   deallocate(work%wk)
  !   deallocate(work%ddx)
  !   deallocate(work%lap)
  !   call fftw_destroy_plan(work%ffft)
  !   call fftw_destroy_plan(work%ifft)
  !   deallocate(work)
  ! end subroutine feval_destroy_workspace

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
       nbox = 1 + ceiling( sqrt( -(4.0*t00)*log((4.0*pi*(t00))**(0.5)*tol) ))

!       do ii = -nbox, nbox
          do i = 1, nvars
             x = Lx*dble(i-nvars/2-1)/dble(nvars) + ii*Lx - t*v
             yex(i) = yex(i) + 1.0/(4.0*pi*t00)**(0.5)*dexp(-x**2/(4.0*t00))
             x = Lx*dble(i)/dble(nvars)  - t*v
             yex(i) =  dsin(two_pi*kfreq*x)
          end do
 !      end do
    end if


  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine eval_f1(this, y, t, level, f1)
    class(ad_sweeper_t), intent(inout)     :: this
    type(c_ptr),         intent(in), value :: y, f1
    real(pfdp),          intent(in)        :: t
    integer,             intent(in)        :: level

    real(pfdp),      pointer :: yvec(:), f1vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    f1vec => array1(f1)
    wk => this%work%wk

    wk = yvec
    call fftw_execute_dft(this%work%ffft, wk, wk)
    wk = -v * this%work%ddx * wk / size(wk)
    call fftw_execute_dft(this%work%ifft, wk, wk)

    f1vec = real(wk)

  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(this, y, t, level, f2)
    class(ad_sweeper_t), intent(inout)     :: this
    type(c_ptr),         intent(in), value :: y, f2
    real(pfdp),          intent(in)        :: t
    integer,             intent(in)        :: level

    real(pfdp),      pointer :: yvec(:), f2vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    f2vec => array1(f2)
    wk => this%work%wk

    wk = yvec
    call fftw_execute_dft(this%work%ffft, wk, wk)
    wk = nu * this%work%lap * wk / size(wk)
    call fftw_execute_dft(this%work%ifft, wk, wk)

    f2vec = real(wk)

  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine comp_f2(this, y, t, dt, rhs, level, f2)
    class(ad_sweeper_t), intent(inout)     :: this
    type(c_ptr),         intent(in), value :: y, rhs, f2
    real(pfdp),          intent(in)        :: t, dt
    integer,             intent(in)        :: level

    real(pfdp),      pointer :: yvec(:), rhsvec(:), f2vec(:)
    complex(pfdp),   pointer :: wk(:)

    yvec  => array1(y)
    rhsvec => array1(rhs)
    f2vec => array1(f2)
    wk => this%work%wk

    wk = rhsvec
    call fftw_execute_dft(this%work%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*this%work%lap) / size(wk)
    call fftw_execute_dft(this%work%ifft, wk, wk)

    yvec  = real(wk)
    f2vec = (yvec - rhsvec) / dt

  end subroutine comp_f2

  subroutine interpolate(levelF, levelG, qFp, qGp, t)
    class(ad_sweeper_t), intent(inout) :: levelF
    class(pf_sweeper_t), intent(inout) :: levelG
    type(c_ptr), intent(in), value :: qFp, qGp
    real(pfdp),  intent(in) :: t

    real(pfdp),      pointer :: qF(:), qG(:)
    complex(kind=8), pointer :: wkF(:), wkG(:)

    integer :: nvarF, nvarG, xrat

    type(ad_sweeper_t), pointer :: levelGad

    select type(levelG)
    type is (ad_sweeper_t)
       levelGad => levelG
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

    wkF => levelF%work%wk
    wkG => levelGad%work%wk

    wkG = qG
    call fftw_execute_dft(levelGad%work%ffft, wkG, wkG)
    wkG = wkG / nvarG

    wkF = 0.0d0
    wkF(1:nvarG/2) = wkG(1:nvarG/2)
    wkF(nvarF-nvarG/2+2:nvarF) = wkG(nvarG/2+2:nvarG)

    call fftw_execute_dft(levelF%work%ifft, wkF, wkF)

    qF = real(wkF)
  end subroutine interpolate

  subroutine restrict(levelF, levelG, qFp, qGp, t)
    class(ad_sweeper_t), intent(inout) :: levelF
    class(pf_sweeper_T), intent(inout) :: levelG
    type(c_ptr), intent(in), value :: qFp, qGp
    real(pfdp),  intent(in) :: t

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
