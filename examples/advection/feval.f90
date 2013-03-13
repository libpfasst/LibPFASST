!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use iso_c_binding
  use encap
  implicit none
  include 'fftw3.f03'

  real(pfdp), parameter :: &
       Lx     = 1.0_pfdp, &        ! domain size
       v      = 1.0_pfdp, &        ! velocity
       nu     = 0.02_pfdp, &       ! viscosity
       t00    = 0.25_pfdp           ! initial time for exact solution

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type :: ad_level_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), pointer :: ddx(:), lap(:)     ! operators
  end type ad_level_t

  type(ad_level_t), pointer :: levels(:)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine feval_init(nlevels, nvars)

    integer, intent(in) :: nlevels, nvars(nlevels)

    integer :: i, l
    type(c_ptr) :: wk
    real(pfdp) :: kx

    allocate(levels(nlevels))
    do l = 1, nlevels

       ! create in-place, complex fft plans
       wk = fftw_alloc_complex(int(nvars(l), c_size_t))
       call c_f_pointer(wk, levels(l)%wk, [nvars(l)])

       levels(l)%ffft = fftw_plan_dft_1d(nvars(l), &
            levels(l)%wk, levels(l)%wk, FFTW_FORWARD, FFTW_ESTIMATE)
       levels(l)%ifft = fftw_plan_dft_1d(nvars(l), &
            levels(l)%wk, levels(l)%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

       ! create operators
       allocate(levels(l)%ddx(nvars(l)))
       allocate(levels(l)%lap(nvars(l)))
       do i = 1, nvars(l)
          if (i <= nvars(l)/2+1) then
             kx = two_pi / Lx * dble(i-1)
          else
             kx = two_pi / Lx * dble(-nvars(l) + i - 1)
          end if

          levels(l)%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

          if (kx**2 < 1e-13) then
             levels(l)%lap(i) = 0.0_pfdp
          else
             levels(l)%lap(i) = -kx**2
          end if
       end do
    end do
  end subroutine feval_init

  subroutine feval_finalize()

    integer :: l

    do l = 1, size(levels)
       deallocate(levels(l)%wk)
       deallocate(levels(l)%ddx)
       deallocate(levels(l)%lap)
    end do

    deallocate(levels)

  end subroutine feval_finalize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    type(pf_encap_t), intent(inout) :: q0

    call exact(0.0_pfdp, size(q0%array), q0%array)
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, nvars, yex)
    real(pfdp), intent(in)  :: t
    integer,      intent(in)  :: nvars
    real(pfdp), intent(out) :: yex(nvars)

    integer :: i, ii, nbox
    real(pfdp) :: x,tol

    yex = 0.0_pfdp

    ! decide how many images so that contribution is neglible
    tol = 1e-16
    if (nu .gt. ZERO) then
       nbox = 1+ceiling( sqrt( -(4.0_pfdp*nu*(t+t00))*log((4.0*pi*nu*(t+t00))**(HALF)*tol) ))

!       do ii = -nbox, nbox
          ii = 0
          do i = 1, nvars
             x = Lx*dble(i-nvars/2-1)/dble(nvars) + dble(ii)*Lx - t*v
!             yex(i) = yex(i) + ONE/(4.0_pfdp*pi*nu*(t+t00))**(0.5)*dexp(-x**2/(4.0_pfdp*nu*(t+t00)))
             yex(i) = yex(i) + dcos(2.0_pfdp*pi*x)*dexp(-4.0_pfdp*pi*pi*nu*t)
          end do
          
 !      end do
    else
       nbox = 1+ceiling( sqrt( -(4.0*t00)*log((4.0*pi*(t00))**(0.5)*tol) ))
       do ii = -nbox, nbox
          do i = 1, nvars
             x = Lx*dble(i-nvars/2-1)/dble(nvars) + ii*Lx - t*v
             yex(i) = yex(i) + 1.0/(4.0*pi*t00)**(0.5)*dexp(-x**2/(4.0*t00))
          end do
          
       end do
    end if


  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine eval_f1(y, t, level, ctx, f1)

    type(pf_encap_t), intent(in)    :: y
    real(pfdp),     intent(in)    :: t
    integer,          intent(in)    :: level
    type(pf_encap_t), intent(inout) :: f1
    type(c_ptr),      intent(in)    :: ctx

    complex(pfdp), pointer :: wk(:)

    wk => levels(level)%wk

    wk = y%array
    call fftw_execute_dft(levels(level)%ffft, wk, wk)
    wk = -v * levels(level)%ddx * wk / size(wk)
    call fftw_execute_dft(levels(level)%ifft, wk, wk)

    f1%array = real(wk)

  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(y, t, level, ctx, f2)

    type(pf_encap_t), intent(in)    :: y
    real(pfdp),     intent(in)    :: t
    integer,          intent(in)    :: level
    type(pf_encap_t), intent(inout) :: f2
    type(c_ptr),      intent(in)    :: ctx


    complex(pfdp), pointer :: wk(:)

    wk => levels(level)%wk

    wk = y%array
    call fftw_execute_dft(levels(level)%ffft, wk, wk)
    wk = nu * levels(level)%lap * wk / size(wk)
    call fftw_execute_dft(levels(level)%ifft, wk, wk)

    f2%array = real(wk)

  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine comp_f2(y, t, dt, rhs, level, ctx, f2)

    type(pf_encap_t), intent(inout) :: y
    type(pf_encap_t), intent(in)    :: rhs
    real(pfdp),     intent(in)    :: t, dt
    integer,          intent(in)    :: level
    type(pf_encap_t), intent(inout) :: f2
    type(c_ptr),      intent(in)    :: ctx

    complex(pfdp), pointer :: wk(:)

    wk => levels(level)%wk

    wk = rhs%array
    call fftw_execute_dft(levels(level)%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*levels(level)%lap) / size(wk)
    call fftw_execute_dft(levels(level)%ifft, wk, wk)

    y%array  = real(wk)
    f2%array = (y%array - rhs%array) / dt

  end subroutine comp_f2

end module feval
