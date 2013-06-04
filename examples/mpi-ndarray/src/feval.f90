!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use iso_c_binding
  use encap
  use probin
  implicit none
  include 'fftw3.f03'

  type, extends(pf_context_t) :: ad_work_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), pointer :: ddx(:), lap(:)     ! operators
  end type ad_work_t

contains

  ! 
  ! Evaluation routines
  !

  subroutine f1eval(yabs, t, level, ctx, f1abs)
    class(pf_encap_t),   intent(inout), target :: yabs, f1abs
    class(pf_context_t), intent(inout), target :: ctx
    real(pfdp),          intent(in)            :: t
    integer,             intent(in)            :: level

    real(pfdp),       pointer :: y(:), f1(:)
    class(ad_work_t), pointer :: work
    complex(pfdp),    pointer :: wk(:)

    work => get_work(ctx)
    y    => get_array(yabs)
    f1   => get_array(f1abs)
    wk   => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)

    select case(problem)
    case (PROB_AD)
       wk = -v * work%ddx * wk / size(wk)
    case (PROB_VB)
       wk = -wk * work%ddx * wk / size(wk)
    case (PROB_KS)
       wk = -wk * work%ddx * wk / size(wk)
    case default
       stop "ERROR: Unknown problem type in f1eval."
    end select

    call fftw_execute_dft(work%ifft, wk, wk)

    f1 = real(wk)
  end subroutine f1eval

  subroutine f2eval(yabs, t, level, ctx, f2abs)
    class(pf_encap_t),   intent(inout), target :: yabs, f2abs
    class(pf_context_t), intent(inout), target :: ctx
    real(pfdp),          intent(in)            :: t
    integer,             intent(in)            :: level

    real(pfdp),       pointer :: y(:), f2(:)
    class(ad_work_t), pointer :: work
    complex(pfdp),    pointer :: wk(:)

    work => get_work(ctx)
    y    => get_array(yabs)
    f2   => get_array(f2abs)
    wk   => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    
    select case(problem)
    case (PROB_AD)
       wk = nu * work%lap * wk / size(wk)
    case (PROB_VB)
       wk = nu * work%lap * wk / size(wk)
    case default
       stop "ERROR: Unknown problem type in f2eval."
    end select

    call fftw_execute_dft(work%ifft, wk, wk)

    f2 = real(wk)
  end subroutine f2eval

  subroutine f2comp(yabs, t, dt, rhsabs, level, ctx, f2abs)
    class(pf_encap_t),   intent(inout), target :: yabs, rhsabs, f2abs
    class(pf_context_t), intent(inout), target :: ctx
    real(pfdp),          intent(in)            :: t, dt
    integer,             intent(in)            :: level

    real(pfdp),       pointer :: y(:), rhs(:), f2(:)
    class(ad_work_t), pointer :: work
    complex(pfdp),    pointer :: wk(:)

    work => get_work(ctx)
    y    => get_array(yabs)
    rhs  => get_array(rhsabs)
    f2   => get_array(f2abs)
    wk   => work%wk

    wk = rhs
    call fftw_execute_dft(work%ffft, wk, wk)

    select case(problem)
    case (PROB_AD)
       wk = wk / (1.0_pfdp - nu*dt*work%lap) / size(wk)
    case (PROB_VB)
       wk = wk / (1.0_pfdp - nu*dt*work%lap) / size(wk)
    case default
       stop "ERROR: Unknown problem type in f2comp."
    end select

    call fftw_execute_dft(work%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp


  !
  ! Helpers
  !

  function get_work(ctx) result(r)
    class(pf_context_t), intent(in), target :: ctx
    class(ad_work_t), pointer :: r
    select type(ctx)
    class is (ad_work_t)
       r => ctx
    end select
  end function get_work


  !
  ! Create/destroy work space
  !

  subroutine feval_create_workspace(ctx, nvars)
    class(pf_context_t), intent(out), pointer :: ctx
    integer,             intent(in)           :: nvars

    class(ad_work_t), pointer :: work
    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    allocate(ad_work_t::ctx)
    work => get_work(ctx)

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

  subroutine feval_destroy_workspace(ctx)
    class(pf_context_t), intent(inout) :: ctx

    class(ad_work_t), pointer :: work
    work => get_work(ctx)

    deallocate(work%wk)
    deallocate(work%ddx)
    deallocate(work%lap)
    call fftw_destroy_plan(work%ffft)
    call fftw_destroy_plan(work%ifft)
    deallocate(work)
  end subroutine feval_destroy_workspace

end module feval
