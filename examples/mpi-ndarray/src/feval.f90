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

  type :: ad_work_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), pointer :: ddx(:), lap(:)     ! operators
  end type ad_work_t

contains

  ! 
  ! Evaluation routines
  !

  subroutine f1eval(yptr, t, level, ctx, f1ptr)
    type(c_ptr),    intent(in), value  :: yptr, f1ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),      pointer :: y(:), f1(:)
    type(ad_work_t), pointer :: work
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f1 => array1(f1ptr)
    wk => work%wk

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

  subroutine f2eval(yptr, t, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),      pointer :: y(:), f2(:)
    type(ad_work_t), pointer :: work
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f2 => array1(f2ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    
    select case(problem)
    case (PROB_HEAT)
       wk = nu * work%lap * wk / size(wk)
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

  subroutine f2comp(yptr, t, dt, rhsptr, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, rhsptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t, dt
    integer(c_int), intent(in)         :: level

    real(pfdp),      pointer :: y(:), rhs(:), f2(:)
    type(ad_work_t), pointer :: work
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y   => array1(yptr)
    rhs => array1(rhsptr)
    f2  => array1(f2ptr)
    wk  => work%wk

    wk = rhs
    call fftw_execute_dft(work%ffft, wk, wk)

    select case(problem)
    case (PROB_HEAT)
       wk = wk / (1.0_pfdp - nu*dt*work%lap) / size(wk)
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
  ! Create/destroy work space
  !

  subroutine feval_create_workspace(ctx, nvars)
    type(c_ptr), intent(out) :: ctx
    integer,     intent(in)  :: nvars

    type(ad_work_t), pointer :: work
    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    allocate(work)
    ctx = c_loc(work)

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
          kx = two_pi * dble(i-1) / Lx
       else
          kx = two_pi * dble(-nvars + i - 1) / Lx
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
    type(c_ptr), intent(in) :: ctx

    type(ad_work_t), pointer :: work

    call c_f_pointer(ctx, work)

    deallocate(work%wk)
    deallocate(work%ddx)
    deallocate(work%lap)
    call fftw_destroy_plan(work%ffft)
    call fftw_destroy_plan(work%ifft)
    deallocate(work)
  end subroutine feval_destroy_workspace

end module feval
