module feval
  use pf_mod_ndarray
  use probin
  implicit none
  include 'fftw3.f03'

  type :: ad_work_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk1(:)                   ! work space
     complex(pfdp), pointer :: ddx1(:), lap1(:), ks1(:) ! operators
  end type ad_work_t

contains

  subroutine f1eval1(yptr, t, level, ctx, f1ptr)
    type(c_ptr),    intent(in), value  :: yptr, f1ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),      pointer :: y(:), f1(:)
    type(ad_work_t), pointer :: work
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f1 => array1(f1ptr)
    wk => work%wk1

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = work%ddx1 * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    select case(problem)
    case (PROB_AD)
       f1 = -v * real(wk)
    case (PROB_VB)
       f1 = -y * real(wk)
    case (PROB_KS)
       f1 = -y * real(wk)
    case default
       stop "ERROR: Unknown problem type in f1eval."
    end select
  end subroutine f1eval1

  subroutine f2eval1(yptr, t, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),      pointer :: y(:), f2(:)
    type(ad_work_t), pointer :: work
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f2 => array1(f2ptr)
    wk => work%wk1

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    
    select case(problem)
    case (PROB_HEAT)
       wk = nu * work%lap1 * wk / size(wk)
    case (PROB_AD)
       wk = nu * work%lap1 * wk / size(wk)
    case (PROB_VB)
       wk = nu * work%lap1 * wk / size(wk)
    case (PROB_KS)
       wk = work%ks1 * wk / size(wk)
    case default
       stop "ERROR: Unknown problem type in f2eval."
    end select

    call fftw_execute_dft(work%ifft, wk, wk)

    f2 = real(wk)
  end subroutine f2eval1

  subroutine f2comp1(yptr, t, dt, rhsptr, level, ctx, f2ptr)
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
    wk  => work%wk1

    wk = rhs
    call fftw_execute_dft(work%ffft, wk, wk)

    select case(problem)
    case (PROB_HEAT)
       wk = wk / (1.0_pfdp - nu*dt*work%lap1) / size(wk)
    case (PROB_AD)
       wk = wk / (1.0_pfdp - nu*dt*work%lap1) / size(wk)
    case (PROB_VB)
       wk = wk / (1.0_pfdp - nu*dt*work%lap1) / size(wk)
    case (PROB_KS)
       wk = wk / (1.0_pfdp - dt*work%ks1) / size(wk)
    case default
       stop "ERROR: Unknown problem type in f2comp."
    end select

    call fftw_execute_dft(work%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp1

  subroutine f1eval2(yptr, t, level, ctx, f1ptr)
    type(c_ptr),    intent(in), value  :: yptr, f1ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),      pointer :: y(:,:), f1(:,:)
    type(ad_work_t), pointer :: work
    complex(pfdp),   pointer :: wk(:,:)

    call c_f_pointer(ctx, work)

    y  => array2(yptr)
    f1 => array2(f1ptr)

  end subroutine f1eval2

  subroutine f2eval2(yptr, t, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),      pointer :: y(:), f2(:)
    type(ad_work_t), pointer :: work
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f2 => array1(f2ptr)

  end subroutine f2eval2

  subroutine f2comp2(yptr, t, dt, rhsptr, level, ctx, f2ptr)
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

  end subroutine f2comp2

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
    call c_f_pointer(wk, work%wk1, [nvars])

    work%ffft = fftw_plan_dft_1d(nvars, &
         work%wk1, work%wk1, FFTW_FORWARD, FFTW_ESTIMATE)
    work%ifft = fftw_plan_dft_1d(nvars, &
         work%wk1, work%wk1, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
    allocate(work%ddx1(nvars))
    allocate(work%lap1(nvars))
    allocate(work%ks1(nvars))
    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx = two_pi * dble(i-1) / Lx
       else
          kx = two_pi * dble(-nvars + i - 1) / Lx
       end if

       work%ddx1(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          work%lap1(i) = 0.0_pfdp
          work%ks1(i)  = 0.0_pfdp
       else
          work%lap1(i) = -kx**2
          work%ks1(i)  = kx**2*(1.0 - kx**2)
       end if
    end do
  end subroutine feval_create_workspace

  subroutine feval_destroy_workspace(ctx)
    type(c_ptr), intent(in) :: ctx

    type(ad_work_t), pointer :: work

    call c_f_pointer(ctx, work)

    deallocate(work%wk1)
    deallocate(work%ddx1)
    deallocate(work%lap1)
    deallocate(work%ks1)
    call fftw_destroy_plan(work%ffft)
    call fftw_destroy_plan(work%ifft)
    deallocate(work)
  end subroutine feval_destroy_workspace

end module feval
