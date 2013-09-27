module feval
  use pf_mod_ndarray
  use probin
  implicit none
  include 'fftw3.f03'

  type :: work1
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)
     complex(pfdp), pointer :: ddx(:), lap(:), ks(:)
  end type work1

  type :: work2
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:,:)
     complex(pfdp), pointer :: ddx(:,:), ddy(:,:), psidx(:,:), psidy(:,:), k2(:,:), w(:,:)
     real(pfdp), pointer :: psi_x(:,:), psi_y(:,:), w_x(:,:), w_y(:,:)
  end type work2


contains

  subroutine f1eval1(yptr, t, level, ctx, f1ptr)
    type(c_ptr),    intent(in), value  :: yptr, f1ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:), f1(:)
    type(work1),   pointer :: work
    complex(pfdp), pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f1 => array1(f1ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = work%ddx * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    select case(problem)
    case (PROB_AD)
       f1 = -v * real(wk)
    case (PROB_VB)
       f1 = -y * real(wk)
    case (PROB_KS)
       f1 = -y * real(wk)
    case default
       stop "ERROR: Unknown problem type (f1eval1)."
    end select
  end subroutine f1eval1

  subroutine f2eval1(yptr, t, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:), f2(:)
    type(work1),   pointer :: work
    complex(pfdp), pointer :: wk(:)

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
    case (PROB_KS)
       wk = work%ks * wk / size(wk)
    case default
       stop "ERROR: Unknown problem type (f2eval1)."
    end select

    call fftw_execute_dft(work%ifft, wk, wk)

    f2 = real(wk)
  end subroutine f2eval1

  subroutine f2comp1(yptr, t, dt, rhsptr, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, rhsptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t, dt
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:), rhs(:), f2(:)
    type(work1),   pointer :: work
    complex(pfdp), pointer :: wk(:)

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
    case (PROB_KS)
       wk = wk / (1.0_pfdp - dt*work%ks) / size(wk)
    case default
       stop "ERROR: Unknown problem type (f2comp1)."
    end select

    call fftw_execute_dft(work%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp1

  subroutine f1eval2(yptr, t, level, ctx, f1ptr)
    type(c_ptr),    intent(in), value  :: yptr, f1ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:,:), f1(:,:)
    type(work2),   pointer :: work
    complex(pfdp), pointer :: wk(:,:)

    call c_f_pointer(ctx, work)

    y  => array2(yptr)
    f1 => array2(f1ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    work%w = wk / size(wk)

    wk = real(work%psidx * work%w)
    call fftw_execute_dft(work%ifft, wk, wk)
    work%psi_x = real(wk)

    wk = real(work%psidy * work%w)
    call fftw_execute_dft(work%ifft, wk, wk)
    work%psi_y = real(wk)

    wk = real(work%ddx * work%w)
    call fftw_execute_dft(work%ifft, wk, wk)
    work%w_x = real(wk)

    wk = real(work%ddy * work%w)
    call fftw_execute_dft(work%ifft, wk, wk)
    work%w_y = real(wk)

    f1 = -(work%psi_y * work%w_x - work%psi_x * work%w_y)
  end subroutine f1eval2

  subroutine f2eval2(yptr, t, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:,:), f2(:,:)
    type(work2),   pointer :: work
    complex(pfdp), pointer :: wk(:,:)

    call c_f_pointer(ctx, work)

    y  => array2(yptr)
    f2 => array2(f2ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = -nu * work%k2 * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)
    f2 = real(wk)
  end subroutine f2eval2

  subroutine f2comp2(yptr, t, dt, rhsptr, level, ctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, rhsptr, f2ptr, ctx
    real(pfdp),     intent(in)         :: t, dt
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:,:), rhs(:,:), f2(:,:)
    type(work2),   pointer :: work
    complex(pfdp), pointer :: wk(:,:)

    call c_f_pointer(ctx, work)

    y   => array2(yptr)
    rhs => array2(rhsptr)
    f2  => array2(f2ptr)
    wk  => work%wk

    wk = rhs
    call fftw_execute_dft(work%ffft, wk, wk)

    wk = wk / (1.0_pfdp + nu*dt*work%k2) / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp2

  subroutine create_work1(ctx, nvars)
    type(c_ptr), intent(out) :: ctx
    integer,     intent(in)  :: nvars

    type(work1), pointer :: work
    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx(nvars)

    allocate(work)
    ctx = c_loc(work)

    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx(i) = two_pi * dble(i-1) / Lx
       else
          kx(i) = two_pi * dble(-nvars + i - 1) / Lx
       end if
    end do

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
    allocate(work%ks(nvars))

    work%ddx = (0.0_pfdp, 1.0_pfdp) * kx
    work%lap = -kx**2
    work%ks  = kx**2*(1.0 - kx**2)

    where (kx**2 < 1e-13)
       work%lap = 0.0_pfdp
       work%ks  = 0.0_pfdp
    end where

  end subroutine create_work1

  subroutine create_work2(ctx, nvars)
    type(c_ptr), intent(out) :: ctx
    integer,     intent(in)  :: nvars

    type(work2), pointer :: work
    integer     :: i, j
    type(c_ptr) :: wk
    real(pfdp)  :: kx(nvars)

    allocate(work)
    ctx = c_loc(work)

    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx(i) = two_pi * dble(i-1) / Lx
       else
          kx(i) = two_pi * dble(-nvars + i - 1) / Lx
       end if
    end do

    ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nvars**2, c_size_t))
    call c_f_pointer(wk, work%wk, [nvars,nvars])
    
    work%ffft = fftw_plan_dft_2d(nvars, nvars, &
         work%wk, work%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    work%ifft = fftw_plan_dft_2d(nvars, nvars, &
         work%wk, work%wk, FFTW_BACKWARD, FFTW_ESTIMATE)
    
    ! create operators
    allocate(work%k2(nvars,nvars))
    allocate(work%ddx(nvars,nvars))
    allocate(work%ddy(nvars,nvars))
    allocate(work%psidx(nvars,nvars))
    allocate(work%psidy(nvars,nvars))
    allocate(work%psi_x(nvars,nvars))
    allocate(work%psi_y(nvars,nvars))
    allocate(work%w_x(nvars,nvars))
    allocate(work%w_y(nvars,nvars))
    allocate(work%w(nvars,nvars))
    
    work%psidx = 0
    work%psidy = 0
    
    do j = 1, nvars
       do i = 1, nvars
          work%k2(i,j) = kx(i)**2 + kx(j)**2
          work%ddx(i,j) = (0.0_pfdp, 1.0_pfdp) * kx(i)
          work%ddy(i,j) = (0.0_pfdp, 1.0_pfdp) * kx(j)
    
          if (work%k2(i,j) /= 0.0_pfdp) then
             work%psidx(i,j) = (0.0_pfdp, 1.0_pfdp) * kx(i) / work%k2(i,j)
             work%psidy(i,j) = (0.0_pfdp, 1.0_pfdp) * kx(j) / work%k2(i,j)
          end if
       end do
    end do

  end subroutine create_work2


  subroutine destroy_work1(ctx)
    type(c_ptr), intent(in) :: ctx
    type(work1), pointer :: work
    call c_f_pointer(ctx, work)
    deallocate(work%wk)
    deallocate(work%ddx)
    deallocate(work%lap)
    deallocate(work%ks)
    call fftw_destroy_plan(work%ffft)
    call fftw_destroy_plan(work%ifft)
    deallocate(work)
  end subroutine destroy_work1

end module feval
