module feval
  use pf_mod_ndarray
  use probin
  implicit none
  include 'fftw3.f03'

  type :: work1
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)
     complex(pfdp), pointer :: ddx(:), lap(:)
  end type work1

  type :: work2
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:,:)
     complex(pfdp), pointer :: ddx(:,:), ddy(:,:), lap(:,:)
  end type work2

  type :: work3
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:,:,:)
     complex(pfdp), pointer :: ddx(:,:,:), ddy(:,:,:), ddz(:,:,:), lap(:,:,:)
  end type work3

contains

  !
  ! one dimension
  !

  subroutine f1eval1(yptr, t, level, levelctx, f1ptr)
    use probin, only: nprob, do_spec
    use mg, only: fill_bc_1d
    type(c_ptr), intent(in), value :: yptr, f1ptr, levelctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level

    type(work1),   pointer :: work
    real(pfdp),    pointer :: y(:), f1(:)
    complex(pfdp), pointer :: wk(:)
    real(pfdp), allocatable :: ybc(:)
    real(pfdp) :: dx
    integer :: Nx,i


    y  => array1(yptr)
    f1 => array1(f1ptr)
    Nx = size(y)
    if (do_spec==1) then
       wk => work%wk
       call c_f_pointer(levelctx, work)
       wk = y
       call fftw_execute_dft(work%ffft, wk, wk)
       if (nprob .eq. 0) then  ! spectral advection diffusion
          wk = -v * work%ddx * wk / size(wk) 
       else if (nprob .eq. 1) then  ! Backward Euler
          wk = 0.0_pfdp*nu * work%lap * wk / size(wk)
       else if (nprob .eq. 2) then  ! trapezoid rule
          wk = 0.5_pfdp*nu * work%lap * wk / size(wk)
       end if
       call fftw_execute_dft(work%ifft, wk, wk)
       f1 = real(wk)
    else
       dx = 1.0_pfdp/dble(Nx)
       allocate(ybc(1-2:Nx+2))
       call fill_bc_1d(y,ybc,Nx,2)
       if (spatial_order == 2) then
          do i = 1,Nx
             f1(i)=-v*(ybc(i+1)-ybc(i-1))/(2.0_pfdp*dx)
          end do
       else
          do i = 1,Nx
             f1(i)=-v*(-(ybc(i+2)-ybc(i-2)) + 16.0_pfdp*(ybc(i+1)-ybc(i-1))/(12.0_pfdp*dx))
          end do
       endif
       deallocate(ybc)
    endif

  end subroutine f1eval1

  subroutine f2eval1(yptr, t, level, levelctx, f2ptr)
    use probin, only: nprob, do_spec
    use mg, only: fill_bc_1d
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, levelctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:), f2(:)
    type(work1),   pointer :: work
    complex(pfdp), pointer :: wk(:)
    real(pfdp), allocatable :: ybc(:)
    real(pfdp) :: dx
    integer :: Nx,i

    y  => array1(yptr)
    f2 => array1(f2ptr)
    Nx = size(y)
    if (do_spec==1) then
       call c_f_pointer(levelctx, work)
       wk => work%wk
       
       wk = y
       call fftw_execute_dft(work%ffft, wk, wk)
       if (nprob .eq. 0) then  ! spectral advection diffusion
          wk = nu * work%lap * wk / size(wk)
       else if (nprob .eq. 1) then  ! backward
          wk = nu * work%lap * wk / size(wk)
       else if (nprob .eq. 2) then  ! trapezoid rule
          wk = 0.5_pfdp*nu * work%lap * wk / size(wk)
       endif
       call fftw_execute_dft(work%ifft, wk, wk)
       
       f2 = real(wk)
    else
       dx = 1.0_pfdp/dble(Nx)
       allocate(ybc(1-2:Nx+2))
       call fill_bc_1d(y,ybc,Nx,2)
       if (spatial_order == 2) then
          do i = 1,Nx
             f2(i)=nu*(-2.0_pfdp*ybc(i) + ybc(i-1)+ybc(i+1))/(dx*dx)
          end do
       else
          do i = 1,Nx
             f2(i)=nu*(-30.0_pfdp*ybc(i) + 16.0_pfdp*(ybc(i-1)+ybc(i+1))-ybc(i-2)-ybc(i+2))/(12.0_pfdp*dx*dx)
!             f2(i)=nu*(-2.0_pfdp*ybc(i) + ybc(i-1)+ybc(i+1))/(dx*dx)
          end do
       endif
       deallocate(ybc)

    endif
  end subroutine f2eval1

  subroutine f2comp1(yptr, t, dt, rhsptr, level, levelctx, f2ptr)
    use mg
    use probin, only: N_Vcycles,do_spec,nprob,Nrelax
    type(c_ptr),    intent(in), value  :: yptr, rhsptr, f2ptr, levelctx
    real(pfdp),     intent(in)         :: t, dt
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:), rhs(:), f2(:)
    real(pfdp)   :: nudt,maxresid,dx
    type(work1),   pointer :: work
    complex(pfdp), pointer :: wk(:)
    integer :: k,N_V,Nx,i,j
    type(ndarray),pointer  :: rhsmg
    type(ndarray),pointer  :: ymg
    real(pfdp), allocatable :: y0(:)
    real(pfdp), allocatable :: rhs0(:)
    real(pfdp), allocatable :: ybc(:)

    Nx = size(y)

    call c_f_pointer(levelctx, work)

    y   => array1(yptr)
    rhs => array1(rhsptr)
    f2  => array1(f2ptr)
    if (nprob .eq. 2) then
       nudt = nu*dt*0.5_pfdp
    else
       nudt = nu*dt
    endif

    if (do_spec .eq. 1) then
       wk  => work%wk
       wk = rhs
       call fftw_execute_dft(work%ffft, wk, wk)
       wk = wk / (1.0_pfdp - nudt*work%lap) / size(wk)
       call fftw_execute_dft(work%ifft, wk, wk)
       
       y  = real(wk)
       f2 = (y - rhs) / dt    
!       call f2eval1(yptr, t, level, levelctx, f2ptr)
   else
      if (nprob .eq. 2) then
         nudt = nu*dt*0.5_pfdp
      else
         nudt = nu*dt
      endif
      call c_f_pointer(yptr,ymg)
      call c_f_pointer(rhsptr,rhsmg)

      allocate(y0(1:Nx))
      dx = 1.0_pfdp/dble(Nx)
      y0=rhs
      if (spatial_order == 4) then
         allocate(rhs0(1:Nx))
         rhs0=rhs + dt*f2
         y = rhs + dt*f2

         !  Subtract L2 from rhs
         allocate(ybc(1-2:Nx+2))
         do j = 1,10
            rhs = rhs0 - dt*f2
            call fill_bc_1d(y,ybc,Nx,2)
            do i = 1,Nx
               rhs(i)=rhs0(i)-nu*dt*(-2.0_pfdp*ybc(i) + ybc(i-1)+ybc(i+1))/(dx*dx)
               rhs(i)=rhs(i)+nu*dt*(-30.0_pfdp*ybc(i) + 16.0_pfdp*(ybc(i-1)+ybc(i+1))-ybc(i-2)-ybc(i+2))/(12.0_pfdp*dx*dx)
            end do
            do k=1,N_Vcycles
               call Vcycle(ymg, t, nudt, rhsmg, level, c_null_ptr, Nrelax,maxresid) 
            end do
            call f2eval1(yptr, t, level, levelctx, f2ptr)
         end do
         deallocate(ybc)
       else
          do k=1,N_Vcycles
             call Vcycle(ymg, t, nudt, rhsmg, level, c_null_ptr, Nrelax,maxresid) 
          end do
          call f2eval1(yptr, t, level, levelctx, f2ptr)
       endif

       deallocate(y0)          

   endif

  end subroutine f2comp1

  !
  ! two dimension
  !

  subroutine f1eval2(yptr, t, level, levelctx, f1ptr)
    type(c_ptr), intent(in), value :: yptr, f1ptr, levelctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level

    type(work2),   pointer :: work
    real(pfdp),    pointer :: y(:,:), f1(:,:)
    complex(pfdp), pointer :: wk(:,:)

    call c_f_pointer(levelctx, work)

    y  => array2(yptr)
    f1 => array2(f1ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = -v * work%ddx * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    f1 = real(wk)
  end subroutine f1eval2

  subroutine f2eval2(yptr, t, level, levelctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, levelctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:,:), f2(:,:)
    type(work2),   pointer :: work
    complex(pfdp), pointer :: wk(:,:)

    call c_f_pointer(levelctx, work)

    y  => array2(yptr)
    f2 => array2(f2ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = nu * work%lap * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    f2 = real(wk)
  end subroutine f2eval2

  subroutine f2comp2(yptr, t, dt, rhsptr, level, levelctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, rhsptr, f2ptr, levelctx
    real(pfdp),     intent(in)         :: t, dt
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:,:), rhs(:,:), f2(:,:)
    type(work2),   pointer :: work
    complex(pfdp), pointer :: wk(:,:)

    call c_f_pointer(levelctx, work)

    y   => array2(yptr)
    rhs => array2(rhsptr)
    f2  => array2(f2ptr)
    wk  => work%wk

    wk = rhs
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = wk / (1.0_pfdp - 0.5_pfdp*nu*dt*work%lap) / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp2

  !
  ! three dimension
  !

  subroutine f1eval3(yptr, t, level, levelctx, f1ptr)
    type(c_ptr), intent(in), value :: yptr, f1ptr, levelctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level

    type(work3),   pointer :: work
    real(pfdp),    pointer :: y(:,:,:), f1(:,:,:)
    complex(pfdp), pointer :: wk(:,:,:)

    call c_f_pointer(levelctx, work)

    y  => array3(yptr)
    f1 => array3(f1ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = -v * work%ddx * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    f1 = real(wk)
  end subroutine f1eval3

  subroutine f2eval3(yptr, t, level, levelctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, f2ptr, levelctx
    real(pfdp),     intent(in)         :: t
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:,:,:), f2(:,:,:)
    type(work3),   pointer :: work
    complex(pfdp), pointer :: wk(:,:,:)

    call c_f_pointer(levelctx, work)

    y  => array3(yptr)
    f2 => array3(f2ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = nu * work%lap * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    f2 = real(wk)
  end subroutine f2eval3

  subroutine f2comp3(yptr, t, dt, rhsptr, level, levelctx, f2ptr)
    type(c_ptr),    intent(in), value  :: yptr, rhsptr, f2ptr, levelctx
    real(pfdp),     intent(in)         :: t, dt
    integer(c_int), intent(in)         :: level

    real(pfdp),    pointer :: y(:,:,:), rhs(:,:,:), f2(:,:,:)
    type(work3),   pointer :: work
    complex(pfdp), pointer :: wk(:,:,:)

    call c_f_pointer(levelctx, work)

    y   => array3(yptr)
    rhs => array3(rhsptr)
    f2  => array3(f2ptr)
    wk  => work%wk

    wk = rhs
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*work%lap) / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp3

  !
  ! work stuff
  !

  subroutine create_work1(levelctx, nvars)
    type(c_ptr), intent(out) :: levelctx
    integer,     intent(in)  :: nvars

    type(work1), pointer :: work
    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx(nvars), dx

    allocate(work)
    levelctx = c_loc(work)

    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx(i) = pi * dble(i-1) / Lx
       else
          kx(i) = pi * dble(-nvars + i - 1) / Lx
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

    work%ddx = (0.0_pfdp, 1.0_pfdp) * kx
     work%lap = -kx**2
    dx = Lx / nvars
    work%lap = (-2 + 2*cos(kx * dx)) / dx**2

    where (kx**2 < 1e-13)
       work%lap = 0.0_pfdp
    end where

  end subroutine create_work1

  subroutine create_work2(levelctx, nvars)
    type(c_ptr), intent(out) :: levelctx
    integer,     intent(in)  :: nvars

    type(work2), pointer :: work
    integer     :: i, j
    type(c_ptr) :: wk
    real(pfdp)  :: kx(nvars)

    allocate(work)
    levelctx = c_loc(work)

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
    allocate(work%ddx(nvars,nvars))
    allocate(work%ddy(nvars,nvars))
    allocate(work%lap(nvars,nvars))

    work%lap = 0

    do j = 1, nvars
       do i = 1, nvars
          work%ddx(i,j) = (0.0_pfdp, 1.0_pfdp) * kx(i)
          work%ddy(i,j) = (0.0_pfdp, 1.0_pfdp) * kx(j)

          if (kx(i)**2 + kx(j)**2 > 0.d0) then
             work%lap = -(kx(i)**2 + kx(j)**2)
          end if
       end do
    end do
  end subroutine create_work2

  subroutine create_work3(levelctx, nvars)
    type(c_ptr), intent(out) :: levelctx
    integer,     intent(in)  :: nvars

    type(work3), pointer :: work
    integer     :: i, j, k
    type(c_ptr) :: wk
    real(pfdp)  :: kx(nvars)

    allocate(work)
    levelctx = c_loc(work)

    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx(i) = two_pi * dble(i-1) / Lx
       else
          kx(i) = two_pi * dble(-nvars + i - 1) / Lx
       end if
    end do

    ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nvars**3, c_size_t))
    call c_f_pointer(wk, work%wk, [nvars,nvars,nvars])

    work%ffft = fftw_plan_dft_3d(nvars, nvars, nvars, &
         work%wk, work%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    work%ifft = fftw_plan_dft_3d(nvars, nvars, nvars, &
         work%wk, work%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
    allocate(work%ddx(nvars,nvars,nvars))
    allocate(work%ddy(nvars,nvars,nvars))
    allocate(work%ddz(nvars,nvars,nvars))
    allocate(work%lap(nvars,nvars,nvars))

    work%lap = 0

    do k = 1, nvars
       do j = 1, nvars
          do i = 1, nvars
             work%ddx(i,j,k) = (0.0_pfdp, 1.0_pfdp) * kx(i)
             work%ddy(i,j,k) = (0.0_pfdp, 1.0_pfdp) * kx(j)
             work%ddy(i,j,k) = (0.0_pfdp, 1.0_pfdp) * kx(k)

             if (kx(i)**2 + kx(j)**2 + kx(k)**2> 0.d0) then
                work%lap = -(kx(i)**2 + kx(j)**2 + kx(k)**2)
             end if
          end do
       end do
    end do
  end subroutine create_work3

end module feval
