! Copyright (c) 2012, Matthew Emmett.  All rights reserved.

module feval
  use iso_c_binding
  use encap
  use fftwpp
  include 'fftw3.f03'

  type :: feval_t
     integer :: n
     real(8), allocatable :: k(:)       ! wave numbers
     real(8) :: nu
     type(c_ptr) :: conv, u1, v1, w1, v2, w2, w3, uu, uv, vv, uw, vw, ww
  end type feval_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine feval_create(n, l, nu, nthreads, cptr)
    use omp_lib
    implicit none

    integer(c_int), intent(in   ), value :: n, nthreads
    real(c_double), intent(in   ), value :: l, nu
    type(c_ptr),    intent(  out)        :: cptr

    type(feval_t), pointer :: ctx
    integer :: k

    allocate(ctx)
    cptr = c_loc(ctx)

    ctx%n = n
    ctx%nu = nu

    allocate(ctx%k(n))
    do k = 1, n
       if (k <= n/2) then
          ctx%k(k) = 6.28318530718d0 * dble(k-1) / l
       else
          ctx%k(k) = 6.28318530718d0 * dble(-n+k-1) / l
       end if
    end do

    ctx%conv = cconv3d_create(n, n, n)

    call omp_set_num_threads(nthreads)
    call set_fftwpp_maxthreads(nthreads)

    ctx%u1 = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%v1 = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%w1 = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%v2 = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%w2 = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%w3 = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%uu = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%uv = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%uw = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%vv = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%vw = fftw_alloc_complex(int(n*n*n, c_size_t))
    ctx%ww = fftw_alloc_complex(int(n*n*n, c_size_t))

  end subroutine feval_create

  subroutine feval_destroy(cptr)
    implicit none
    type(c_ptr), intent(in) :: cptr

    type(feval_t), pointer :: ctx
    call c_f_pointer(cptr, ctx)

    call delete_cconv3d(ctx%conv)

    call fftw_free(ctx%u1)
    call fftw_free(ctx%v1)
    call fftw_free(ctx%w1)
    call fftw_free(ctx%v2)
    call fftw_free(ctx%w2)
    call fftw_free(ctx%w3)
    call fftw_free(ctx%uu)
    call fftw_free(ctx%uv)
    call fftw_free(ctx%uw)
    call fftw_free(ctx%vv)
    call fftw_free(ctx%vw)
    call fftw_free(ctx%ww)

    deallocate(ctx%k)
    deallocate(ctx)

  end subroutine feval_destroy

  subroutine feval_finalize()
    call fftw_cleanup()
  end subroutine feval_finalize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine project(cptr, n1, n2, n3, ustar, u)
    ! Project ustar to divergence free u
    !
    !   u = ustar - grad(inv_lap(div(ustar)))
    !
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: n1, n2, n3
    complex(c_double_complex), intent(in)        :: ustar(n3, n2, n1, 3)
    complex(c_double_complex), intent(out)       :: u(n3, n2, n1, 3)

    type(feval_t), pointer    :: fptr
    integer                   :: i1, i2, i3
    complex(c_double_complex) :: phi
    real(c_double)            :: kk(n1)

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    !$omp parallel do private(i1, i2, i3, phi)
    do i1 = 1,  n1
       do  i2 = 1,  n2
          do  i3 = 1,  n3

             ! phi = div(ustar)
             phi = &
                  kk(i1) * (0.0d0,1.0d0) * ustar(i3, i2, i1, 1) + &
                  kk(i2) * (0.0d0,1.0d0) * ustar(i3, i2, i1, 2) + &
                  kk(i3) * (0.0d0,1.0d0) * ustar(i3, i2, i1, 3)

             ! phi = inv_lap(phi)
             if (i1 > 1 .or. i2 > 1 .or. i3 > 1) then
                phi = -phi / (kk(i1)**2 + kk(i2)**2 + kk(i3)**2)
             else
                phi = 0.0d0
             end if

             ! u = ustar - grad(phi)
             u(i3, i2, i1, 1) = ustar(i3, i2, i1, 1) - kk(i1) * (0.0d0,1.0d0) * phi
             u(i3, i2, i1, 2) = ustar(i3, i2, i1, 2) - kk(i2) * (0.0d0,1.0d0) * phi
             u(i3, i2, i1, 3) = ustar(i3, i2, i1, 3) - kk(i3) * (0.0d0,1.0d0) * phi

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine project

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copy_for_convolve(y, ctx)
    type(carray4), intent(in   ) :: y
    type(feval_t), intent(in   ) :: ctx

    complex(c_double_complex), dimension(:,:,:), pointer :: &
         u1, v1, w1, v2, w2, w3, uu, uv, uw, vv, vw, ww

    call c_f_pointer(ctx%u1, u1, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%v1, v1, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%w1, w1, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%v2, v2, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%w2, w2, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%w3, w3, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%uu, uu, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%uv, uv, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%uw, uw, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%vv, vv, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%vw, vw, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%ww, ww, [ ctx%n, ctx%n, ctx%n ])

    !$omp parallel workshare
    u1 = y%array(:,:,:,1)
    v1 = y%array(:,:,:,2)
    w1 = y%array(:,:,:,3)
    v2 = y%array(:,:,:,2)
    w2 = y%array(:,:,:,3)
    w3 = y%array(:,:,:,3)
    uu = y%array(:,:,:,1)
    uv = y%array(:,:,:,1)
    vv = y%array(:,:,:,2)
    uw = y%array(:,:,:,1)
    vw = y%array(:,:,:,2)
    ww = y%array(:,:,:,3)
    !$omp end parallel workshare
  end subroutine copy_for_convolve

  subroutine eval_f1(yptr, t, level, ctxptr, f1ptr)
    type(c_ptr), intent(in   ), value :: yptr, ctxptr, f1ptr
    real(pfdp),  intent(in   )        :: t
    integer,     intent(in   )        :: level

    complex(c_double_complex), dimension(:,:,:), pointer :: uu, uv, uw, vv, vw, ww
    type(carray4)                                        :: u

    type(feval_t), pointer :: ctx
    type(carray4), pointer :: y, f1

    call c_f_pointer(yptr, y)
    call c_f_pointer(f1ptr, f1)
    call c_f_pointer(ctxptr, ctx)

    call carray4_create(u, y%shape)
    call project(ctxptr, ctx%n, ctx%n, ctx%n, y%array, u%array)

    call copy_for_convolve(u, ctx)

    call cconv3d_convolve(ctx%conv, ctx%uu, ctx%u1)
    call cconv3d_convolve(ctx%conv, ctx%uv, ctx%v1)
    call cconv3d_convolve(ctx%conv, ctx%uw, ctx%w1)
    call cconv3d_convolve(ctx%conv, ctx%vv, ctx%v2)
    call cconv3d_convolve(ctx%conv, ctx%vw, ctx%w2)
    call cconv3d_convolve(ctx%conv, ctx%ww, ctx%w3)

    call c_f_pointer(ctx%uu, uu, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%uv, uv, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%uw, uw, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%vv, vv, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%vw, vw, [ ctx%n, ctx%n, ctx%n ])
    call c_f_pointer(ctx%ww, ww, [ ctx%n, ctx%n, ctx%n ])

    call f1eval(ctxptr, ctx%n, ctx%n, ctx%n, uu, uv, vv, uw, vw, ww, f1%array)

    deallocate(u%array)

  end subroutine eval_f1

  subroutine f1eval(cptr, n1, n2, n3,  uu, uv, vv, uw, vw, ww , f1)
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: n1, n2, n3
    complex(c_double_complex), intent(in)        :: &
         uu(n3, n2, n1), uv(n3, n2, n1), vv(n3, n2, n1), &
         uw(n3, n2, n1), vw(n3, n2, n1), ww(n3, n2, n1)
    complex(c_double_complex), intent(out)       :: f1(n3, n2, n1, 3)

    type(feval_t), pointer :: fptr
    integer                :: i1, i2, i3
    real(c_double)         :: kk(n1)

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    !$omp parallel do private(i1, i2, i3)
    do i1 = 1,  n1
       do  i2 = 1,  n2
          do  i3 = 1,  n3

             f1(i3, i2, i1, 1) = -( &
                  kk(i1) * (0.0d0,1.0d0) * uu(i3, i2, i1) + &
                  kk(i2) * (0.0d0,1.0d0) * uv(i3, i2, i1) + &
                  kk(i3) * (0.0d0,1.0d0) * uw(i3, i2, i1)  )
             f1(i3, i2, i1, 2) = -( &
                  kk(i1) * (0.0d0,1.0d0) * uv(i3, i2, i1) + &
                  kk(i2) * (0.0d0,1.0d0) * vv(i3, i2, i1) + &
                  kk(i3) * (0.0d0,1.0d0) * vw(i3, i2, i1)  )
             f1(i3, i2, i1, 3) = -( &
                  kk(i1) * (0.0d0,1.0d0) * uw(i3, i2, i1) + &
                  kk(i2) * (0.0d0,1.0d0) * vw(i3, i2, i1) + &
                  kk(i3) * (0.0d0,1.0d0) * ww(i3, i2, i1)  )

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine f1eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(yptr, t, level, ctxptr, f2ptr)
    type(c_ptr),   intent(in   ), value :: yptr, ctxptr, f2ptr
    real(pfdp),    intent(in   )        :: t
    integer,       intent(in   )        :: level

    type(feval_t), pointer :: ctx
    type(carray4), pointer :: y, f2
    call c_f_pointer(ctxptr, ctx)
    call c_f_pointer(yptr, y)
    call c_f_pointer(f2ptr, f2)

    call f2eval(ctxptr, ctx%n, ctx%n, ctx%n, y%array, ctx%nu, f2%array)
  end subroutine eval_f2

  subroutine f2eval(cptr, n1, n2, n3, ustar, nu, f2)
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: n1, n2, n3
    real(c_double),            intent(in), value :: nu
    complex(c_double_complex), intent(in)        :: ustar(n3, n2, n1, 3)
    complex(c_double_complex), intent(out)       :: f2(n3, n2, n1, 3)

    type(feval_t), pointer :: fptr
    integer                :: i1, i2, i3
    real(c_double)         :: kk(n1), lap

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    !$omp parallel do private(i1, i2, i3, lap)
    do i1 = 1,  n1
       do  i2 = 1,  n2
          do  i3 = 1,  n3

             lap = - (kk(i1)**2 + kk(i2)**2 + kk(i3)**2)
             f2(i3, i2, i1, 1) = nu * lap * ustar(i3, i2, i1, 1)
             f2(i3, i2, i1, 2) = nu * lap * ustar(i3, i2, i1, 2)
             f2(i3, i2, i1, 3) = nu * lap * ustar(i3, i2, i1, 3)

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_f2(yptr, t, dt, rhsptr, level, ctxptr, f2ptr)
    type(c_ptr), intent(in   ), value :: yptr, rhsptr, ctxptr, f2ptr
    real(pfdp),  intent(in   )        :: t, dt
    integer,     intent(in   )        :: level

    type(carray4), pointer :: y, rhs, f2
    type(feval_t), pointer :: ctx
    call c_f_pointer(ctxptr, ctx)
    call c_f_pointer(yptr, y)
    call c_f_pointer(f2ptr, f2)
    call c_f_pointer(rhsptr, rhs)

    call f2solv(ctxptr, ctx%n, ctx%n, ctx%n, rhs%array, y%array, ctx%nu, dt, f2%array)
  end subroutine comp_f2

  subroutine f2solv(cptr, n1, n2, n3, rhs, ustar, nu, dt, f2)
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: n1, n2, n3
    real(c_double),            intent(in), value :: nu, dt
    complex(c_double_complex), intent(in)        :: rhs(n3, n2, n1, 3)
    complex(c_double_complex), intent(out)       :: ustar(n3, n2, n1, 3), f2(n3, n2, n1, 3)

    type(feval_t), pointer :: fptr
    integer                :: i1, i2, i3
    real(c_double)         :: kk(n1), lap

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    !$omp parallel do private(i1, i2, i3, lap)
    do i1 = 1,  n1
       do  i2 = 1,  n2
          do  i3 = 1,  n3

             lap = -(kk(i1)**2 + kk(i2)**2 + kk(i3)**2)
             ustar(i3, i2, i1, 1) = rhs(i3, i2, i1, 1) / (1.0d0 - nu*dt*lap)
             ustar(i3, i2, i1, 2) = rhs(i3, i2, i1, 2) / (1.0d0 - nu*dt*lap)
             ustar(i3, i2, i1, 3) = rhs(i3, i2, i1, 3) / (1.0d0 - nu*dt*lap)
             f2(i3, i2, i1, 1)    = nu * lap * ustar(i3, i2, i1, 1)
             f2(i3, i2, i1, 2)    = nu * lap * ustar(i3, i2, i1, 2)
             f2(i3, i2, i1, 3)    = nu * lap * ustar(i3, i2, i1, 3)

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine f2solv

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module feval
