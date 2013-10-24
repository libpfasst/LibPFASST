! Copyright (c) 2012, Matthew Emmett.  All rights reserved.

module feval
  use iso_c_binding

  type :: feval_t
     integer :: n
     real(8), allocatable :: k(:)       ! wave numbers
  end type feval_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine create(n, l, nthreads, cptr) bind(c, name='create')
    use omp_lib
    implicit none

    integer(c_int), intent(in), value :: n, nthreads
    real(c_double), intent(in), value :: l
    type(c_ptr),    intent(out)       :: cptr

    type(feval_t), pointer :: fptr
    integer              :: k


    ! setup

    allocate(fptr)
    cptr = c_loc(fptr)

    ! init n and k (wave numbers)

    fptr%n = n

    allocate(fptr%k(n))
    do k = 1, n
       if (k <= n/2) then
          fptr%k(k) = 6.28318530718d0 * dble(k-1) / l
       else
          fptr%k(k) = 6.28318530718d0 * dble(-n+k-1) / l
       end if
    end do

    ! set threads

    call omp_set_num_threads(nthreads)

  end subroutine create

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine divergence(cptr, <{n}>, u, div) bind(c, name='divergence')
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: <{n}>
    complex(c_double_complex), intent(in)        :: u({{n;.}})
    complex(c_double_complex), intent(out)       :: div({{n}})

    type(feval_t), pointer      :: fptr
    integer                   :: <{i}>
    real(c_double)            :: kk(n1)

    ! setup

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    !$omp parallel do private(<{i}>)
    do multi(<{i}>; <{n}>)
       ! phi = div(ustar)
       div({{i}}) = &
            [{ kk(i1) * (0.0d0,1.0d0) * u({{i;1}}); + &
               kk(i2) * (0.0d0,1.0d0) * u({{i;2}}); + &
               kk(i3) * (0.0d0,1.0d0) * u({{i;3}}) }]
    end do
    !$omp end parallel do

  end subroutine divergence

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine project(cptr, <{n}>, ustar, u) bind(c, name='project')
    ! Project ustar to divergence free u
    !
    !   u = ustar - grad(inv_lap(div(ustar)))
    !
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: <{n}>
    complex(c_double_complex), intent(in)        :: ustar({{n;.}})
    complex(c_double_complex), intent(out)       :: u({{n;.}})

    type(feval_t), pointer      :: fptr
    integer                   :: <{i}>
    complex(c_double_complex) :: phi
    real(c_double)            :: kk(n1)

    ! setup

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    !$omp parallel do private(<{i}>, phi)
    do multi(<{i}>; <{n}>)
       ! phi = div(ustar)
       phi = &
            [{ kk(i1) * (0.0d0,1.0d0) * ustar({{i;1}}); + &
               kk(i2) * (0.0d0,1.0d0) * ustar({{i;2}}); + &
               kk(i3) * (0.0d0,1.0d0) * ustar({{i;3}}) }]

       ! phi = inv_lap(phi)
       if ([{i1 > 1; .or. i2 > 1; .or. i3 > 1}]) then
          phi = -phi / ([{kk(i1)**2; + kk(i2)**2; + kk(i3)**2}])
       else
          phi = 0.0d0
       end if

       ! u = ustar - grad(phi)
       [{
       u({{i;1}}) = ustar({{i;1}}) - kk(i1) * (0.0d0,1.0d0) * phi;
       u({{i;2}}) = ustar({{i;2}}) - kk(i2) * (0.0d0,1.0d0) * phi;
       u({{i;3}}) = ustar({{i;3}}) - kk(i3) * (0.0d0,1.0d0) * phi
       }]
    end do
    !$omp end parallel do

  end subroutine project

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine f1eval(cptr, <{n}>, [{ uu;, uv, vv;, uw, vw, ww }], f1) bind(c, name='f1eval')
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: <{n}>
    complex(c_double_complex), intent(in)        :: &
      [{ uu({{n}});, uv({{n}}), vv({{n}});, uw({{n}}), vw({{n}}), ww({{n}}) }]
    complex(c_double_complex), intent(out)       :: f1({{n;.}})

    type(feval_t), pointer      :: fptr
    integer                   :: <{i}>
    real(c_double)            :: kk(n1)

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    ! compute divergence

    !$omp parallel do private(<{i}>)
    do multi(<{i}>; <{n}>)
       [{
       f1({{i;1}}) = -( &
            [{ kk(i1) * (0.0d0,1.0d0) * uu({{i}}); + &
               kk(i2) * (0.0d0,1.0d0) * uv({{i}}); + &
               kk(i3) * (0.0d0,1.0d0) * uw({{i}}) }] );
       f1({{i;2}}) = -( &
            [{ kk(i1) * (0.0d0,1.0d0) * uv({{i}}); + &
               kk(i2) * (0.0d0,1.0d0) * vv({{i}}); + &
               kk(i3) * (0.0d0,1.0d0) * vw({{i}}) }] );
       f1({{i;3}}) = -( &
            [{ kk(i1) * (0.0d0,1.0d0) * uw({{i}}); + &
               kk(i2) * (0.0d0,1.0d0) * vw({{i}}); + &
               kk(i3) * (0.0d0,1.0d0) * ww({{i}}) }] )
       }]
    end do
    !$omp end parallel do

  end subroutine f1eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine f2eval(cptr, <{n}>, ustar, nu, f2) bind(c, name='f2eval')
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: <{n}>
    real(c_double),            intent(in), value :: nu
    complex(c_double_complex), intent(in)        :: ustar({{n;.}})
    complex(c_double_complex), intent(out)       :: f2({{n;.}})

    type(feval_t), pointer :: fptr
    integer              :: <{i}>, c
    real(c_double)       :: kk(n1), lap

    ! setup

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    ! evaluate

    do c = 1, ({1;2;3})
       !$omp parallel do private(<{i}>, lap)
       do multi(<{i}>; <{n}>)
          lap         = - ([{kk(i1)**2; + kk(i2)**2; + kk(i3)**2}])
          f2({{i;c}}) = nu * lap * ustar({{i;c}})
       end do
       !$omp end parallel do
    end do

  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine f2solv(cptr, <{n}>, rhs, ustar, nu, dt, f2) bind(c, name='f2solv')
    implicit none

    type(c_ptr),               intent(in), value :: cptr
    integer(c_int),            intent(in), value :: <{n}>
    real(c_double),            intent(in), value :: nu, dt
    complex(c_double_complex), intent(in)        :: rhs({{n;.}})
    complex(c_double_complex), intent(out)       :: ustar({{n;.}}), f2({{n;.}})

    type(feval_t), pointer :: fptr
    integer              :: <{i}>, c
    real(c_double)       :: kk(n1), lap

    ! setup

    call c_f_pointer(cptr, fptr)
    kk = fptr%k

    ! solve

    do c = 1, ({1;2;3})
       !$omp parallel do private(<{i}>, lap)
       do multi(<{i}>; <{n}>)
          lap            = -([{kk(i1)**2; + kk(i2)**2; + kk(i3)**2}])
          ustar({{i;c}}) = rhs({{i;c}}) / (1.0d0 - nu*dt*lap)
          f2({{i;c}})    = nu * lap * ustar({{i;c}})
       end do
       !$omp end parallel do
    end do

  end subroutine f2solv

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module feval
