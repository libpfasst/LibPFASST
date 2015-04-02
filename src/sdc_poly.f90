!
! Copyright (C) 2012, 2013 Matthew Emmett.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

!
! Polynomial manipulation routines.
!
! A polynomial p
!
!   p(x) = a_n x^n + ... + a_2 x^2 + a_1 x + a_0
!
! is stored as a Fortran array p(0:n) according to
!
!   p = [ a_0, a_1, ..., a_n ].
!

module sdc_mod_poly
  use iso_c_binding
  implicit none

  integer,  parameter :: qp = c_long_double
  integer,  parameter :: dp = c_double
  real(qp), parameter :: eps = 1.0e-23_qp

  private :: qsort_partition

  interface poly_eval
     module procedure poly_eval
     module procedure poly_eval_complex
  end interface

contains

  !
  ! Evaluate polynomial.
  !
  real(qp) function poly_eval(p, n, x) result(v) bind(c)
    integer(c_int), intent(in), value :: n
    real(qp),       intent(in)        :: p(0:n), x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function

  complex(qp) function poly_eval_complex(p, n, x) result(v)
    integer(c_int), intent(in), value :: n
    real(qp),       intent(in)        :: p(0:n)
    complex(qp),    intent(in)        :: x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function


  !
  ! Differentiate polynomial (in place)
  !
  subroutine poly_diff(p, n) bind(c)
    integer(c_int), intent(in),   value :: n
    real(qp),       intent(inout) :: p(0:n)

    integer  :: j
    real(qp) :: pp(0:n)

    pp = 0.0_qp

    do j = 1, n
       pp(j-1) = j * p(j)
    end do

    p = pp
  end subroutine poly_diff


  !
  ! Integrate polynomial (in place)
  !
  subroutine poly_int(p, n) bind(c)
    integer(c_int), intent(in),   value :: n
    real(qp),       intent(inout) :: p(0:n)

    integer  :: j
    real(qp) :: pp(0:n)

    pp = 0.0_qp

    do j = 0, n-1
       pp(j+1) = p(j) / (j+1)
    end do

    p = pp
  end subroutine poly_int


  !
  ! Compute Legendre polynomial coefficients using Bonnet's recursion
  ! formula.
  !
  subroutine poly_legendre(p, n) bind(c)
    integer(c_int), intent(in), value :: n
    real(qp),       intent(out)       :: p(0:n)

    real(qp), dimension(0:n) :: p0, p1, p2
    integer :: j, m

    if (n == 0) then
       p = [ 1.0_qp ]
       return
    end if

    if (n == 1) then
       p = [ 0.0_qp, 1.0_qp ]
       return
    end if

    p0 = 0.0_qp; p1 = 0.0_qp; p2 = 0.0_qp

    p0(0) = 1.0_qp
    p1(1) = 1.0_qp

    ! (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
    do m = 1, n-1
       do j = 1, n
          p2(j) = ( (2*m + 1) * p1(j-1) - m * p0(j) ) / (m + 1)
       end do
       p2(0) = - m * p0(0) / (m + 1)

       p0 = p1
       p1 = p2
    end do

    p = p2
  end subroutine poly_legendre


  !
  ! Compute polynomial roots using the Durand-Kerner algorithm.
  !
  ! The roots are assumed to be real.
  !
  subroutine poly_roots(roots, p0, n) bind(c)
    integer(c_int),  intent(in), value   :: n
    real(qp),        intent(out)  :: roots(n)
    real(qp),        intent(in)   :: p0(0:n)

    integer     :: i, j, k
    complex(qp) :: num, den, z0(n), z1(n)
    real(qp)    :: p(0:n)

    p = p0 / p0(n)

    ! initial guess
    do i = 1, n
       z0(i) = (0.4_qp, 0.9_qp)**i
    end do

    ! durand-kerner-weierstrass iterations
    z1 = z0
    do k = 1, 100
       do i = 1, n

          ! evaluate poly at z0(i)
          num = poly_eval(p, n, z0(i))

          ! evaluate denominator
          den = 1.0_qp
          do j = 1, n
             if (j == i) cycle
             den = den * (z0(i) - z0(j))
          end do

          ! update
          z0(i) = z0(i) - num / den
       end do

       ! converged?
       if (sum(abs(z0 - z1)) < eps) exit

       z1 = z0
    end do

    roots = real(z0)
    where (abs(roots) < eps) roots = 0.0_qp
    call qsort(roots)

  end subroutine poly_roots


  !
  ! Sort (inplace) using the quick sort algorithm.
  !
  ! Adapted from http://www.fortran.com/qsort_c.f95.
  !
  recursive subroutine qsort(a)
    real(qp), intent(inout) :: a(:)
    integer :: iq

    if (size(a) > 1) then
       call qsort_partition(a, iq)
       call qsort(a(:iq-1))
       call qsort(a(iq:))
    end if
  end subroutine qsort

  subroutine qsort_partition(a, marker)
    real(qp), intent(inout) :: a(:)
    integer,  intent(out)   :: marker

    integer  :: i, j
    real(qp) :: temp, x

    x = a(1)
    i = 0
    j = size(a) + 1

    do
       j = j-1
       do
          if (a(j) <= x) exit
          j = j-1
       end do

       i = i+1
       do
          if (a(i) >= x) exit
          i = i+1
       end do

       if (i < j) then
          temp = a(i)
          a(i) = a(j)
          a(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
  end subroutine qsort_partition

end module sdc_mod_poly
