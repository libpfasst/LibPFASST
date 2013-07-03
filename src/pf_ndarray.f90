!
! Copyright (C) 2013 Matthew Emmett.
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
! N-dimensional array encapsulation.
!
! When a new solution is created by a PFASST level, this encapsulation
! uses the levels 'shape' attribute to create a new array with that
! shape.  Thus, the 'shape' attributes of the PFASST levels should be
! set appropriately.  For example, before calling pf_pfasst_run we can
! set the shape of the coarsest level by doing:
!
!   allocate(pf%levels(1)%shape(2))
!   pf%levels(1)%shape = [ 3, 10 ]
!
! The helper routines array1, array2, array3, and array4 can be used
! to extract pointers to the encapsulated array from a C pointer
! without performing any copies.
!

module pf_mod_ndarray
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  type :: ndarray
     integer             :: dim
     integer,    pointer :: shape(:)
     real(pfdp), pointer :: flatarray(:)
     type(c_ptr)         :: aptr
  end type ndarray

contains

  subroutine ndarray_create_simple(q, shape)
    type(ndarray), intent(inout) :: q
    integer,       intent(in)    :: shape(:)
    allocate(q%shape(size(shape)))
    allocate(q%flatarray(product(shape)))
    q%dim   = size(shape)
    q%shape = shape
    q%aptr  = c_loc(q%flatarray(1))
  end subroutine ndarray_create_simple

  ! Allocate/create solution (spatial data set).
  subroutine ndarray_create(solptr, level, kind, nvars, shape, levelctx, encapctx)
    type(c_ptr),       intent(inout)     :: solptr
    integer,           intent(in)        :: level, nvars, shape(:)
    integer,           intent(in)        :: kind
    type(c_ptr),       intent(in), value :: levelctx, encapctx

    type(ndarray),     pointer :: sol

    allocate(sol)
    call ndarray_create_simple(sol, shape)
    solptr = c_loc(sol)
  end subroutine ndarray_create

  ! Deallocate/destroy solution.
  subroutine ndarray_destroy(solptr)
    type(c_ptr), intent(in), value :: solptr

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    deallocate(sol%flatarray)
    deallocate(sol%shape)
    deallocate(sol)
  end subroutine ndarray_destroy

  ! Set solution value.
  subroutine ndarray_setval(solptr, val, flags)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp),  intent(in)        :: val
    integer,     intent(in), optional :: flags

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    sol%flatarray = val
  end subroutine ndarray_setval

  ! Copy solution value.
  subroutine ndarray_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in), value    :: dstptr, srcptr
    integer,     intent(in), optional :: flags

    type(ndarray), pointer :: dst, src
    call c_f_pointer(dstptr, dst)
    call c_f_pointer(srcptr, src)

    dst%flatarray = src%flatarray
  end subroutine ndarray_copy

  ! Pack solution q into a flat array.
  subroutine ndarray_pack(z, ptr)
    type(c_ptr), intent(in), value  :: ptr
    real(pfdp),  intent(out)        :: z(:)
    type(ndarray), pointer :: y
    call c_f_pointer(ptr, y)
    z = y%flatarray
  end subroutine ndarray_pack

  ! Unpack solution from a flat array.
  subroutine ndarray_unpack(ptr, z)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp),  intent(in)        :: z(:)
    type(ndarray), pointer :: y
    call c_f_pointer(ptr, y)
    y%flatarray = z
  end subroutine ndarray_unpack

  ! Compute norm of solution
  function ndarray_norm(ptr) result (norm)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp) :: norm
    type(ndarray), pointer :: y
    call c_f_pointer(ptr, y)
    norm = maxval(abs(y%flatarray))
  end function ndarray_norm

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine ndarray_saxpy(yptr, a, xptr, flags)
    type(c_ptr), intent(in), value    :: yptr, xptr
    real(pfdp),  intent(in)           :: a
    integer,     intent(in), optional :: flags

    type(ndarray), pointer :: y, x
    call c_f_pointer(yptr, y)
    call c_f_pointer(xptr, x)

    y%flatarray = a * x%flatarray + y%flatarray
  end subroutine ndarray_saxpy

  ! Helpers
  function dims(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    integer :: r

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    r = sol%dim
  end function dims

  function array1(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:)

    integer                :: shp(1)
    type(ndarray), pointer :: sol

    call c_f_pointer(solptr, sol)
    if (sol%dim == 1) then
       shp = sol%shape
       call c_f_pointer(sol%aptr, r, shp)
    else
       stop "ERROR: array1 dimension mismatch."
    end if
  end function array1

  function array2(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:,:)

    integer                :: shp(2)
    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    if (sol%dim == 2) then
       shp = sol%shape
       call c_f_pointer(sol%aptr, r, shp)
    else
       stop "ERROR: array2 dimension mismatch."
    end if
  end function array2

  function array3(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:,:,:)

    integer                :: shp(3)
    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    if (sol%dim == 3) then
       shp = sol%shape
       call c_f_pointer(sol%aptr, r, shp)
    else
       stop "ERROR: array3 dimension mismatch."
    end if
  end function array3

  function array4(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:,:,:,:)

    integer                :: shp(4)
    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    if (sol%dim == 4) then
       shp = sol%shape
       call c_f_pointer(sol%aptr, r, shp)
    else
       stop "ERROR: array4 dimension mismatch."
    end if
  end function array4

  subroutine ndarray_encap_create(encap)
    type(pf_encap_t), intent(out) :: encap

    encap%create  => ndarray_create
    encap%destroy => ndarray_destroy
    encap%setval  => ndarray_setval
    encap%copy    => ndarray_copy
    encap%norm    => ndarray_norm
    encap%pack    => ndarray_pack
    encap%unpack  => ndarray_unpack
    encap%axpy    => ndarray_saxpy
  end subroutine ndarray_encap_create

end module pf_mod_ndarray
