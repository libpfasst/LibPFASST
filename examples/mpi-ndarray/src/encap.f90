!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module encap
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  type :: ndarray
     integer             :: dim
     integer,    pointer :: shape(:)
     real(pfdp), pointer :: array(:)
  end type ndarray

contains

  subroutine ndarray_create_simple(q, shape)
    type(ndarray), intent(inout) :: q
    integer,       intent(in)    :: shape(:)
    allocate(q%shape(size(shape)))
    allocate(q%array(product(shape)))
    q%dim   = size(shape)
    q%shape = shape
  end subroutine ndarray_create_simple

  ! Allocate/create solution (spatial data set).
  subroutine ndarray_create(solptr, level, kind, nvars, shape, levelctx, encapctx)
    type(c_ptr),       intent(inout)     :: solptr
    integer,           intent(in)        :: level, nvars, shape(:)
    integer,           intent(in)        :: kind
    type(c_ptr),       intent(in), value :: levelctx, encapctx

    type(ndarray),     pointer :: sol

    allocate(sol)
    allocate(sol%shape(size(shape)))
    allocate(sol%array(product(shape)))

    sol%dim   = size(shape)
    sol%shape = shape
    solptr    = c_loc(sol)
  end subroutine ndarray_create

  ! Deallocate/destroy solution.
  subroutine ndarray_destroy(solptr)
    type(c_ptr), intent(in), value :: solptr

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    deallocate(sol%array)
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

    sol%array = val
  end subroutine ndarray_setval

  ! Copy solution value.
  subroutine ndarray_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in), value    :: dstptr, srcptr
    integer,     intent(in), optional :: flags

    type(ndarray), pointer :: dst, src
    call c_f_pointer(dstptr, dst)
    call c_f_pointer(srcptr, src)

    dst%array = src%array
  end subroutine ndarray_copy

  ! Pack solution q into a flat array.
  subroutine ndarray_pack(z, ptr)
    type(c_ptr), intent(in), value  :: ptr
    real(pfdp),  intent(out)        :: z(:)
    real(pfdp), pointer :: q(:)
    q => array1(ptr)
    z = q
  end subroutine ndarray_pack

  ! Unpack solution from a flat array.
  subroutine ndarray_unpack(ptr, z)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp),  intent(in)        :: z(:)
    real(pfdp), pointer :: q(:)
    q => array1(ptr)
    q = z
  end subroutine ndarray_unpack

  ! Compute norm of solution
  function ndarray_norm(ptr) result (norm)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp) :: norm
    real(pfdp), pointer :: q(:)
    q => array1(ptr)
    norm = maxval(abs(q))
  end function ndarray_norm

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine ndarray_saxpy(yptr, a, xptr, flags)
    type(c_ptr), intent(in), value    :: yptr, xptr
    real(pfdp),  intent(in)           :: a
    integer,     intent(in), optional :: flags

    type(ndarray), pointer :: y, x
    call c_f_pointer(yptr, y)
    call c_f_pointer(xptr, x)

    y%array = a * x%array + y%array
  end subroutine ndarray_saxpy

  ! Helpers
  function array1(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:)

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    if (sol%dim == 1) then
       call c_f_pointer(c_loc(sol%array(1)), r, sol%shape)
    else
       stop
    end if
  end function array1

  function array2(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:,:)

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    if (sol%dim == 2) then
       call c_f_pointer(c_loc(sol%array(1)), r, sol%shape)
    else
       stop
    end if
  end function array2

  function array3(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:,:,:)

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    if (sol%dim == 3) then
       call c_f_pointer(c_loc(sol%array(1)), r, sol%shape)
    else
       stop
    end if
  end function array3

  function array4(solptr) result(r)
    type(c_ptr), intent(in), value :: solptr
    real(pfdp), pointer :: r(:,:,:,:)

    type(ndarray), pointer :: sol
    call c_f_pointer(solptr, sol)

    if (sol%dim == 4) then
       call c_f_pointer(c_loc(sol%array(1)), r, sol%shape)
    else
       stop
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

end module encap
