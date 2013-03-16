module encap_array1d
  use iso_c_binding
  use pfasst
  implicit none

  type :: array1d
     integer :: nvars
     real(pfdp), pointer :: array(:)
  end type array1d

contains

  function array(ptr) result(r)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp), pointer :: r(:)

    type(array1d), pointer :: q
    call c_f_pointer(ptr, q)

    r => q%array
  end function array

  subroutine array1d_encap_create(encap)
    type(pf_encap_t), intent(out) :: encap

    encap%create  => encap_create
    encap%destroy => encap_destroy
    encap%setval  => encap_setval
    encap%copy    => encap_copy
    encap%pack    => encap_pack
    encap%unpack  => encap_unpack
    encap%axpy    => encap_axpy
  end subroutine array1d_encap_create

  ! Allocate/create solution (spatial data set) for the given level.
  !
  ! This is called for each SDC node.
  subroutine encap_create(sol, level, feval, nvars, shape, levelctx, encapctx)
    type(c_ptr),       intent(inout)     :: sol
    integer,           intent(in)        :: level, nvars, shape(:)
    logical,           intent(in)        :: feval
    type(c_ptr),       intent(in), value :: levelctx, encapctx

    type(array1d), pointer :: q

    allocate(q)
    allocate(q%array(nvars))

    sol = c_loc(q)
  end subroutine encap_create

  ! Deallocate/destroy solution.
  subroutine encap_destroy(ptr)
    type(c_ptr), intent(in), value :: ptr

    type(array1d), pointer :: q
    call c_f_pointer(ptr, q)

    deallocate(q%array)
    deallocate(q)
  end subroutine encap_destroy

  ! Set solution value.
  subroutine encap_setval(ptr, val, flags)
    type(c_ptr), intent(in), value    :: ptr
    real(pfdp),  intent(in)           :: val
    integer,     intent(in), optional :: flags

    real(pfdp), pointer :: q(:)

    q => array(ptr)
    q = val
  end subroutine encap_setval

  ! Copy solution value.
  subroutine encap_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in), value    :: dstptr, srcptr
    integer,     intent(in), optional :: flags

    real(pfdp), pointer :: dst(:), src(:)

    dst => array(dstptr)
    src => array(srcptr)

    dst = src
  end subroutine encap_copy

  ! Pack solution q into a flat array.
  subroutine encap_pack(z, ptr)
    type(c_ptr), intent(in), value  :: ptr
    real(pfdp),  intent(out)        :: z(:)

    real(pfdp), pointer :: q(:)
    q => array(ptr)

    z = q
  end subroutine encap_pack

  ! Unpack solution from a flat array.
  subroutine encap_unpack(ptr, z)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp),  intent(in)        :: z(:)

    real(pfdp), pointer :: q(:)
    q => array(ptr)

    q = z
  end subroutine encap_unpack

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine encap_axpy(yptr, a, xptr, flags)
    type(c_ptr), intent(in), value    :: xptr, yptr
    real(pfdp),  intent(in)           :: a
    integer,     intent(in), optional :: flags

    real(pfdp), pointer :: x(:), y(:)
    x => array(xptr)
    y => array(yptr)

    y = a * x + y
  end subroutine encap_axpy

end module encap_array1d
