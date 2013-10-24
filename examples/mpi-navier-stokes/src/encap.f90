!
! Copyright (c) 2013, Matthew Emmett.  All rights reserved.
!

module encap
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  ! Data encapsulation: simple multicomponent 3d array
  type :: carray4
     integer :: nvars, shape(4)
     complex(pfdp), pointer :: array(:,:,:,:)
  end type carray4

contains


  subroutine carray4_create(sol, shape)
    type(carray4), intent(inout) :: sol
    integer,       intent(in   ) :: shape(:)

    ! if (2*product(shape) .ne. nvars) then
    !    stop "NVARS/SHAPE MISMATCH"
    ! end if

    sol%shape = shape
    sol%nvars = 2*product(shape)
    allocate(sol%array(sol%shape(1), sol%shape(2), sol%shape(3), sol%shape(4)))
  end subroutine carray4_create


  ! Allocate/create solution (spatial data set) for the given level.
  !
  ! This is called for each SDC node.
  subroutine encap_create(solptr, level, kind, nvars, shape, levelctx, encapctx)
    type(c_ptr), intent(inout)        :: solptr
    integer,     intent(in   )        :: level, nvars, shape(:)
    integer,     intent(in   )        :: kind
    type(c_ptr), intent(in   ), value :: levelctx, encapctx

    type(carray4), pointer :: sol

    allocate(sol)
    call carray4_create(sol, shape)
    if (nvars /= sol%nvars) then
       stop "NVARS/SHAPE MISMATCH (encap.f90)"
    end if
    solptr = c_loc(sol)
  end subroutine encap_create

  ! Deallocate/destroy solution.
  subroutine encap_destroy(solptr)
    type(c_ptr), intent(in   ), value :: solptr

    type(carray4), pointer :: sol
    call c_f_pointer(solptr, sol)

    deallocate(sol%array)
    deallocate(sol)
  end subroutine encap_destroy

  ! Set solution value.
  subroutine encap_setval(solptr, val, flags)
    type(c_ptr), intent(in   ), value    :: solptr
    real(pfdp),  intent(in   )           :: val
    integer,     intent(in   ), optional :: flags

    type(carray4), pointer :: sol
    call c_f_pointer(solptr, sol)

    !$omp parallel workshare
    sol%array = val
    !$omp end parallel workshare
  end subroutine encap_setval

  ! Copy solptrution value.
  subroutine encap_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in   ), value    :: dstptr, srcptr
    integer,     intent(in   ), optional :: flags

    type(carray4), pointer :: src, dst
    call c_f_pointer(dstptr, dst)
    call c_f_pointer(srcptr, src)

    !$omp parallel workshare
    dst%array = src%array
    !$omp end parallel workshare
  end subroutine encap_copy

  ! Pack solution into a flat array.
  subroutine encap_pack(q, solptr)
    type(c_ptr), intent(in   ), value :: solptr
    real(pfdp),  intent(  out)        :: q(:)

    type(c_ptr)            :: qptr
    type(carray4), pointer :: sol
    real(pfdp), pointer    :: qz(:)

    call c_f_pointer(solptr, sol)
    qptr = c_loc(sol%array(1,1,1,1))
    call c_f_pointer(qptr, qz, [ sol%nvars ])

    !$omp parallel workshare
    q = qz
    !$omp end parallel workshare
  end subroutine encap_pack

  ! Unpack solution from a flat array.
  subroutine encap_unpack(solptr, q)
    type(c_ptr), intent(in   ), value :: solptr
    real(pfdp),  intent(  out)        :: q(:)

    type(c_ptr)            :: qptr
    type(carray4), pointer :: sol
    real(pfdp), pointer    :: qz(:)

    call c_f_pointer(solptr, sol)
    qptr = c_loc(sol%array(1,1,1,1))
    call c_f_pointer(qptr, qz, [ sol%nvars ])

    !$omp parallel workshare
    qz = q
    !$omp end parallel workshare
  end subroutine encap_unpack

  ! Compute norm of solution
  function encap_norm(solptr) result (norm)
    type(c_ptr), intent(in   ), value :: solptr
    real(pfdp) :: norm

    type(carray4), pointer :: sol
    call c_f_pointer(solptr, sol)

    norm = maxval(abs(sol%array))
  end function encap_norm

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine encap_axpy(yptr, a, xptr, flags)
    type(c_ptr), intent(in   ), value    :: yptr, xptr
    real(pfdp),  intent(in   )           :: a
    integer,     intent(in   ), optional :: flags

    type(carray4), pointer :: x, y
    call c_f_pointer(xptr, x)
    call c_f_pointer(yptr, y)

    !$omp parallel workshare
    y%array = a * x%array + y%array
    !$omp end parallel workshare
  end subroutine encap_axpy

  subroutine carray4_encap_create(encap)
    type(pf_encap_t), intent(  out) :: encap

    encap%create  => encap_create
    encap%destroy => encap_destroy
    encap%setval  => encap_setval
    encap%copy    => encap_copy
    encap%norm    => encap_norm
    encap%pack    => encap_pack
    encap%unpack  => encap_unpack
    encap%axpy    => encap_axpy
  end subroutine carray4_encap_create

end module encap
