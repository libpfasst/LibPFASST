!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module encap
  use iso_c_binding
  implicit none
  ! Data precision
  integer, parameter :: pfdp = kind(0.0d0)
  real(pfdp), parameter :: ZERO = 0.0_pfdp
  real(pfdp), parameter :: ONE = 1.0_pfdp
  real(pfdp), parameter :: TWO = 2.0_pfdp
  real(pfdp), parameter :: HALF = 0.5_pfdp


  ! Data encapsulation: simple 1d array
  type :: pf_encap_t
     real(pfdp), pointer :: array(:)
  end type pf_encap_t

  interface create
     module procedure encap_create
  end interface create

  interface destroy
     module procedure encap_destroy
  end interface destroy

  interface setval
     module procedure encap_setval
  end interface setval

  interface copy
     module procedure encap_copy
  end interface copy

  interface axpy
     module procedure encap_axpy
  end interface axpy
 
  interface pack
     module procedure encap_pack
  end interface pack

  interface unpack
     module procedure encap_unpack
  end interface unpack
 
contains

  ! Allocate/create solution (spatial data set) for the given level.
  !
  ! This is called for each SDC node.
  subroutine encap_create(sol, level, feval, nvars, shape, ctx)
    type(pf_encap_t),  intent(inout) :: sol
    integer,           intent(in)    :: level, nvars, shape(:)
    logical,           intent(in)    :: feval
    type(c_ptr),       intent(in)    :: ctx

    allocate(sol%array(nvars))
  end subroutine encap_create

  ! Deallocate/destroy solution.
  subroutine encap_destroy(sol)
    type(pf_encap_t), intent(inout) :: sol

    deallocate(sol%array)
  end subroutine encap_destroy

  ! Set solution value.
  subroutine encap_setval(sol, val)
    type(pf_encap_t), intent(inout) :: sol
    real(pfdp),       intent(in)    :: val

    sol%array = val
  end subroutine encap_setval

  ! Copy solution value.
  subroutine encap_copy(dst, src)
    type(pf_encap_t), intent(inout) :: dst
    type(pf_encap_t), intent(in)    :: src

    dst%array = src%array
  end subroutine encap_copy

  ! Pack solution q into a flat array.
  subroutine encap_pack(z, q)
    type(pf_encap_t), intent(in)  :: q
    real(pfdp),       intent(out) :: z(:)

    z = q%array
  end subroutine encap_pack

  ! Unpack solution from a flat array.
  subroutine encap_unpack(q, z)
    type(pf_encap_t), intent(inout) :: q
    real(pfdp),       intent(in)    :: z(:)

    q%array = z
  end subroutine encap_unpack

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine encap_axpy(y, a, x)
    real(pfdp),       intent(in)    :: a
    type(pf_encap_t), intent(in)    :: x
    type(pf_encap_t), intent(inout) :: y

    y%array = a * x%array + y%array
  end subroutine encap_axpy

end module encap
