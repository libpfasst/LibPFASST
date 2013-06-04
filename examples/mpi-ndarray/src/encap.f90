module encap
  use pf_mod_dtype, only: pf_encap_t, pf_context_t, pfdp
  implicit none

  type, extends(pf_encap_t) :: array1d
     integer :: nvars
     real(pfdp), pointer :: array(:)
   contains
     procedure :: create  => encap_create
     procedure :: destroy => encap_destroy
     procedure :: setval  => encap_setval
     procedure :: norm    => encap_norm
     procedure :: pack    => encap_pack
     procedure :: unpack  => encap_unpack
     procedure :: copy    => encap_copy
     procedure :: axpy    => encap_axpy
  end type array1d

contains

  function get_array(q) result(r)
    class(pf_encap_t), intent(in), target :: q
    real(pfdp),        pointer    :: r(:)
    select type(q)
    class is (array1d)
       r => q%array
    end select
  end function get_array

  ! Instantiate a new encapsulation class.
  subroutine encap_new(r)
    class(pf_encap_t), intent(out), allocatable :: r
    allocate(array1d::r)
  end subroutine encap_new

  ! Allocate/create solution (spatial data set) for the given level.
  !
  ! This is called for each SDC node.
  subroutine encap_create(sol, level, kind, nvars, shape, ctx)
    class(array1d),      intent(inout) :: sol
    integer,             intent(in)    :: kind, level, nvars, shape(:)
    class(pf_context_t), intent(in)    :: ctx

    allocate(sol%array(nvars))
    sol%nvars = nvars
  end subroutine encap_create

  ! Deallocate/destroy solution.
  subroutine encap_destroy(sol)
    class(array1d), intent(inout) :: sol
    deallocate(sol%array)
  end subroutine encap_destroy

  ! Set solution value.
  subroutine encap_setval(sol, val, flags)
    class(array1d), intent(inout) :: sol
    real(pfdp),     intent(in)    :: val
    integer,        intent(in), optional :: flags
    sol%array = val
  end subroutine encap_setval

  ! Copy solution value.
  subroutine encap_copy(dst, src, flags)
    class(array1d),    intent(inout) :: dst
    class(pf_encap_t), intent(in)    :: src
    integer,           intent(in), optional :: flags
    real(pfdp), pointer :: src_array(:)
    src_array => get_array(src)
    dst%array = src_array
  end subroutine encap_copy

  ! Compute norm of solution
  function encap_norm(sol) result (norm)
    class(array1d), intent(in) :: sol
    real(pfdp) :: norm
    norm = maxval(abs(sol%array))
  end function encap_norm

  ! Pack solution q into a flat array.
  subroutine encap_pack(sol, z)
    class(array1d), intent(in)  :: sol
    real(pfdp),     intent(out) :: z(:)
    z = sol%array
  end subroutine encap_pack

  ! Unpack solution from a flat array.
  subroutine encap_unpack(sol, z)
    class(array1d), intent(inout) :: sol
    real(pfdp),     intent(in)    :: z(:)
    sol%array = z
  end subroutine encap_unpack

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine encap_axpy(y, a, x, flags)
    class(array1d),    intent(inout) :: y
    real(pfdp),        intent(in)    :: a
    class(pf_encap_t), intent(in)    :: x
    integer,           intent(in), optional :: flags
    real(pfdp), pointer :: x_array(:)
    x_array => get_array(x)
    y%array = a * x_array + y%array
  end subroutine encap_axpy

end module encap
