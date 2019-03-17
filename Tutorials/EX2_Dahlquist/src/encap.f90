!
! This file is part of LIBPFASST.
!

!> Simple scalar encapsulation example
!!

module encap
  use pf_mod_dtype
  implicit none

  !>  Type to create and destroy the arrays
  type, extends(pf_factory_t) :: scalar_factory
   contains
     procedure :: create_single  => scalar_create_single
     procedure :: create_array  => scalar_create_array
     procedure :: destroy_single => scalar_destroy_single
     procedure :: destroy_array => scalar_destroy_array
  end type scalar_factory

  !>  Type to extend the abstract encap and set procedure pointers
  type, extends(pf_encap_t) :: scalar_encap
     real(pfdp) :: y   !  The scalar value
   contains
     procedure :: setval => scalar_setval
     procedure :: copy => scalar_copy
     procedure :: norm => scalar_norm
     procedure :: pack => scalar_pack
     procedure :: unpack => scalar_unpack
     procedure :: axpy => scalar_axpy
     procedure :: eprint => scalar_eprint
  end type scalar_encap

contains
  !>  Subroutine to allocate the array and set the size parameters
  subroutine scalar_create_single(this, x, level, shape)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, shape(:)
    integer :: i

    allocate(scalar_encap::x)
  end subroutine scalar_create_single

  !> Subroutine to create an array of arrays
  subroutine scalar_create_array(this, x, n, level,  shape)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, shape(:)
    integer :: i
    allocate(scalar_encap::x(n))
  end subroutine scalar_create_array

  !>  Subroutine to destroy array
  subroutine scalar_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(scalar_encap), pointer :: scalar_obj

    scalar_obj => cast_as_scalar(encap)

    nullify(scalar_obj)

  end subroutine scalar_destroy

  !> Subroutine to destroy an single array
  subroutine scalar_destroy_single(this, x)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x

    deallocate(x)
  end subroutine scalar_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine scalar_destroy_array(this, x)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer                                            :: i

    deallocate(x)

  end subroutine scalar_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalare  value.
  subroutine scalar_setval(this, val, flags)
    class(scalar_encap), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%y = val
  end subroutine scalar_setval

  !> Subroutine to copy an array
  subroutine scalar_copy(this, src, flags)
    class(scalar_encap),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (scalar_encap)
    this%y = src%y
    class default
       stop "TYPE ERROR"
    end select

  end subroutine scalar_copy

  !> Subroutine to pack into a flat array for sending
  subroutine scalar_pack(this, z, flags)
    class(scalar_encap), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    z = this%y
  end subroutine scalar_pack

  !> Subroutine to unpack  after receiving
  subroutine scalar_unpack(this, z, flags)
    class(scalar_encap), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags
    this%y = z(1)
  end subroutine scalar_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function scalar_norm(this, flags) result (norm)
    class(scalar_encap), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = abs(this%y)
  end function scalar_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine scalar_axpy(this, a, x, flags)
    class(scalar_encap),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags

    select type(x)
    type is (scalar_encap)
       this%y = a * x%y + this%y
    class default
       stop "TYPE ERROR"
    end select
  end subroutine scalar_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine scalar_eprint(this,flags)
    class(scalar_encap), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
    print *, this%y
  end subroutine scalar_eprint


  !  Helper function to cast an abstract encap to the scalar_encap
  function cast_as_scalar(encap_polymorph) result(scalar_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(scalar_encap), pointer :: scalar_obj
    
    select type(encap_polymorph)
    type is (scalar_encap)
       scalar_obj => encap_polymorph
    end select
  end function cast_as_scalar



end module encap
