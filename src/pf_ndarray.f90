!!  N-dimensional array encapsulation.
!
! This file is part of LIBPFASST.
!

!> N-dimensional array encapsulation.
!!
!! When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'shape' attribute to create a new array with that
!! shape.  Thus, the 'shape' attributes of the PFASST levels should be
!! set appropriately.  For example, before calling pf_pfasst_run we can
!! set the shape of the coarsest level by doing:
!!
!!   allocate(pf%levels(1)%shape(2))
!!   pf%levels(1)%shape = [ 3, 10 ]
!!
!! The helper routines array1, array2, array3, etc can be used to
!! extract pointers to the encapsulated array without
!! performing any copies.
!!
module pf_mod_ndarray
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_utils  
  implicit none

  !>  Type to create and destroy N-dimenstional arrays
  type, extends(pf_factory_t) :: ndarray_factory
   contains
     procedure :: create_single  => ndarray_create_single
     procedure :: create_array  => ndarray_create_array
     procedure :: destroy_single => ndarray_destroy_single
     procedure :: destroy_array => ndarray_destroy_array
  end type ndarray_factory

  !>  N-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: ndarray
     integer             :: dim
     integer,    allocatable :: shape(:)
     real(pfdp), allocatable :: flatarray(:)
   contains
     procedure :: setval => ndarray_setval
     procedure :: copy => ndarray_copy
     procedure :: norm => ndarray_norm
     procedure :: pack => ndarray_pack
     procedure :: unpack => ndarray_unpack
     procedure :: axpy => ndarray_axpy
     procedure :: eprint => ndarray_eprint
  end type ndarray

  !> Interfaces to output routines in pf_numpy.c
  interface
     !>  Subroutine to make a directory for output
     subroutine ndarray_mkdir(dname, dlen) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname
       integer,    intent(in   ), value :: dlen
     end subroutine ndarray_mkdir
     !>  Subroutine to write an the array to a file
     subroutine ndarray_dump_numpy(dname, fname, endian, dim, mpibuflen, shape, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname, fname, endian(5)
       integer,    intent(in   ), value :: dim, mpibuflen
       integer,    intent(in   )        :: shape(dim)
       real(c_double),    intent(in   )        :: array(mpibuflen)
     end subroutine ndarray_dump_numpy
  end interface

contains
  function cast_as_ndarray(encap_polymorph) result(ndarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(ndarray), pointer :: ndarray_obj
    
    select type(encap_polymorph)
    type is (ndarray)
       ndarray_obj => encap_polymorph
    end select
  end function cast_as_ndarray

  !>  Subroutine to allocate the array and set the size parameters
  subroutine ndarray_build(q, shape)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape(:)

    select type (q)
    class is (ndarray)
       allocate(q%shape(size(shape)))
       allocate(q%flatarray(product(shape)))
       q%dim   = size(shape)
       q%shape = shape
    end select
  end subroutine ndarray_build

  !> Subroutine to  create a single array
  subroutine ndarray_create_single(this, x, level, shape)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, shape(:)
    integer :: i
    allocate(ndarray::x)
    call ndarray_build(x, shape)
  end subroutine ndarray_create_single

  !> Subroutine to create an array of arrays
  subroutine ndarray_create_array(this, x, n, level,  shape)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, shape(:)
    integer :: i
    allocate(ndarray::x(n))
    do i = 1, n
       call ndarray_build(x(i), shape)
    end do
  end subroutine ndarray_create_array

  !>  Subroutine to destroy array
  subroutine ndarray_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(ndarray), pointer :: ndarray_obj

    ndarray_obj => cast_as_ndarray(encap)

    deallocate(ndarray_obj%shape)
    deallocate(ndarray_obj%flatarray)

    nullify(ndarray_obj)

  end subroutine ndarray_destroy

  !> Subroutine to destroy an single array
  subroutine ndarray_destroy_single(this, x)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x

    select type (x)
    class is (ndarray)
       deallocate(x%shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine ndarray_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine ndarray_destroy_array(this, x)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout),allocatable :: x(:)
    integer                                            :: i

    select type(x)
    class is (ndarray)
       do i = 1,size(x)
          deallocate(x(i)%shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
  end subroutine ndarray_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalare  value.
  subroutine ndarray_setval(this, val, flags)
    class(ndarray), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%flatarray = val
  end subroutine ndarray_setval

  !> Subroutine to copy an array
  subroutine ndarray_copy(this, src, flags)
    class(ndarray),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (ndarray)
       this%flatarray = src%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine ndarray_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine ndarray_pack(this, z, flags)
    class(ndarray), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    z = this%flatarray
  end subroutine ndarray_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine ndarray_unpack(this, z, flags)
    class(ndarray), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags
    this%flatarray = z
  end subroutine ndarray_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function ndarray_norm(this, flags) result (norm)
    class(ndarray), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%flatarray))
  end function ndarray_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine ndarray_axpy(this, a, x, flags)
    class(ndarray),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags


    select type(x)
    type is (ndarray)
       this%flatarray = a * x%flatarray + this%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine ndarray_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine ndarray_eprint(this,flags)
    class(ndarray), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
    if (product(this%shape) < 10) then
       print *, this%flatarray
    else
       print *, this%flatarray(1:10)
    endif
  end subroutine ndarray_eprint



  !>  Helper function to return the array part
  function get_array1d(x,flags) result(r)
    class(pf_encap_t), target, intent(in) :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:)
    select type (x)
    type is (ndarray)
       r => x%flatarray
    end select
  end function get_array1d
  

  function get_array2d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:)

    select type (x)
    type is (ndarray)
       r(1:x%shape(1),1:x%shape(2)) => x%flatarray
    end select
  end function get_array2d
  

  function get_array3d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (ndarray)
       r(1:x%shape(1),1:x%shape(2),1:x%shape(3)) => x%flatarray
    end select
  end function get_array3d
  






end module pf_mod_ndarray
