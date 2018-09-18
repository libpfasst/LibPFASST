!!  N-dimensional complex array encapsulation.
!
! This file is part of LIBPFASST.
!

!> System of complex N-dimensional arrays encapsulation.
!!
!! When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'arr_shape' attribute to create a new multi-component array with that
!! shape.  Thus, the 'arr_shape' attributes of the PFASST levels should be
!! set appropriately.  The last component of arr_shape is the number of components in the system
!!
!! For example, before calling pf_pfasst_run we can
!! set the arr_shape of the coarsest level by doing:
!!
!!   allocate(pf%levels(1)%arr_shape(3))
!!   pf%levels(1)%arr_shape = [ nx, ny, 3 ]
!!
!! Which would imply that a 3 component system of two-dimensional solutions.
!!
!! The helper routines array1, array2, array3, etc can be used to
!! extract pointers to a component of  encapsulated system
!! performing any copies.
!!
module pf_mod_zndsysarray
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  !>  Type to create and destroy the arrays
  type, extends(pf_factory_t) :: zndsysarray_factory
   contains
     procedure :: create_single  => zndsysarray_create_single
     procedure :: create_array  => zndsysarray_create_array
     procedure :: destroy_single => zndsysarray_destroy_single
     procedure :: destroy_array => zndsysarray_destroy_array
  end type zndsysarray_factory

  !>  Type to extend the abstract encap and set procedure pointers
  type, extends(pf_encap_t) :: zndsysarray
     integer             :: dim    !  The spatial dimension of each component in system
     integer             :: ncomp  !  The number of components in the system
     integer             :: ndof   !  The number of variables in each component
     integer,    allocatable :: arr_shape(:)
     complex(pfdp), allocatable :: flatarray(:)
   contains
     procedure :: setval => zndsysarray_setval
     procedure :: copy => zndsysarray_copy
     procedure :: norm => zndsysarray_norm
     procedure :: pack => zndsysarray_pack
     procedure :: unpack => zndsysarray_unpack
     procedure :: axpy => zndsysarray_axpy
     procedure :: eprint => zndsysarray_eprint
  end type zndsysarray

  !> Interfaces to output routines in pf_numpy.c
  interface
     !>  Subroutine to make a directory for output
     subroutine zndsysarray_mkdir(dname, dlen) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname
       integer,    intent(in   ), value :: dlen
     end subroutine zndsysarray_mkdir
     !>  Subroutine to write an the array to a file
     subroutine zndsysarray_dump_numpy(dname, fname, endian, dim, mpibuflen, arr_shape, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname, fname, endian(5)
       integer,    intent(in   ), value :: dim, mpibuflen
       integer,    intent(in   )        :: arr_shape(dim)
       real(c_double),    intent(in   )        :: array(mpibuflen)
     end subroutine zndsysarray_dump_numpy
  end interface

contains
  !>  Subroutine to allocate the array and set the size parameters
  subroutine zndsysarray_build(q, arr_shape)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: arr_shape(:)

    select type (q)
    class is (zndsysarray)
       allocate(q%arr_shape(size(arr_shape)))
       q%dim   = size(arr_shape)-1
       q%ncomp = arr_shape(q%dim+1)
       q%ndof = product(arr_shape(1:q%dim))
       q%arr_shape = arr_shape

       allocate(q%flatarray(product(arr_shape)))
    end select
  end subroutine zndsysarray_build

  !> Subroutine to  create a single array
  subroutine zndsysarray_create_single(this, x, level, shape)
    class(zndsysarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, shape(:)
    integer :: i
    allocate(zndsysarray::x)
    call zndsysarray_build(x, shape)
  end subroutine zndsysarray_create_single

  !> Subroutine to create an array of arrays
  subroutine zndsysarray_create_array(this, x, n, level,  shape)
    class(zndsysarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, shape(:)
    integer :: i
    allocate(zndsysarray::x(n))
    do i = 1, n
       call zndsysarray_build(x(i), shape)
    end do
  end subroutine zndsysarray_create_array

  !>  Subroutine to destroy array
  subroutine zndsysarray_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(zndsysarray), pointer :: zndsysarray_obj

    zndsysarray_obj => cast_as_zndsysarray(encap)

    deallocate(zndsysarray_obj%arr_shape)
    deallocate(zndsysarray_obj%flatarray)

    nullify(zndsysarray_obj)

  end subroutine zndsysarray_destroy

  !> Subroutine to destroy an single array
  subroutine zndsysarray_destroy_single(this, x, level, shape)
    class(zndsysarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, shape(:)

    select type (x)
    class is (zndsysarray)
       deallocate(x%arr_shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine zndsysarray_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine zndsysarray_destroy_array(this, x, n, level,  shape)
    class(zndsysarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, shape(:)
    integer                                            :: i

    select type(x)
    class is (zndsysarray)
       do i = 1,n
          deallocate(x(i)%arr_shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
  end subroutine zndsysarray_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalare  value.
  subroutine zndsysarray_setval(this, val, flags)
    class(zndsysarray), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%flatarray = val
  end subroutine zndsysarray_setval

  !> Subroutine to copy an array
  subroutine zndsysarray_copy(this, src, flags)
    class(zndsysarray),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (zndsysarray)
       this%flatarray = src%flatarray
    class default
       stop "TYPE ERROR"
    end select
  end subroutine zndsysarray_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine zndsysarray_pack(this, z, flags)
    class(zndsysarray), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: ntot
    ntot = this%ndof*this%ncomp
    
    z(1:ntot) = real(this%flatarray,pfdp)
    z(ntot+1:2*ntot) = aimag(this%flatarray)
  end subroutine zndsysarray_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine zndsysarray_unpack(this, z, flags)
    class(zndsysarray), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags

    integer :: ntot    
    ntot = this%ndof*this%ncomp
    this%flatarray = z(1:ntot)
    this%flatarray = this%flatarray + cmplx(0.0,1.0,pfdp)*z(ntot+1:2*ntot)
    this%flatarray =  cmplx(z(1:ntot),z(ntot+1:2*ntot))        
  end subroutine zndsysarray_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function zndsysarray_norm(this, flags) result (norm)
    class(zndsysarray), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%flatarray))
  end function zndsysarray_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine zndsysarray_axpy(this, a, x, flags)
    class(zndsysarray),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags

    select type(x)
    type is (zndsysarray)
       this%flatarray = a * x%flatarray + this%flatarray
    class default
       stop "TYPE ERROR"
    end select
  end subroutine zndsysarray_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine zndsysarray_eprint(this,flags)
    class(zndsysarray), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
        print *, this%flatarray(1:10)
    !print *, this%flatarray
  end subroutine zndsysarray_eprint


  function cast_as_zndsysarray(encap_polymorph) result(zndsysarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(zndsysarray), pointer :: zndsysarray_obj
    
    select type(encap_polymorph)
    type is (zndsysarray)
       zndsysarray_obj => encap_polymorph
    end select
  end function cast_as_zndsysarray

  !>  Helper function to return the array part
  function get_array1d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:)
    select type (x)
    type is (zndsysarray)
       r => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
    end select
  end function get_array1d
  

  function get_array2d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n    
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:)


    select type (x)
    type is (zndsysarray)
       r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
    end select
  end function get_array2d
  

  function get_array3d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (zndsysarray)
       r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
    end select
  end function get_array3d
  






end module pf_mod_zndsysarray
