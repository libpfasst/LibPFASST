!!  N-dimensional array encapsulation using petsc.
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
module pf_mod_my_petscVec
#include <petsc/finclude/petscvec.h>
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_utils  
  use petscvec
  implicit none

  !>  Type to create and destroy N-dimenstional arrays
  type, extends(pf_factory_t) :: my_petscVec_factory
   contains
     procedure :: create_single  => my_petscVec_create_single
     procedure :: create_array  => my_petscVec_create_array
     procedure :: destroy_single => my_petscVec_destroy_single
     procedure :: destroy_array => my_petscVec_destroy_array
  end type my_petscVec_factory

  !>  N-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: my_petscVec
     integer             :: dim
     integer,    allocatable :: shape(:)
     real(pfdp), allocatable :: flatarray(:)
     type(tVec),pointer ::  foo
   contains
     procedure :: setval => my_petscVec_setval
     procedure :: copy => my_petscVec_copy
     procedure :: norm => my_petscVec_norm
     procedure :: pack => my_petscVec_pack
     procedure :: unpack => my_petscVec_unpack
     procedure :: axpy => my_petscVec_axpy
     procedure :: eprint => my_petscVec_eprint
  end type my_petscVec

  !> Interfaces to output routines in pf_numpy.c
  interface
     !>  Subroutine to write an the array to a file
     subroutine my_petscVec_dump_numpy(dname, fname, endian, dim, mpibuflen, shape, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname, fname, endian(5)
       integer,    intent(in   ), value :: dim, mpibuflen
       integer,    intent(in   )        :: shape(dim)
       real(c_double),    intent(in   )        :: array(mpibuflen)
     end subroutine my_petscVec_dump_numpy
  end interface

contains
  function cast_as_my_petscVec(encap_polymorph) result(my_petscVec_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(my_petscVec), pointer :: my_petscVec_obj
    
    select type(encap_polymorph)
    type is (my_petscVec)
       my_petscVec_obj => encap_polymorph
    end select
  end function cast_as_my_petscVec

  !>  Subroutine to allocate the array and set the size parameters
  subroutine my_petscVec_build(q, shape)
    use pf_mod_comm_mpi
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape(:)

    integer ierr,nn
    select type (q)
    class is (my_petscVec)
       !       call VecCreateSeq(PETSC_COMM_SELF,shape(1),q%foo,ierr);CHKERRA(ierr)
       call VecCreateSeq(0,shape(1),q%foo,ierr)

       call VecSetSizes(q%foo,PETSC_DECIDE,shape(1),ierr)
       call VecSetFromOptions(q%foo,ierr)

       call VecGetSize(q%foo,nn,ierr)
       print *,nn
       allocate(q%shape(size(shape)))
       allocate(q%flatarray(product(shape)))
       q%dim   = size(shape)
       q%shape = shape
    end select
  end subroutine my_petscVec_build

  !> Subroutine to  create a single array
  subroutine my_petscVec_create_single(this, x, level, shape)
    class(my_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, shape(:)
    integer :: i
    allocate(my_petscVec::x)
    call my_petscVec_build(x, shape)
  end subroutine my_petscVec_create_single

  !> Subroutine to create an array of arrays
  subroutine my_petscVec_create_array(this, x, n, level,  shape)
    class(my_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, shape(:)
    integer :: i
    allocate(my_petscVec::x(n))
    do i = 1, n
       call my_petscVec_build(x(i), shape)
    end do
  end subroutine my_petscVec_create_array

  !>  Subroutine to destroy array
  subroutine my_petscVec_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(my_petscVec), pointer :: my_petscVec_obj

    my_petscVec_obj => cast_as_my_petscVec(encap)

    deallocate(my_petscVec_obj%shape)
    deallocate(my_petscVec_obj%flatarray)

    nullify(my_petscVec_obj)

  end subroutine my_petscVec_destroy

  !> Subroutine to destroy an single array
  subroutine my_petscVec_destroy_single(this, x)
    class(my_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x

    select type (x)
    class is (my_petscVec)
       deallocate(x%shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine my_petscVec_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine my_petscVec_destroy_array(this, x)
    class(my_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout),allocatable :: x(:)
    integer                                            :: i

    select type(x)
    class is (my_petscVec)
       do i = 1,size(x)
          deallocate(x(i)%shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
  end subroutine my_petscVec_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalare  value.
  subroutine my_petscVec_setval(this, val, flags)
    class(my_petscVec), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%flatarray = val
  end subroutine my_petscVec_setval

  !> Subroutine to copy an array
  subroutine my_petscVec_copy(this, src, flags)
    class(my_petscVec),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (my_petscVec)
       this%flatarray = src%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine my_petscVec_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine my_petscVec_pack(this, z, flags)
    class(my_petscVec), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    z = this%flatarray
  end subroutine my_petscVec_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine my_petscVec_unpack(this, z, flags)
    class(my_petscVec), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags
    this%flatarray = z
  end subroutine my_petscVec_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function my_petscVec_norm(this, flags) result (norm)
    class(my_petscVec), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%flatarray))
  end function my_petscVec_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine my_petscVec_axpy(this, a, x, flags)
    class(my_petscVec),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags


    select type(x)
    type is (my_petscVec)
       this%flatarray = a * x%flatarray + this%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine my_petscVec_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine my_petscVec_eprint(this,flags)
    class(my_petscVec), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
    if (product(this%shape) < 10) then
       print *, this%flatarray
    else
       print *, this%flatarray(1:10)
    endif
  end subroutine my_petscVec_eprint



  !>  Helper function to return the array part
  function get_array1d(x,flags) result(r)
    class(pf_encap_t), target, intent(in) :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:)
    select type (x)
    type is (my_petscVec)
       r => x%flatarray
    end select
  end function get_array1d
  

  function get_array2d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:)

    select type (x)
    type is (my_petscVec)
       r(1:x%shape(1),1:x%shape(2)) => x%flatarray
    end select
  end function get_array2d
  

  function get_array3d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (my_petscVec)
       r(1:x%shape(1),1:x%shape(2),1:x%shape(3)) => x%flatarray
    end select
  end function get_array3d
  






end module pf_mod_my_petscVec
