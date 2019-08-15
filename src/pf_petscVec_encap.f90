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
module pf_mod_pf_petscVec
#include <petsc/finclude/petscvec.h>
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_utils  
  use petscvec
  implicit none

  !>  Type to create and destroy N-dimenstional arrays
  type, extends(pf_factory_t) :: pf_petscVec_factory
   contains
     procedure :: create_single  => pf_petscVec_create_single
     procedure :: create_array  => pf_petscVec_create_array
     procedure :: destroy_single => pf_petscVec_destroy_single
     procedure :: destroy_array => pf_petscVec_destroy_array
  end type pf_petscVec_factory

  !>  N-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: pf_petscVec
     integer             :: dim
     integer,    allocatable :: shape(:)
     real(pfdp), allocatable :: flatarray(:)
     type(tVec) ::  foo
   contains
     procedure :: setval => pf_petscVec_setval
     procedure :: copy => pf_petscVec_copy
     procedure :: norm => pf_petscVec_norm
     procedure :: pack => pf_petscVec_pack
     procedure :: unpack => pf_petscVec_unpack
     procedure :: axpy => pf_petscVec_axpy
     procedure :: eprint => pf_petscVec_eprint
  end type pf_petscVec

  !> Interfaces to output routines in pf_numpy.c
  interface
     !>  Subroutine to write an the array to a file
     subroutine pf_petscVec_dump_numpy(dname, fname, endian, dim, mpibuflen, shape, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname, fname, endian(5)
       integer,    intent(in   ), value :: dim, mpibuflen
       integer,    intent(in   )        :: shape(dim)
       real(c_double),    intent(in   )        :: array(mpibuflen)
     end subroutine pf_petscVec_dump_numpy
  end interface

contains
  function cast_as_pf_petscVec(encap_polymorph) result(pf_petscVec_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pf_petscVec), pointer :: pf_petscVec_obj
    
    select type(encap_polymorph)
    type is (pf_petscVec)
       pf_petscVec_obj => encap_polymorph
    end select
  end function cast_as_pf_petscVec

  !>  Subroutine to allocate the array and set the size parameters
  subroutine pf_petscVec_build(q, shape)
    use pf_mod_comm_mpi
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape(:)

    integer ierr,nn
    select type (q)
    class is (pf_petscVec)
       print *,'creating petscVec',   PETSC_COMM_WORLD    
       call VecCreate(PETSC_COMM_WORLD,q%foo,ierr)

       call VecSetSizes(q%foo,PETSC_DECIDE,shape(1),ierr)
       call VecSetFromOptions(q%foo,ierr)
       allocate(q%shape(size(shape)))
       allocate(q%flatarray(product(shape)))
       q%dim   = size(shape)
       q%shape = shape
    end select
  end subroutine pf_petscVec_build

  !> Subroutine to  create a single array
  subroutine pf_petscVec_create_single(this, x, level, shape)
    class(pf_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, shape(:)
    integer :: i
    allocate(pf_petscVec::x)
    call pf_petscVec_build(x, shape)
  end subroutine pf_petscVec_create_single

  !> Subroutine to create an array of arrays
  subroutine pf_petscVec_create_array(this, x, n, level,  shape)
    class(pf_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, shape(:)
    integer :: i
    allocate(pf_petscVec::x(n))
    do i = 1, n
       call pf_petscVec_build(x(i), shape)
    end do
  end subroutine pf_petscVec_create_array

  !>  Subroutine to destroy array
  subroutine pf_petscVec_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_petscVec), pointer :: pf_petscVec_obj
    integer :: ierr

    pf_petscVec_obj => cast_as_pf_petscVec(encap)
    print *,'destroying petscVec'
    call VecDestroy(pf_petscVec_obj%foo,ierr)
    deallocate(pf_petscVec_obj%shape)
    deallocate(pf_petscVec_obj%flatarray)

    nullify(pf_petscVec_obj)

  end subroutine pf_petscVec_destroy

  !> Subroutine to destroy an single array
  subroutine pf_petscVec_destroy_single(this, x)
    class(pf_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer :: ierr
    select type (x)
    class is (pf_petscVec)
       print *,'destroying petscVec'
       call VecDestroy(x%foo,ierr)
       
       deallocate(x%shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine pf_petscVec_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine pf_petscVec_destroy_array(this, x)
    class(pf_petscVec_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout),allocatable :: x(:)
    integer                                            :: i
    integer :: ierr
    select type(x)
    class is (pf_petscVec)
       do i = 1,size(x)
          print *,'destroying  petscVec'
          call VecDestroy(x(i)%foo,ierr)
          deallocate(x(i)%shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
  end subroutine pf_petscVec_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalare  value.
  subroutine pf_petscVec_setval(this, val, flags)
    class(pf_petscVec), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%flatarray = val
  end subroutine pf_petscVec_setval

  !> Subroutine to copy an array
  subroutine pf_petscVec_copy(this, src, flags)
    class(pf_petscVec),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (pf_petscVec)
       this%flatarray = src%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine pf_petscVec_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine pf_petscVec_pack(this, z, flags)
    class(pf_petscVec), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    z = this%flatarray
  end subroutine pf_petscVec_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine pf_petscVec_unpack(this, z, flags)
    class(pf_petscVec), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags
    this%flatarray = z
  end subroutine pf_petscVec_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function pf_petscVec_norm(this, flags) result (norm)
    class(pf_petscVec), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%flatarray))
  end function pf_petscVec_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine pf_petscVec_axpy(this, a, x, flags)
    class(pf_petscVec),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags


    select type(x)
    type is (pf_petscVec)
       this%flatarray = a * x%flatarray + this%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine pf_petscVec_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine pf_petscVec_eprint(this,flags)
    class(pf_petscVec), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
    if (product(this%shape) < 10) then
       print *, this%flatarray
    else
       print *, this%flatarray(1:10)
    endif
  end subroutine pf_petscVec_eprint



  !>  Helper function to return the array part
  function get_array1d(x,flags) result(r)
    class(pf_encap_t), target, intent(in) :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:)
    select type (x)
    type is (pf_petscVec)
       r => x%flatarray
    end select
  end function get_array1d
  

  function get_array2d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:)

    select type (x)
    type is (pf_petscVec)
       r(1:x%shape(1),1:x%shape(2)) => x%flatarray
    end select
  end function get_array2d
  

  function get_array3d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (pf_petscVec)
       r(1:x%shape(1),1:x%shape(2),1:x%shape(3)) => x%flatarray
    end select
  end function get_array3d
  






end module pf_mod_pf_petscVec
