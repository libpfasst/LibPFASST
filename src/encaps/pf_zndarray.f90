!!  N-dimensional complex array encapsulation
!
! This file is part of LIBPFASST.
!

!> N-dimensional complex array encapsulation.
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
!! extract pointers to the encapsulated array  without
!! performing any copies.
!!
module pf_mod_zndarray
  use iso_c_binding  
  use pf_mod_dtype
  use pf_mod_utils  
  implicit none

  !>  Factory for making zndarray
  type, extends(pf_factory_t) :: zndarray_factory
  contains
     procedure :: create_single => zndarray_create_single
     procedure :: create_array => zndarray_create_array
     procedure :: destroy_single => zndarray_destroy_single
     procedure :: destroy_array => zndarray_destroy_array
  end type zndarray_factory

  !>  Complex ndarray
  type, extends(pf_encap_t) :: zndarray
     integer :: dim
     integer,    allocatable :: shape(:)     
    complex(pfdp), allocatable :: flatarray(:)
  contains
    procedure :: setval => zndarray_setval
    procedure :: copy => zndarray_copy
    procedure :: norm => zndarray_norm
    procedure :: pack => zndarray_pack
    procedure :: unpack => zndarray_unpack
    procedure :: axpy => zndarray_axpy
    procedure :: eprint => zndarray_eprint
    procedure :: write_to_disk
  end type zndarray

  contains

  function cast_as_zndarray(encap_polymorph) result(zndarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(zndarray), pointer :: zndarray_obj

    select type(encap_polymorph)
    type is (zndarray)
       zndarray_obj => encap_polymorph
    end select
  end function cast_as_zndarray

  !> Allocates complex ndarray
  subroutine zndarray_build(q, shape)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape(:)

    type(zndarray), pointer :: zndarray_obj

    select type (q)
    class is (zndarray)
       allocate(q%shape(size(shape)))
       allocate(q%flatarray(product(shape)))
       q%dim   = size(shape)
       q%shape = shape
       q%flatarray = cmplx(0.0, 0.0,pfdp)
       
    end select


    nullify(zndarray_obj)
  end subroutine zndarray_build

  subroutine zndarray_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(zndarray), pointer :: zndarray_obj


    zndarray_obj => cast_as_zndarray(encap)
    deallocate(zndarray_obj%shape)
    deallocate(zndarray_obj%flatarray)
    nullify(zndarray_obj)

  end subroutine zndarray_destroy

  !> Wrapper routine for allocation of a single zndarray type array
  subroutine zndarray_create_single(this, x, level,  shape)
    class(zndarray_factory), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x
    integer, intent(in) :: level, shape(:)

    allocate(zndarray::x)
    call zndarray_build(x, shape)

  end subroutine zndarray_create_single

  !> Wrapper routine for looped allocation of many zndarray type arrays
  subroutine zndarray_create_array(this, x, n, level,  shape)
    class(zndarray_factory), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x(:)
    integer, intent(in) :: n, level,  shape(:)
    integer :: i

    allocate(zndarray::x(n))
    do i = 1, n
       call zndarray_build(x(i), shape)
    end do

  end subroutine zndarray_create_array

  subroutine zndarray_destroy_single(this, x)
    class(zndarray_factory), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x

    select type (x)
    class is (zndarray)
       deallocate(x%shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)

  end subroutine zndarray_destroy_single

  !> Wrapper routine for looped allocation of many zndarray type arrays
  subroutine zndarray_destroy_array(this, x)
    class(zndarray_factory), intent(inout)       :: this
    class(pf_encap_t), intent(inout),allocatable :: x(:)
    integer :: i

    select type(x)
    class is (zndarray)
       do i = 1,size(x)
          deallocate(x(i)%shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
    
  end subroutine zndarray_destroy_array


  !> Set solution value.
  subroutine zndarray_setval(this, val, flags)
    class(zndarray), intent(inout) :: this
    real(pfdp), intent(in) :: val
    integer, intent(in), optional :: flags
    complex(pfdp) :: zval

    zval = cmplx(val, 0.0, pfdp)
    this%flatarray = zval
  end subroutine zndarray_setval

  !> Copy solution value.
  subroutine zndarray_copy(this, src, flags)
    class(zndarray), intent(inout) :: this
    class(pf_encap_t), intent(in) :: src
    integer, intent(in), optional :: flags
    class(zndarray), pointer :: zndarray_src

    zndarray_src => cast_as_zndarray(src)

    this%flatarray =  zndarray_src%flatarray
  end subroutine zndarray_copy

  !> Pack solution q into a flat array.
  subroutine zndarray_pack(this, z,flags)
    class(zndarray), intent(in) :: this
    real(pfdp), intent(out) :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: i

    do i = 1,product(this%shape)
       z(2*i-1) = real(this%flatarray(i))
       z(2*i)    = aimag(this%flatarray(i))
    end do
  end subroutine zndarray_pack

  ! Unpack solution from a flat array.
  subroutine zndarray_unpack(this, z,flags)
    class(zndarray), intent(inout) :: this
    real(pfdp), intent(in) :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: i

    do i = 1,product(this%shape)
       this%flatarray(i) = cmplx(z(2*i-1), z(2*i), pfdp)
    enddo
  end subroutine zndarray_unpack

  ! Compute norm of solution
  function zndarray_norm(this,flags) result (norm)
    class(zndarray), intent(in) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm

    norm = maxval(abs(this%flatarray))
  end function zndarray_norm

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine zndarray_axpy(this, a, x, flags)
    class(zndarray), intent(inout) :: this
    class(pf_encap_t), intent(in) :: x
    real(pfdp), intent(in) :: a
    integer, intent(in), optional :: flags
    class(zndarray), pointer :: zndarray_obj

    zndarray_obj => cast_as_zndarray(x)
    this%flatarray = a * zndarray_obj%flatarray + this%flatarray
  end subroutine zndarray_axpy

  subroutine zndarray_eprint(this,flags)
    class(zndarray), intent(inout) :: this
    integer,           intent(in   ), optional :: flags

    print*, this%flatarray(1:10)
  end subroutine zndarray_eprint

  subroutine write_to_disk(this, filename)
    class(zndarray), intent(inout) :: this
    character(len=*), intent(in) :: filename

    complex(pfdp),      pointer :: z_array(:,:)    
    open(unit=1, file=trim(filename), form='unformatted')

    write(1) this%flatarray
    !z_array=>get_array2d(this)
    !write(1) z_array
    close(1)
  end subroutine write_to_disk

    !>  Helper function to return the array part
  function get_array1d(x,flags) result(r)
    class(pf_encap_t), target, intent(in) :: x
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:)
    select type (x)
    type is (zndarray)
       r => x%flatarray
    end select
  end function get_array1d
  

  function get_array2d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:)

    select type (x)
    type is (zndarray)
       r(1:x%shape(1),1:x%shape(2)) => x%flatarray
    end select
  end function get_array2d
  

  function get_array3d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (zndarray)
       r(1:x%shape(1),1:x%shape(2),1:x%shape(3)) => x%flatarray
    end select
  end function get_array3d

end module pf_mod_zndarray
