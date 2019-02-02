!--------------------------------------------------------------------------------
! Copyright (c) 2017, Brandon Krull.  All rights reserved.
!--------------------------------------------------------------------------------
! MODULE: factory
!
!> @author
!> Brandon Krull, Berkeley Lab
!
! Description:
!> This module defines the (previously encap) factory which defines
!! the time-dependent Kohn-Sham equations or simpler model time-dependent problems
module mod_zmkpair
  use pf_mod_dtype

  implicit none

  type, extends(pf_factory_t) :: zmkpair_factory
  contains
     procedure :: create_single => zmkpair_create_single
     procedure :: create_array => zmkpair_create_array
     procedure :: destroy_single => zmkpair_destroy_single
     procedure :: destroy_array => zmkpair_destroy_array
  end type zmkpair_factory

  type, extends(pf_encap_t) :: zmkpair
    integer :: dim
    complex(pfdp), allocatable :: array(:,:), y(:,:)
  contains
    procedure :: setval => zmkpair_setval
    procedure :: copy => zmkpair_copy
    procedure :: norm => zmkpair_norm
    procedure :: pack => zmkpair_pack
    procedure :: unpack => zmkpair_unpack
    procedure :: axpy => zmkpair_axpy
    procedure :: eprint => zmkpair_eprint
    procedure :: write_to_disk
  end type zmkpair

  contains

  function cast_as_zmkpair(encap_polymorph) result(zmkpair_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(zmkpair), pointer :: zmkpair_obj

    select type(encap_polymorph)
    type is (zmkpair)
       zmkpair_obj => encap_polymorph
    end select
  end function cast_as_zmkpair

  !> Allocates complex mkpair
  subroutine zmkpair_build(encap, dim)
    class(pf_encap_t), intent(inout) :: encap
    integer, intent(in) ::  dim

    type(zmkpair), pointer :: zmkpair_obj

    zmkpair_obj => cast_as_zmkpair(encap)
    allocate(zmkpair_obj%array(dim, dim))
    allocate(zmkpair_obj%y(dim, dim))

    zmkpair_obj%dim = dim
    zmkpair_obj%array(:,:) = cmplx(0.0, 0.0,pfdp)
    zmkpair_obj%y(:,:) = cmplx(0.0, 0.0,pfdp)

    nullify(zmkpair_obj)
  end subroutine zmkpair_build

  subroutine zmkpair_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(zmkpair), pointer :: zmkpair_obj

    zmkpair_obj => cast_as_zmkpair(encap)

    deallocate(zmkpair_obj%array)
    deallocate(zmkpair_obj%y)
    nullify(zmkpair_obj)

  end subroutine zmkpair_destroy

  !> Wrapper routine for allocation of a single zmkpair type array
  subroutine zmkpair_create_single(this, x, level, shape)
    class(zmkpair_factory), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x
    integer, intent(in) :: level, shape(:)

    allocate(zmkpair::x)
    call zmkpair_build(x, shape(1))

  end subroutine zmkpair_create_single

  !> Wrapper routine for looped allocation of many zmkpair type arrays
  subroutine zmkpair_create_array(this, x, n, level, shape)
    class(zmkpair_factory), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x(:)
    integer, intent(in) :: n, level, shape(:)
    integer :: i

    allocate(zmkpair::x(n))
    do i = 1, n
       call zmkpair_build(x(i), shape(1))
    end do

  end subroutine zmkpair_create_array

  subroutine zmkpair_destroy_single(this, x, level, shape)
    class(zmkpair_factory), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: x
    integer, intent(in) :: level, shape(:)
    class(zmkpair), pointer :: zmkpair_obj

    zmkpair_obj => cast_as_zmkpair(x)
    deallocate(zmkpair_obj%array)
    deallocate(zmkpair_obj%y)

  end subroutine zmkpair_destroy_single

  !> Wrapper routine for looped allocation of many zmkpair type arrays
  subroutine zmkpair_destroy_array(this, x, n, level, shape)
    class(zmkpair_factory), intent(inout):: this
    class(pf_encap_t), intent(inout), pointer :: x(:)
    integer, intent(in) :: n, level, shape(:)
    class(zmkpair), pointer :: zmkpair_obj
    integer :: i

    do i = 1, n
      zmkpair_obj => cast_as_zmkpair(x(i))
      deallocate(zmkpair_obj%array)
      deallocate(zmkpair_obj%y)
    end do

!    deallocate(x)
  end subroutine zmkpair_destroy_array

  !> Set solution value.
  subroutine zmkpair_setval(this, val, flags)
    class(zmkpair), intent(inout) :: this
    real(pfdp), intent(in) :: val
    integer, intent(in), optional :: flags
    complex(pfdp) :: zval

    zval = cmplx(val, 0.0, pfdp)
    if (present(flags)) then
       this%y =  zval
    else
       this%array =  zval
    endif
  end subroutine zmkpair_setval

  !> Copy solution value.
  subroutine zmkpair_copy(this, src, flags)
    class(zmkpair), intent(inout) :: this
    class(pf_encap_t), intent(in) :: src
    integer, intent(in), optional :: flags
    class(zmkpair), pointer :: zmkpair_src

    zmkpair_src => cast_as_zmkpair(src)

    if (present(flags)) then
        this%y = zmkpair_src%y
    else
        this%array = zmkpair_src%array
    endif
  end subroutine zmkpair_copy

  !> Pack solution q into a flat array.
  subroutine zmkpair_pack(this, z, flags)
    class(zmkpair), intent(in) :: this
    integer, intent(in), optional :: flags
    real(pfdp), intent(out) :: z(:)
    integer :: nx,ny,nxny,i,j,ij
    nx=this%dim
    ny=this%dim
    nxny=nx*ny*2
    do j = 1,ny
       do i = 1,nx
          ij = 2*((j-1)*nx + i)
          z(ij-1) =  real(this%y(i,j))
          z(ij) = aimag(this%y(i,j))
       end do
    end do
  end subroutine zmkpair_pack

  ! Unpack solution from a flat array.
  subroutine zmkpair_unpack(this, z, flags)
    class(zmkpair), intent(inout) :: this
    integer, intent(in), optional :: flags
    real(pfdp), intent(in) :: z(:)
    integer :: nx,ny,nxny,i,j,ij
    nx=this%dim
    ny=this%dim
    nxny=nx*ny*2
    do j = 1,ny
       do i = 1,nx
          ij = 2*((j-1)*nx + i)
          this%y(i,j) = cmplx(z(ij-1), z(ij), pfdp)
      enddo
    enddo
  end subroutine zmkpair_unpack

  ! Compute norm of solution
  function zmkpair_norm(this, flags) result (norm)
    class(zmkpair), intent(in) :: this
    integer, intent(in), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%array))
  end function zmkpair_norm

  ! Compute y = a x + y where a is a scalar and x and y are solutions.
  subroutine zmkpair_axpy(this, a, x, flags)
    class(zmkpair), intent(inout) :: this
    class(pf_encap_t), intent(in) :: x
    real(pfdp), intent(in) :: a
    integer, intent(in), optional :: flags
    class(zmkpair), pointer :: zmkpair_obj
    integer :: this_shape(2), dim, i

    zmkpair_obj => cast_as_zmkpair(x)
    if (present(flags)) then
       this%y = a * zmkpair_obj%y + this%y
    else
       this%array = a * zmkpair_obj%array + this%array
    endif

  end subroutine zmkpair_axpy

  subroutine zmkpair_eprint(this, flags)
    class(zmkpair), intent(inout) :: this
    integer, intent(in), optional :: flags
    integer :: this_shape(2), i, j, dim

    this_shape = shape(this%array)
    dim = this_shape(1)
    print*, 'omega'
    do i = 1, dim-1
       write(*, 100) real(this%array(i,i+1))
    end do
    print*, 'y'
    do i = 1, dim
       write(*, 100) real(this%y(i,i))
    end do
    100 format (*(F18.14))
    call flush
  end subroutine zmkpair_eprint

  subroutine write_to_disk(this, filename)
    class(zmkpair), intent(inout) :: this
    character(len=*), intent(in) :: filename

    open(unit=1, file=trim(filename), form='unformatted')
    write(1) this%y
    close(1)
  end subroutine write_to_disk

end module mod_zmkpair

