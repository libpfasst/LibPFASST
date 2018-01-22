!
! Copyright (C) 2013 Matthew Emmett.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

!
! N-dimensional array encapsulation.
!
! When a new solution is created by a PFASST level, this encapsulation
! uses the levels 'shape' attribute to create a new array with that
! shape.  Thus, the 'shape' attributes of the PFASST levels should be
! set appropriately.  For example, before calling pf_pfasst_run we can
! set the shape of the coarsest level by doing:
!
!   allocate(pf%levels(1)%shape(2))
!   pf%levels(1)%shape = [ 3, 10 ]
!
! The helper routines array1, array2, array3, etc can be used to
! extract pointers to the encapsulated array from a C pointer without
! performing any copies.
!
!>  Module to encapsulate an N-dimensional array
module pf_mod_ndarray
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  !>  Type to create and destroy the arrays
  type, extends(pf_factory_t) :: ndarray_factory
   contains
     procedure :: create_single  => ndarray_create_single
     procedure :: create_array  => ndarray_create_array
     procedure :: destroy_single => ndarray_destroy_single
     procedure :: destroy_array => ndarray_destroy_array
  end type ndarray_factory

  !>  Type to extend the abstract encap and set procedure pointers
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
     subroutine ndarray_dump_numpy(dname, fname, endian, dim, shape, nvars, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname, fname, endian(5)
       integer,    intent(in   ), value :: dim, nvars
       integer,    intent(in   )        :: shape(dim)
       real(c_double),    intent(in   )        :: array(nvars)
     end subroutine ndarray_dump_numpy
  end interface

contains
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
  subroutine ndarray_create_single(this, x, level, kind, nvars, shape)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, kind, nvars, shape(:)
    integer :: i
    allocate(ndarray::x)
    call ndarray_build(x, shape)
  end subroutine ndarray_create_single

  !> Subroutine to create an array of arrays
  subroutine ndarray_create_array(this, x, n, level, kind, nvars, shape)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, kind, nvars, shape(:)
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
  subroutine ndarray_destroy_single(this, x, level, kind, nvars, shape)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, kind, nvars, shape(:)
    
    select type (x)
    class is (ndarray)
       deallocate(x%shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine ndarray_destroy_single

  
  !> Subroutine to destroy an array of arrays
  subroutine ndarray_destroy_array(this, x, n, level, kind, nvars, shape)
    class(ndarray_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, kind, nvars, shape(:)
    integer                                            :: i
    
    select type(x)
    class is (ndarray)
       do i = 1,n
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
       stop "TYPE ERROR"
    end select
  end subroutine ndarray_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine ndarray_pack(this, z)
    class(ndarray), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    z = this%flatarray
  end subroutine ndarray_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine ndarray_unpack(this, z)
    class(ndarray), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    this%flatarray = z
  end subroutine ndarray_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function ndarray_norm(this) result (norm)
    class(ndarray), intent(in   ) :: this
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
       stop "TYPE ERROR"
    end select
  end subroutine ndarray_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine ndarray_eprint(this)
    class(ndarray), intent(inout) :: this
    !  Just print the first few values
    print *, this%flatarray(1:10)
  end subroutine ndarray_eprint

  ! ! Helpers
  ! function dims(solptr) result(r)
  !   type(ndarray), intent(in   ), value :: solptr
  !   integer :: r

  !   type(ndarray), pointer :: sol
  !   call c_f_pointer(solptr, sol)

  !   r = sol%dim
  ! end function dims

  function cast_as_ndarray(encap_polymorph) result(ndarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(ndarray), pointer :: ndarray_obj
    
    select type(encap_polymorph)
    type is (ndarray)
       ndarray_obj => encap_polymorph
    end select
  end function cast_as_ndarray

  function array1(x) result(r)
    class(pf_encap_t), intent(in) :: x
    real(pfdp), pointer :: r(:)
    select type (x)
    type is (ndarray)
       r => x%flatarray
    end select
  end function array1

  ! function array2(solptr) result(r)
  !   type(ndarray), intent(in   ), value :: solptr
  !   real(pfdp), pointer :: r(:,:)

  !   integer                :: shp(2)
  !   type(ndarray), pointer :: sol
  !   call c_f_pointer(solptr, sol)

  !   if (sol%dim == 2) then
  !      shp = sol%shape
  !      call c_f_pointer(sol%aptr, r, shp)
  !   else
  !      stop "ERROR: array2 dimension mismatch."
  !   end if
  ! end function array2

  ! function array3(solptr) result(r)
  !   type(ndarray), intent(in   ), value :: solptr
  !   real(pfdp), pointer :: r(:,:,:)

  !   integer                :: shp(3)
  !   type(ndarray), pointer :: sol
  !   call c_f_pointer(solptr, sol)

  !   if (sol%dim == 3) then
  !      shp = sol%shape
  !      call c_f_pointer(sol%aptr, r, shp)
  !   else
  !      stop "ERROR: array3 dimension mismatch."
  !   end if
  ! end function array3

  ! function array4(solptr) result(r)
  !   type(ndarray), intent(in   ), value :: solptr
  !   real(pfdp), pointer :: r(:,:,:,:)

  !   integer                :: shp(4)
  !   type(ndarray), pointer :: sol
  !   call c_f_pointer(solptr, sol)

  !   if (sol%dim == 4) then
  !      shp = sol%shape
  !      call c_f_pointer(sol%aptr, r, shp)
  !   else
  !      stop "ERROR: array4 dimension mismatch."
  !   end if
  ! end function array4

  ! function array5(solptr) result(r)
  !   type(ndarray), intent(in   ), value :: solptr
  !   real(pfdp), pointer :: r(:,:,:,:,:)

  !   integer                :: shp(5)
  !   type(ndarray), pointer :: sol
  !   call c_f_pointer(solptr, sol)

  !   if (sol%dim == 5) then
  !      shp = sol%shape
  !      call c_f_pointer(sol%aptr, r, shp)
  !   else
  !      stop "ERROR: array5 dimension mismatch."
  !   end if
  ! end function array5

  ! subroutine ndarray_encap_create(encap)
  !   type(pf_encap_t), intent(out) :: encap

  !   encap%create  => ndarray_create
  !   encap%destroy => ndarray_destroy
  !   encap%setval  => ndarray_setval
  !   encap%copy    => ndarray_copy
  !   encap%norm    => ndarray_norm
  !   encap%pack    => ndarray_pack
  !   encap%unpack  => ndarray_unpack
  !   encap%axpy    => ndarray_saxpy
  ! end subroutine ndarray_encap_create

  ! subroutine ndarray_dump_hook(pf, level, state)
  !   type(pf_pfasst_t),   intent(inout) :: pf
  !   type(pf_level_t),    intent(inout) :: level
  !   type(pf_state_t),    intent(in)    :: state

  !   character(len=256)     :: fname
  !   type(ndarray), pointer :: qend

  !   call c_f_pointer(level%qend, qend)

  !   write(fname, "('s',i0.5,'i',i0.3,'l',i0.2,'.npy')") &
  !        state%step, state%iter, level%level

  !   call ndarray_dump_numpy(trim(pf%outdir)//c_null_char, trim(fname)//c_null_char, '<f8'//c_null_char, &
  !        qend%dim, qend%shape, size(qend%flatarray), qend%flatarray)

  ! end subroutine ndarray_dump_hook


end module pf_mod_ndarray
