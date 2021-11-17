!!  N-dimensional array encapsulation.
!
! This file is part of LIBPFASST.
!

!> Module to define and encapsulation for an N-dimensional array 
!!
!! When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'lev_shape' attribute to create a new array with that
!! shape.  Thus, the 'lev_shape' attributes of the PFASST levels should be
!! set appropriately.  For example, before calling pf_pfasst_run we can
!! set the shape of the coarsest level by doing:
!!
!!   allocate(pf%levels(1)%lev_shape(2))
!!   pf%levels(1)%lev_shape = [ 3, 10 ]
!!
!! The helper routine get_array, etc can be used to
!! extract pointers to the encapsulated array without
!! performing any copies.
!!
module pf_mod_ndarray
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_stop
  implicit none

  !>  Type to create and destroy N-dimensional arrays
  type, extends(pf_factory_t) :: pf_ndarray_factory_t
   contains
     procedure :: create_single  => ndarray_create_single
     procedure :: create_array  => ndarray_create_array
     procedure :: destroy_single => ndarray_destroy_single
     procedure :: destroy_array => ndarray_destroy_array
  end type pf_ndarray_factory_t

  !>  N-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: pf_ndarray_t
     integer             :: ndim
     integer,    allocatable :: arr_shape(:)
     real(pfdp), allocatable :: flatarray(:)
   contains
     procedure :: setval => ndarray_setval
     procedure :: copy => ndarray_copy
     procedure :: norm => ndarray_norm
     procedure :: pack => ndarray_pack
     procedure :: unpack => ndarray_unpack
     procedure :: axpy => ndarray_axpy
     procedure :: eprint => ndarray_eprint
     procedure, private  :: get_array_1d  ,get_array_2d,get_array_3d,get_array_4d
     generic :: get_array => get_array_1d ,get_array_2d,get_array_3d,get_array_4d
  end type pf_ndarray_t
  
contains

  function cast_as_ndarray(encap_polymorph) result(ndarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pf_ndarray_t), pointer :: ndarray_obj
    
    select type(encap_polymorph)
    type is (pf_ndarray_t)
       ndarray_obj => encap_polymorph
    end select
  end function cast_as_ndarray
  
  !>  Subroutine to allocate the array and set the size parameters
  subroutine ndarray_build(q, shape_in)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape_in(:)
    
    integer :: ierr
    select type (q)
    class is (pf_ndarray_t)
       allocate(q%arr_shape(SIZE(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
       allocate(q%flatarray(product(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)       
       q%ndim   = SIZE(shape_in)
       q%arr_shape = shape_in
       q%flatarray = 0.0_pfdp       
    end select
  end subroutine ndarray_build
  
  !> Subroutine to  create a single array
  subroutine ndarray_create_single(this, x, level_index, lev_shape)
    class(pf_ndarray_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: ierr
    allocate(pf_ndarray_t::x,stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)       
    
    call ndarray_build(x, lev_shape)
    
  end subroutine ndarray_create_single

  !> Subroutine to create an array of arrays
  subroutine ndarray_create_array(this, x, n, level_index,  lev_shape)
    class(pf_ndarray_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: i,ierr
    
    allocate(pf_ndarray_t::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)       
    
    do i = 1, n
       call ndarray_build(x(i), lev_shape)
    end do
    
  end subroutine ndarray_create_array

  !>  Subroutine to destroy array
  subroutine ndarray_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_ndarray_t), pointer :: ndarray_obj

    ndarray_obj => cast_as_ndarray(encap)
    deallocate(ndarray_obj%arr_shape)
    deallocate(ndarray_obj%flatarray)
    nullify(ndarray_obj)

  end subroutine ndarray_destroy

  !> Subroutine to destroy an single array
  subroutine ndarray_destroy_single(this, x)
    class(pf_ndarray_factory_t), intent(inout)   :: this
    class(pf_encap_t),      intent(inout), allocatable :: x

    select type (x)
    class is (pf_ndarray_t)
       deallocate(x%arr_shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
    
  end subroutine ndarray_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine ndarray_destroy_array(this, x)
    class(pf_ndarray_factory_t), intent(inout)  :: this
    class(pf_encap_t), intent(inout),allocatable :: x(:)
    integer :: i

    select type(x)
    class is (pf_ndarray_t)
       do i = 1,SIZE(x)
          deallocate(x(i)%arr_shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
    
  end subroutine ndarray_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalar  value.
  subroutine ndarray_setval(this, val, flags)
    class(pf_ndarray_t), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%flatarray = val
  end subroutine ndarray_setval

  !> Subroutine to copy an array
  subroutine ndarray_copy(this, src, flags)
    class(pf_ndarray_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (pf_ndarray_t)
       this%flatarray = src%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine ndarray_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine ndarray_pack(this, z, flags)
    class(pf_ndarray_t), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    z = this%flatarray
  end subroutine ndarray_pack

  !> Subroutine to unpack to a flatarray after receiving
  subroutine ndarray_unpack(this, z, flags)
    class(pf_ndarray_t), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags
    this%flatarray = z
  end subroutine ndarray_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function ndarray_norm(this, flags) result (norm)
    class(pf_ndarray_t), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%flatarray))
  end function ndarray_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine ndarray_axpy(this, a, x, flags)
    class(pf_ndarray_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags

    if (a .eq. 0.0_pfdp) return
    select type(x)
    type is (pf_ndarray_t)
       if (a .eq. 1.0_pfdp) then
          this%flatarray = x%flatarray + this%flatarray
       elseif (a .eq. -1.0_pfdp) then
          this%flatarray = -x%flatarray + this%flatarray
       else
          this%flatarray = a*x%flatarray + this%flatarray
       end if
          
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine ndarray_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine ndarray_eprint(this,flags)
    class(pf_ndarray_t), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
    if (product(this%arr_shape) < 10) then
       print *, this%flatarray
    else
       print *, this%flatarray(1:10)
    endif
  end subroutine ndarray_eprint

  !>  Helper function to return the array part, these are called with get_array
  subroutine get_array_1d(this,r,flags) 
    class(pf_ndarray_t), target, intent(in) :: this
    real(pfdp), pointer, intent(inout) :: r(:)
    integer,    intent(in   ), optional :: flags
    r => this%flatarray
  end subroutine get_array_1d
  subroutine get_array_2d(this,r,flags) 
    class(pf_ndarray_t), target, intent(in) :: this
    real(pfdp), pointer, intent(inout) :: r(:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2)) => this%flatarray
  end subroutine get_array_2d
  subroutine get_array_3d(this,r,flags) 
    class(pf_ndarray_t), target, intent(in) :: this
    real(pfdp), pointer, intent(inout) :: r(:,:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2),1:this%arr_shape(3)) => this%flatarray
  end subroutine get_array_3d
  subroutine get_array_4d(this,r,flags) 
    class(pf_ndarray_t), target, intent(in) :: this
    real(pfdp), pointer, intent(inout) :: r(:,:,:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2),1:this%arr_shape(3),1:this%arr_shape(4)) => this%flatarray
  end subroutine get_array_4d

  !>  Helper function to return the array part
  function get_array1d(x,flags) result(r)
    class(pf_encap_t), target, intent(in) :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:)
    select type (x)
    type is (pf_ndarray_t)
       r => x%flatarray
    end select
  end function get_array1d
  

  function get_array2d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:)

    select type (x)
    type is (pf_ndarray_t)
       r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray
    end select
  end function get_array2d
  

  function get_array3d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (pf_ndarray_t)
       r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray
    end select
  end function get_array3d

end module pf_mod_ndarray
