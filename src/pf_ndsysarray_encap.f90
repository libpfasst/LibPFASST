!!  N-dimensional system of array encapsulation.
!
! This file is part of LIBPFASST.
!

!> Module to define and encapsulation for a system of N-dimensional arrays 
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
!! The helper routines get_array1d, get_array2d, get_array3d, etc can be used to
!! extract pointers to a component of  encapsulated system
!! performing any copies or the whole system if desired.
!!
module pf_mod_ndsysarray
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_stop
  implicit none

  !>  Type to create and destroy systems of N-dimensional arrays
  type, extends(pf_factory_t) :: pf_ndsysarray_factory_t
   contains
     procedure :: create_single  => ndsysarray_create_single
     procedure :: create_array  => ndsysarray_create_array
     procedure :: destroy_single => ndsysarray_destroy_single
     procedure :: destroy_array => ndsysarray_destroy_array
  end type pf_ndsysarray_factory_t

  !> Type for system of  N-dimensional arrays,  extends the abstract encap type  
  type, extends(pf_encap_t) :: pf_ndsysarray_t
     integer             :: ndim    !  The spatial dimension of each component in system
     integer             :: ncomp  !  The number of components in the system
     integer             :: ndof   !  The number of variables in each component
     integer,    allocatable :: arr_shape(:)
     real(pfdp), allocatable :: flatarray(:)
   contains
     procedure :: setval => ndsysarray_setval
     procedure :: copy => ndsysarray_copy
     procedure :: norm => ndsysarray_norm
     procedure :: pack => ndsysarray_pack
     procedure :: unpack => ndsysarray_unpack
     procedure :: axpy => ndsysarray_axpy
     procedure :: eprint => ndsysarray_eprint
  end type pf_ndsysarray_t

contains
  !>  Subroutine to allocate the array and set the size parameters
  subroutine ndsysarray_build(q, arr_shape)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: arr_shape(:)

    integer :: ierr
    select type (q)
    class is (pf_ndsysarray_t)
       allocate(q%arr_shape(SIZE(arr_shape)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                         
       
       q%ndim   = SIZE(arr_shape)-1
       q%ncomp = arr_shape(q%ndim+1)
       q%ndof = product(arr_shape(1:q%ndim))
       q%arr_shape = arr_shape

       allocate(q%flatarray(product(arr_shape)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                                
    end select
  end subroutine ndsysarray_build

  !> Subroutine to  create a single array
  subroutine ndsysarray_create_single(this, x, level_index, lev_shape)
    class(pf_ndsysarray_factory_t), intent(inout)      :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)

    integer :: ierr
    allocate(pf_ndsysarray_t::x,stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    call ndsysarray_build(x, lev_shape)
  end subroutine ndsysarray_create_single

  !> Subroutine to create an array of arrays
  subroutine ndsysarray_create_array(this, x, n, level_index,  lev_shape)
    class(pf_ndsysarray_factory_t), intent(inout)      :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: i,ierr
    allocate(pf_ndsysarray_t::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                                 
    do i = 1, n
       call ndsysarray_build(x(i), lev_shape)
    end do
  end subroutine ndsysarray_create_array
!!$
  !>  Subroutine to destroy array (simple)
  subroutine ndsysarray_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_ndsysarray_t), pointer :: ndsysarray_obj

    ndsysarray_obj => cast_as_ndsysarray(encap)

    deallocate(ndsysarray_obj%arr_shape)
    deallocate(ndsysarray_obj%flatarray)

    nullify(ndsysarray_obj)

  end subroutine ndsysarray_destroy

  !> Subroutine to destroy an single array
  subroutine ndsysarray_destroy_single(this, x)
    class(pf_ndsysarray_factory_t), intent(inout)      :: this
    class(pf_encap_t),      intent(inout), allocatable :: x

    select type (x)
    class is (pf_ndsysarray_t)
       deallocate(x%arr_shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine ndsysarray_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine ndsysarray_destroy_array(this, x)
    class(pf_ndsysarray_factory_t), intent(inout)      :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer                                            :: i

    select type(x)
    class is (pf_ndsysarray_t)
       do i = 1,SIZE(x)
          deallocate(x(i)%arr_shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
  end subroutine ndsysarray_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalare  value.
  subroutine ndsysarray_setval(this, val, flags)
    class(pf_ndsysarray_t), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%flatarray = val
  end subroutine ndsysarray_setval

  !> Subroutine to copy an array
  subroutine ndsysarray_copy(this, src, flags)
    class(pf_ndsysarray_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (pf_ndsysarray_t)
       this%flatarray = src%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'invalid type')              
    end select
  end subroutine ndsysarray_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine ndsysarray_pack(this, z, flags)
    class(pf_ndsysarray_t), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    z = this%flatarray
  end subroutine ndsysarray_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine ndsysarray_unpack(this, z, flags)
    class(pf_ndsysarray_t), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags
    this%flatarray = z
  end subroutine ndsysarray_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function ndsysarray_norm(this, flags) result (norm)
    class(pf_ndsysarray_t), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%flatarray))
  end function ndsysarray_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine ndsysarray_axpy(this, a, x, flags)
    class(pf_ndsysarray_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags

    if (a .eq. 0.0_pfdp) return
    select type(x)
    type is (pf_ndsysarray_t)
       if (a .eq. 1.0_pfdp) then
          this%flatarray = x%flatarray + this%flatarray
       elseif (a .eq. -1.0_pfdp) then
          this%flatarray = -x%flatarray + this%flatarray
       else
          this%flatarray = a*x%flatarray + this%flatarray
       end if
    class default
       call pf_stop(__FILE__,__LINE__,'invalid type')              
    end select

  end subroutine ndsysarray_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine ndsysarray_eprint(this,flags)
    class(pf_ndsysarray_t), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
        print *, this%flatarray(1:10)
    !print *, this%flatarray
  end subroutine ndsysarray_eprint


  function cast_as_ndsysarray(encap_polymorph) result(ndsysarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pf_ndsysarray_t), pointer :: ndsysarray_obj
    
    select type(encap_polymorph)
    type is (pf_ndsysarray_t)
       ndsysarray_obj => encap_polymorph
    end select
  end function cast_as_ndsysarray

  !>  Helper function to return the array part or the whole thing
  function get_array1d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:)
    select type (x)
    type is (pf_ndsysarray_t)
       if (n .eq. 0) then    !  Return pointer to whole array
          if (x%ncomp .eq. 1 .and. x%ndim .eq. 1) then
             r => x%flatarray
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
          end if
       else                  !  Return pointer to nth component
          if (x%ndim .eq. 1) then
             r => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
          end if
       end if
    end select
  end function get_array1d
  

  function get_array2d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n    
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:)


    select type (x)
    type is (pf_ndsysarray_t)
       if (n .eq. 0) then
          if (x%ndim .eq. 1) then  !  Return pointer to whole array
             r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
          end if
       else                      !  Return pointer to nth component
          if (x%ndim .eq. 2) then
             r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
          end if
          
       endif
       
    end select
  end function get_array2d
  

  function get_array3d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (pf_ndsysarray_t)
       if (n .eq. 0) then   !  Return pointer to whole array
          if (x%ndim .eq. 2) then
             r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
          end if
       else         !  Return pointer to nth component
          if (x%ndim .eq. 3) then
             r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 3. ndim=',x%ndim)
          end if
       end if
    end select
  end function get_array3d

  function get_array4d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n
    integer,           intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:,:)

    select type (x)
    type is (pf_ndsysarray_t)
       if (n .eq. 0) then  !  Return pointer to whole array
          if (x%ndim .eq. 3) then
             r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3),1:x%arr_shape(4)) => x%flatarray
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 3. ndim=',x%ndim)
          end if
       else               !  Return pointer to nth component
          if (x%ndim .eq. 4) then
             r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3),1:x%arr_shape(4)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 4. ndim=',x%ndim)
          end if
       end if
    end select
  end function get_array4d
  
end module pf_mod_ndsysarray
