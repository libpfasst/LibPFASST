!
! This file is part of LIBPFASST.
!

!> Simple scalar encapsulation example
!!

module encap
   use pfasst
   implicit none
 
   !>  Type to create and destroy the local data encapsulation
   type, extends(pf_factory_t) :: scalar_factory
      type(c_ptr) :: c_factory_ptr = c_null_ptr ! c pointer
    contains
      procedure :: create_single  => scalar_create_single
      procedure :: create_array  => scalar_create_array
      procedure :: destroy_single => scalar_destroy_single
      procedure :: destroy_array => scalar_destroy_array
   end type scalar_factory
 
   !>  Type to extend the abstract encap and set procedure pointers
   type, extends(pf_encap_t) :: scalar_encap
      type(c_ptr) :: c_encap_ptr = c_null_ptr ! c pointer
    contains
      procedure :: setval => scalar_setval
      procedure :: copy => scalar_copy
      procedure :: norm => scalar_norm
      procedure :: pack => scalar_pack
      procedure :: unpack => scalar_unpack
      procedure :: axpy => scalar_axpy
      procedure :: eprint => scalar_eprint
      procedure :: getval
   end type scalar_encap

   interface

      subroutine ScalarCreate(x) bind(c, name="ScalarCreate")
         use iso_c_binding
         type(c_ptr) :: x
      end subroutine ScalarCreate
    
      subroutine ScalarDestroy(x) bind(c, name="ScalarDestroy")
         use iso_c_binding
         type(c_ptr), value :: x
      end subroutine ScalarDestroy
    
      subroutine ScalarSetVal(x, val) bind(c, name="ScalarSetVal")
         use iso_c_binding
         type(c_ptr), value:: x
         real(c_double), value :: val
      end subroutine ScalarSetVal
    
      subroutine ScalarCopy(dest, src) bind(c, name="ScalarCopy")
         use iso_c_binding
         type(c_ptr), value :: dest, src
      end subroutine ScalarCopy
   
      function ScalarPack(x) result(z) bind(c, name="ScalarPack")
         use iso_c_binding
         type(c_ptr), value :: x
         real(c_double) :: z
      end function
 
      subroutine ScalarUnpack(x, z) bind(c, name="ScalarUnpack")
         use iso_c_binding
         type(c_ptr), value :: x
         real(c_double), value :: z
      end subroutine ScalarUnpack
    
      function ScalarNorm(x) result(norm) bind(c, name="ScalarNorm")
        use iso_c_binding
        type(c_ptr), value :: x
        real(c_double) :: norm
      end function
    
      subroutine ScalarAxpy(y, a, x) bind(c, name="ScalarAxpy")
        use iso_c_binding
        type(c_ptr), value :: x, y
        real(c_double), value  :: a
      end subroutine ScalarAxpy
    
      subroutine ScalarPrint(x) bind(c, name="ScalarPrint")
        use iso_c_binding
        type(c_ptr), value :: x
      end subroutine ScalarPrint

      function ScalarGetVal(x) result(val) bind(c, name="ScalarGetVal")
        use iso_c_binding
        type(c_ptr), value :: x
        real(c_double) :: val
      end function
   end interface

contains

  !>  The following are the base subroutines that encapsulation factories need to provide
  
  !>  Subroutine to allocate one encap
  subroutine scalar_create_single(this, x, level_index, lev_shape)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x
    integer,               intent(in   ) ::  level_index ! passed by default,  not needed here
    integer,               intent(in   ) ::  lev_shape(:) ! passed by default, not needed here
    integer :: ierr

    allocate(scalar_encap::x, stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)

    select type(x)
    type is (scalar_encap)
       call ScalarCreate(x%c_encap_ptr)
    end select
  end subroutine scalar_create_single

  !> Subroutine to create an array of encaps
  subroutine scalar_create_array(this, x, n, level_index,lev_shape)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x(:)
    integer,               intent(in   )              :: n  ! size of array to build
    integer,               intent(in   ) ::  level_index ! passed by default,  not needed here
    integer,               intent(in   ) ::  lev_shape(:) ! passed by default, not needed here
    integer :: i, ierr

    allocate(scalar_encap::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',n)

    select type(x)
    type is (scalar_encap)
       do i = 1, n
           call ScalarCreate(x(i)%c_encap_ptr)
       end do
    end select
  end subroutine scalar_create_array

  !> Subroutine to destroy a single array encap
  subroutine scalar_destroy_single(this, x)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x
    integer                                           :: ierr
   
    select type(x)
    type is (scalar_encap)
       call ScalarDestroy(x%c_encap_ptr)
    end select
    deallocate(x,stat=ierr)
  end subroutine scalar_destroy_single

  !> Subroutine to destroy an array of arrays
  subroutine scalar_destroy_array(this, x)
    class(scalar_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x(:)
    integer                                           :: i, ierr
    class(scalar_encap), pointer :: x_ptr

    select type(x)
    type is (scalar_encap)
       do i = 1, size(x)
          call ScalarDestroy(x(i)%c_encap_ptr)
       end do
    end select
    deallocate(x,stat=ierr)
  end subroutine scalar_destroy_array

  !>  The following are the base subroutines that all encapsulations must provide

  !> Subroutine to set array to a scalar  value.
  subroutine scalar_setval(this, val, flags)
    class(scalar_encap), intent(inout)      :: this
    real(c_double),     intent(in   )       :: val
    integer,        intent(in   ), optional :: flags
    call ScalarSetVal(this%c_encap_ptr, val)
  end subroutine scalar_setval

  !> Subroutine to copy an array
  subroutine scalar_copy(this, src, flags)
    class(scalar_encap),    intent(inout)      :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    class(scalar_encap), pointer               :: src_scalar_encap

    select type(src)
    type is (scalar_encap)
       call ScalarCopy(this%c_encap_ptr, src%c_encap_ptr)
    end select
  end subroutine scalar_copy

  !> Subroutine to pack into a flat array for sending
  subroutine scalar_pack(this, z, flags)
    class(scalar_encap), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags

    z(1) = ScalarPack(this%c_encap_ptr)
  end subroutine scalar_pack

  !> Subroutine to unpack  after receiving
  subroutine scalar_unpack(this, z, flags)
     class(scalar_encap), intent(inout) :: this
     real(pfdp),     intent(in   ) :: z(:)
     integer,     intent(in   ), optional :: flags

     call ScalarUnpack(this%c_encap_ptr, z(1));
  end subroutine scalar_unpack

  !> Subroutine to define the norm of the array (here the abs value)
  function scalar_norm(this, flags) result (norm)
    class(scalar_encap), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(c_double) :: norm
    ! call ScalarNorm(this%c_encap_ptr, norm)
    norm = ScalarNorm(this%c_encap_ptr)
  end function scalar_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine scalar_axpy(this, a, x, flags)
    class(scalar_encap), intent(inout) :: this
    class(pf_encap_t), intent(in) :: x
    real(c_double), intent(in) :: a
    integer, intent(in), optional :: flags

    select type(x)
    type is (scalar_encap) 
       call ScalarAxpy(this%c_encap_ptr, a, x%c_encap_ptr)
    end select
  end subroutine scalar_axpy

  !> Jordi stopped here
  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine scalar_eprint(this,flags)
    class(scalar_encap), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Print the  value
    call ScalarPrint(this%c_encap_ptr)
  end subroutine scalar_eprint

  function getval(this) result(val)
     class(scalar_encap), intent(inout) :: this
     real(pfdp) :: val
     val = ScalarGetVal(this%c_encap_ptr)
  end function

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
