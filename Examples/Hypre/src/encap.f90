!
! This file is part of LIBPFASST.
!

!> Simple hypre_vector encapsulation example
!!

module encap
   use pfasst
   implicit none
 
   !>  Type to create and destroy the local data encapsulation
   type, extends(pf_factory_t) :: hypre_vector_factory
    contains
      procedure :: create_single  => hypre_vector_create_single
      procedure :: create_array  => hypre_vector_create_array
      procedure :: destroy_single => hypre_vector_destroy_single
      procedure :: destroy_array => hypre_vector_destroy_array
   end type hypre_vector_factory
 
   !>  Type to extend the abstract encap and set procedure pointers
   type, extends(pf_encap_t) :: hypre_vector_encap
      integer :: vector_size
      type(c_ptr) :: c_hypre_vector_ptr = c_null_ptr ! c pointer
    contains
      procedure :: setval => hypre_vector_setval
      procedure :: copy => hypre_vector_copy
      procedure :: norm => hypre_vector_norm
      procedure :: pack => hypre_vector_pack
      procedure :: unpack => hypre_vector_unpack
      procedure :: axpy => hypre_vector_axpy
      procedure :: eprint => hypre_vector_eprint
   end type hypre_vector_encap

   interface

      subroutine HypreVectorCreate(x_ptr, &
                                   level_index, &
                                   nx, &
                                   comm_color, &
                                   space_dim) bind(c, name="HypreVectorCreate")
         use iso_c_binding
         type(c_ptr) :: x_ptr
         integer, value :: level_index, nx, comm_color, space_dim
      end subroutine HypreVectorCreate
    
      subroutine HypreVectorDestroy(x) bind(c, name="HypreVectorDestroy")
         use iso_c_binding
         type(c_ptr), value :: x
      end subroutine HypreVectorDestroy
    
      subroutine HypreVectorSetVal(x, val) bind(c, name="HypreVectorSetVal")
         use iso_c_binding
         type(c_ptr), value :: x
         real(c_double), value :: val
      end subroutine HypreVectorSetVal
    
      subroutine HypreVectorCopy(dest, src) bind(c, name="HypreVectorCopy")
         use iso_c_binding
         type(c_ptr), value :: dest, src
      end subroutine HypreVectorCopy
   
      function HypreVectorPack(x) result(z) bind(c, name="HypreVectorPack")
         use iso_c_binding
         type(c_ptr), value :: x
         type(c_ptr) :: z
      end function
  
      subroutine HypreVectorUnpack(x, z) bind(c, name="HypreVectorUnpack")
         use iso_c_binding
         type(c_ptr), value :: x
         type(c_ptr), value :: z
      end subroutine HypreVectorUnpack
    
      function HypreVectorNorm(x) result(norm) bind(c, name="HypreVectorNorm")
        use iso_c_binding
        type(c_ptr), value :: x
        real(c_double) :: norm
      end function

      subroutine HypreVectorAxpy(y, a, x) bind(c, name="HypreVectorAxpy")
        use iso_c_binding
        type(c_ptr), value :: x, y
        real(c_double), value  :: a
      end subroutine HypreVectorAxpy
    
      subroutine HypreVectorPrint(x) bind(c, name="HypreVectorPrint")
        use iso_c_binding
        type(c_ptr), value :: x
      end subroutine HypreVectorPrint

      subroutine HypreVectorSetInitCond(x, val) bind(c, name="HypreVectorSetInitCond")
         use iso_c_binding
         type(c_ptr), value :: x
         real(c_double), value :: val
      end subroutine HypreVectorSetInitCond
      
      subroutine GetHypreStats() bind(c, name="GetHypreStats")
        use iso_c_binding
      end subroutine GetHypreStats
   end interface

contains

  !>  The following are the base subroutines that encapsulation factories need to provide
  
  !>  Subroutine to allocate one encap
  subroutine hypre_vector_create_single(this, x, level_index, lev_shape)
    class(hypre_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x
    integer,               intent(in   ) ::  level_index ! passed by default,  not needed here
    integer,               intent(in   ) ::  lev_shape(:) ! passed by default, not needed here
    integer :: ierr
    integer :: nx, comm_color, n_space, space_dim, max_space_v_cycles
    integer :: nrows, ilower0, ilower1, iupper0, iupper1

    nx = lev_shape(1)
    comm_color = lev_shape(2)
    space_dim = lev_shape(3)
    max_space_v_cycles = lev_shape(4)
    nrows = lev_shape(5)

    allocate(hypre_vector_encap::x, stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)

    select type(x)
    type is (hypre_vector_encap)
       call HypreVectorCreate(x%c_hypre_vector_ptr, level_index, nx, comm_color, space_dim)
       x%vector_size = nrows
    end select
  end subroutine hypre_vector_create_single

  !> Subroutine to create an array of encaps
  subroutine hypre_vector_create_array(this, x, n, level_index,lev_shape)
    class(hypre_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x(:)
    integer,               intent(in   )              :: n  ! size of array to build
    integer,               intent(in   ) ::  level_index ! passed by default,  not needed here
    integer,               intent(in   ) ::  lev_shape(:) ! passed by default, not needed here
    integer :: i, ierr
    integer :: nx, comm_color, n_space, space_dim, max_space_v_cycles
    integer :: nrows, ilower0, ilower1, iupper0, iupper1

    nx = lev_shape(1)
    comm_color = lev_shape(2)
    space_dim = lev_shape(3)
    max_space_v_cycles = lev_shape(4)
    nrows = lev_shape(5)

    allocate(hypre_vector_encap::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',n)

    select type(x)
    type is (hypre_vector_encap)
       do i = 1, n
           call HypreVectorCreate(x(i)%c_hypre_vector_ptr, level_index, nx, comm_color, space_dim)
           x(i)%vector_size = nrows
       end do
    end select
  end subroutine hypre_vector_create_array

  !> Subroutine to destroy a single array encap
  subroutine hypre_vector_destroy_single(this, x)
    class(hypre_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x
    integer                                           :: ierr
   
    select type(x)
    type is (hypre_vector_encap)
       call HypreVectorDestroy(x%c_hypre_vector_ptr)
    end select
    deallocate(x,stat=ierr)
  end subroutine hypre_vector_destroy_single

  !> Subroutine to destroy an array of arrays
  subroutine hypre_vector_destroy_array(this, x)
    class(hypre_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x(:)
    integer                                           :: i, ierr
    class(hypre_vector_encap), pointer :: x_ptr

    select type(x)
    type is (hypre_vector_encap)
       do i = 1, size(x)
          call HypreVectorDestroy(x(i)%c_hypre_vector_ptr)
       end do
    end select
    deallocate(x,stat=ierr)
  end subroutine hypre_vector_destroy_array

  !>  The following are the base subroutines that all encapsulations must provide

  !> Subroutine to set array to a hypre_vector  value.
  subroutine hypre_vector_setval(this, val, flags)
    class(hypre_vector_encap), intent(inout) :: this
    real(c_double), intent(in) :: val
    integer, intent(in), optional :: flags
    call HypreVectorSetVal(this%c_hypre_vector_ptr, val)
  end subroutine hypre_vector_setval

  !> Subroutine to copy an array
  subroutine hypre_vector_copy(this, src, flags)
    class(hypre_vector_encap),    intent(inout)      :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    class(hypre_vector_encap), pointer               :: src_hypre_vector_encap

    select type(src)
    type is (hypre_vector_encap)
       call HypreVectorCopy(this%c_hypre_vector_ptr, src%c_hypre_vector_ptr)
    end select
  end subroutine hypre_vector_copy

  !> Subroutine to pack into a flat array for sending
  subroutine hypre_vector_pack(this, z, flags)
    class(hypre_vector_encap), intent(in) :: this
    real(pfdp), intent(out) :: z(:)
    integer, intent(in), optional :: flags
    real(pfdp), pointer :: z_ptr(:)
    type(c_ptr) :: z_c_ptr

    z_c_ptr = HypreVectorPack(this%c_hypre_vector_ptr)
    call c_f_pointer(z_c_ptr, z_ptr, [this%vector_size])
    z = z_ptr
  end subroutine hypre_vector_pack

  !> Subroutine to unpack  after receiving
  subroutine hypre_vector_unpack(this, z, flags)
     class(hypre_vector_encap), intent(inout) :: this
     real(pfdp), intent(in) :: z(:)
     integer, intent(in), optional :: flags
     real(pfdp), target :: z2(size(z))
     type(c_ptr) :: z_c_ptr

     z2 = z
     z_c_ptr = c_loc(z2(1))
     call HypreVectorUnpack(this%c_hypre_vector_ptr, z_c_ptr);
  end subroutine hypre_vector_unpack

  !> Subroutine to define the norm of the array (here the abs value)
  function hypre_vector_norm(this, flags) result (norm)
    class(hypre_vector_encap), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(c_double) :: norm
    norm = HypreVectorNorm(this%c_hypre_vector_ptr)
  end function hypre_vector_norm

  !> Subroutine to compute y = a x + y where a is a hypre_vector and x and y are arrays
  subroutine hypre_vector_axpy(this, a, x, flags)
    class(hypre_vector_encap), intent(inout) :: this
    class(pf_encap_t), intent(in) :: x
    real(c_double), intent(in) :: a
    integer, intent(in), optional :: flags

    select type(x)
    type is (hypre_vector_encap) 
       call HypreVectorAxpy(this%c_hypre_vector_ptr, a, x%c_hypre_vector_ptr)
    end select
  end subroutine hypre_vector_axpy

  !> Jordi stopped here
  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine hypre_vector_eprint(this,flags)
    class(hypre_vector_encap), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Print the  value
    call HypreVectorPrint(this%c_hypre_vector_ptr)
  end subroutine hypre_vector_eprint

  !  Helper function to cast an abstract encap to the hypre_vector_encap
  function cast_as_hypre_vector(encap_polymorph) result(hypre_vector_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(hypre_vector_encap), pointer :: hypre_vector_obj
    
    select type(encap_polymorph)
    type is (hypre_vector_encap)
       hypre_vector_obj => encap_polymorph
    end select
  end function cast_as_hypre_vector

end module encap
