!! N-dimensional array encapsulation for optimal control problems.
!
! This file is part of LIBPFASST.
!
!!
!> When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'shape' attribute to create a new array with that
!! shape.  Thus, the 'shape' attributes of the PFASST levels should be
!! set appropriately.  For example, before calling pf_pfasst_run we can
!! set the shape of the coarsest level by doing:
!!
!!   allocate(pf%levels(1)%shape(2))
!!   pf%levels(1)%shape = [ 3, 10 ]
!!
!! The helper routines array1, array2, array3, etc can be used to
!! extract pointers to the encapsulated array from a C pointer without
!! performing any copies.
module pf_mod_ndarray_oc
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  type, extends(pf_factory_t) :: ndarray_oc_factory
   contains
     procedure :: create_single  => ndarray_oc_create_single
     procedure :: create_array   => ndarray_oc_create_array
     procedure :: destroy_single => ndarray_oc_destroy_single
     procedure :: destroy_array  => ndarray_oc_destroy_array
  end type ndarray_oc_factory
  
  type, extends(pf_encap_t) :: ndarray_oc
     integer                 :: dim
     integer,    allocatable :: shape(:)    
     real(pfdp), allocatable :: yflatarray(:)
     real(pfdp), allocatable :: pflatarray(:) 

   contains
     procedure :: setval => ndarray_oc_setval
     procedure :: copy   => ndarray_oc_copy
     procedure :: norm   => ndarray_oc_norm
     procedure :: pack   => ndarray_oc_pack
     procedure :: unpack => ndarray_oc_unpack
     procedure :: axpy   => ndarray_oc_axpy
     procedure :: eprint => ndarray_oc_eprint
  end type ndarray_oc

  ! interfaces to routines in pf_numpy.c
  interface
     subroutine ndarray_mkdir(dname, dlen) bind(c)
       use pf_mod_dtype       
       use iso_c_binding
       character(c_char), intent(in   )        :: dname
       integer(c_int),    intent(in   ), value :: dlen
     end subroutine ndarray_mkdir

     subroutine ndarray_dump_numpy(dname, fname, endian, dim, shape, nvars, array) bind(c)
       use iso_c_binding
       use pf_mod_dtype
       character(c_char), intent(in   )        :: dname, fname, endian(5)
       integer(c_int),    intent(in   ), value :: dim, nvars
       integer(c_int),    intent(in   )        :: shape(dim)
       real(pfdp),    intent(in   )        :: array(nvars)
     end subroutine ndarray_dump_numpy
  end interface

contains
  !>  Subroutine to allocate the array and set the size parameters
  subroutine ndarray_oc_build(q, shape)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape(:)

    select type (q)
    class is (ndarray_oc)
       allocate(q%shape(size(shape)))
       allocate(q%yflatarray(product(shape)))
       allocate(q%pflatarray(product(shape)))
       q%dim   = size(shape)
       q%shape = shape
    class default
        print *, "wrong class in ndarray_oc_build!"
        stop
    end select
  end subroutine ndarray_oc_build

  !> Subroutine to  create a single array
  subroutine ndarray_oc_create_single(this, x, level, shape)
    class(ndarray_oc_factory), intent(inout)           :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level, shape(:)
    allocate(ndarray_oc::x)
    call ndarray_oc_build(x, shape)
  end subroutine ndarray_oc_create_single

  !> Subroutine to create an array of arrays
  subroutine ndarray_oc_create_array(this, x, n, level, shape)
    class(ndarray_oc_factory), intent(inout)              :: this
    class(pf_encap_t),         intent(inout), allocatable :: x(:)
    integer,                   intent(in   )              :: n, level, shape(:)
    integer :: i
    allocate(ndarray_oc::x(n))
    do i = 1, n
       call ndarray_oc_build(x(i), shape)
    end do
  end subroutine ndarray_oc_create_array
  
  !>  Subroutine to destroy array
  subroutine ndarray_oc_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(ndarray_oc), pointer :: ndarray_oc_obj

    ndarray_oc_obj => cast_as_ndarray_oc(encap) !??

    deallocate(ndarray_oc_obj%pflatarray)
    deallocate(ndarray_oc_obj%yflatarray)
    deallocate(ndarray_oc_obj%shape)
    nullify(ndarray_oc_obj)
  end subroutine ndarray_oc_destroy

    !> Subroutine to destroy an single array
  subroutine ndarray_oc_destroy_single(this, x, level, shape)
    class(ndarray_oc_factory), intent(inout)              :: this
    class(pf_encap_t),         intent(inout), allocatable :: x
    integer,                   intent(in   )              :: level, shape(:)
    
    select type (x)
    class is (ndarray_oc)
       deallocate(x%pflatarray)
       deallocate(x%yflatarray)
       deallocate(x%shape)
    class default
      stop "TYPE ERROR in ndarray_oc_destroy_single"
    end select
    deallocate(x)
  end subroutine ndarray_oc_destroy_single
  
    !> Subroutine to destroy an array of arrays
  subroutine ndarray_oc_destroy_array(this, x, n, level, shape)
    class(ndarray_oc_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level, shape(:)
    integer                                            :: i
    
    select type(x)
    class is (ndarray_oc)
       do i = 1,n
          deallocate(x(i)%pflatarray)
          deallocate(x(i)%yflatarray)
          deallocate(x(i)%shape)
       end do
    class default
      stop "TYPE ERROR in ndarray_oc_destroy_array"
    end select
    deallocate(x)
  end subroutine ndarray_oc_destroy_array
  
  
  !> Subroutine to set array to a scalar  value.
  subroutine ndarray_oc_setval(this, val, flags)
    class(ndarray_oc), intent(inout)     :: this
    real(pfdp),  intent(in   )           :: val
    integer,     intent(in   ), optional :: flags
    integer :: which

    which = 0
    if (present(flags)) which = flags
!     if(.not.present(flags)) print *, "setval without flags"

    select case (which)
    case (0)
      this%yflatarray = val
      this%pflatarray = val
    case (1)
      this%yflatarray = val
    case (2)
      this%pflatarray = val
    case default
       stop "ERROR in ndarray_oc_setval: only 0, 1, 2 allowed as flags"
    end select

  end subroutine ndarray_oc_setval

  !> Subroutine to copy an array
  subroutine ndarray_oc_copy(this, src, flags)
    class(ndarray_oc),    intent(inout)  :: this
    class(pf_encap_t),    intent(in   )  :: src
    integer,     intent(in   ), optional :: flags
    integer :: which

    which = 0
    if (present(flags)) which = flags
!     if(.not.present(flags)) print *, "copy without flags"
    
    select type(src)
    type is (ndarray_oc)
      select case (which)
      case (0)
        this%yflatarray = src%yflatarray
        this%pflatarray = src%pflatarray
      case (1)
        this%yflatarray = src%yflatarray
      case (2)
        this%pflatarray = src%pflatarray
      case default
        stop "ERROR in ndarray_oc_copy: only 0, 1, 2 allowed as flags"
      end select
    class default
      stop "TYPE ERROR in ndarray_oc_copy"
    end select
  end subroutine ndarray_oc_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine ndarray_oc_pack(this, z, flags)
    class(ndarray_oc), intent(in   ) :: this
    real(pfdp),  intent(  out)       :: z(:)
    integer,     intent(in   ), optional   :: flags
    integer :: which
    
    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (0)
       !z = [sol%yflatarray, sol%pflatarray] 
       !z has to be right size? initialized to nvars, so it can hold either y or p
       !is it ever needed to pack y and p simultaneously?
       stop "ERROR in ndarray_oc_pack: only 1, 2 allowed as flags"
    case (1)
       z = this%yflatarray
    case (2)
       z = this%pflatarray
    case default
       stop "ERROR in ndarray_oc_pack: only 1, 2 allowed as flags"
    end select
  end subroutine ndarray_oc_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine ndarray_oc_unpack(this, z, flags)
    class(ndarray_oc), intent(inout) :: this
    real(pfdp),  intent(in   )        :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: which

    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (0)
       stop "ERROR in ndarray_oc_unpack: only 1, 2 allowed as flags"
    case (1)
       this%yflatarray = z
    case (2)
       this%pflatarray = z 
    case default
       stop "ERROR in ndarray_oc_upack: only 1, 2 allowed as flags"
    end select
  end subroutine ndarray_oc_unpack
  

  !> Subroutine to define the norm of the array (here the max norm)
  function ndarray_oc_norm(this, flags) result (norm)
    class(ndarray_oc), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp)   :: norm
    integer :: which

    which = 0
    if (present(flags)) which = flags
    if(.not.present(flags)) print *, "norm without flags"
    
    select case (which)
    case (0)
       norm = max(maxval(abs(this%yflatarray)), maxval(abs(this%pflatarray)))
    case (1)
       norm = maxval(abs(this%yflatarray))
    case (2)
       norm = maxval(abs(this%pflatarray))  
    case default
       stop "ERROR in ndarray_oc_norm: only 0, 1, 2 allowed as flags"
    end select
  end function ndarray_oc_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine ndarray_oc_axpy(this, a, x, flags)
    class(ndarray_oc), intent(inout) :: this    
    class(pf_encap_t), intent(in   )     :: x
    real(pfdp),  intent(in   )           :: a
    integer,     intent(in   ), optional :: flags
    integer :: which

    which = 0
    if (present(flags)) which = flags
!     if (.not.present(flags)) stop "axpy without flags" 
    
    select type(x)
    type is (ndarray_oc)
      select case (which)
      case (0)
        this%yflatarray = a * x%yflatarray + this%yflatarray
        this%pflatarray = a * x%pflatarray + this%pflatarray
      case (1)
        this%yflatarray = a * x%yflatarray + this%yflatarray
      case (2)
        this%pflatarray = a * x%pflatarray + this%pflatarray
      case default
        stop "ERROR in ndarray_oc_axpy: only 0, 1, 2 allowed as flags"
      end select
    class default
      stop "TYPE ERROR in ndarray_oc_axpy"
    end select  
  end subroutine ndarray_oc_axpy

  
  ! Helpers
!   function dims(solptr) result(r)
!     type(c_ptr), intent(in   ), value :: solptr
!     integer :: r
! 
!     type(ndarray_oc), pointer :: sol
!     call c_f_pointer(solptr, sol)
! 
!     r = sol%dim
!   end function dims

  function cast_as_ndarray_oc(encap_polymorph) result(ndarray_oc_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(ndarray_oc), pointer :: ndarray_oc_obj
    
    select type(encap_polymorph)
    type is (ndarray_oc)
       ndarray_oc_obj => encap_polymorph
    class default
      stop "TYPE ERROR in cast_as_ndarray_oc"
    end select
  end function cast_as_ndarray_oc

  
  function get_array1d_oc(x, flags) result(r)
    class(pf_encap_t), intent(in) :: x
    integer,     intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:)
    integer                :: which

    which = 0
    if (present(flags)) which = flags
!     if(.not.present(flags)) print *, "array1d_oc without flags"
    
    select type (x)
    type is (ndarray_oc)
      select case (which)
        case (1)
          r => x%yflatarray
        case (2)
          r => x%pflatarray
       case default
          stop "ERROR in get_array1d_oc: only 1, 2 allowed as flags"
       end select
    class default
       stop "ERROR: get_array1d_oc type mismatch."
    end select
  end function get_array1d_oc

  
  function get_array2d_oc(x, flags) result(r)
    class(pf_encap_t), intent(in) :: x
    integer,     intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:)
    integer                :: which

    which = 0
    if (present(flags)) which = flags
!     if(.not.present(flags)) print *, "array2d_oc without flags"
    
    select type (x)
    type is (ndarray_oc)
      select case (which)
        case (1)
          r(1:x%shape(1),1:x%shape(2)) => x%yflatarray
        case (2)
          r(1:x%shape(1),1:x%shape(2)) => x%pflatarray
        case default
          stop "ERROR in get_array2d_oc: only 1, 2 allowed as flags"
       end select
    class default
       stop "ERROR: get_array2d_oc type mismatch."
    end select
  end function get_array2d_oc
  
  function get_array3d_oc(x, flags) result(r)
    class(pf_encap_t), intent(in) :: x
    integer,     intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:)
    integer                :: which

    which = 0
    if (present(flags)) which = flags
    select type (x)
    type is (ndarray_oc)
      select case (which)
        case (1)
          r(1:x%shape(1),1:x%shape(2),1:x%shape(3)) => x%yflatarray
        case (2)
          r(1:x%shape(1),1:x%shape(2),1:x%shape(3)) => x%pflatarray
        case default
          stop "ERROR in get_array2d_oc: only 1, 2 allowed as flags"
       end select
    class default
       stop "ERROR: get_array2d_oc type mismatch."
    end select
  end function get_array3d_oc
  
  
  subroutine ndarray_oc_dump_hook(pf, lev, state)
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: lev
    type(pf_state_t),    intent(in)    :: state

    character(len=256)     :: fnamey, fnamep
    type(ndarray_oc), pointer :: qend
    
    qend => cast_as_ndarray_oc(lev%qend)
    
    write(fnamey, "('y_s',i0.2,'i',i0.3,'l',i0.2,'.npy')") &
         state%step, state%iter, lev%index

    call ndarray_dump_numpy(trim(pf%outdir)//c_null_char, trim(fnamey)//c_null_char, '<f8'//c_null_char//c_null_char, &
         qend%dim, qend%shape, size(qend%yflatarray), qend%yflatarray)

    write(fnamep, "('p_s',i0.2,'i',i0.3,'l',i0.2,'.npy')") &
         state%step, state%iter, lev%index

    call ndarray_dump_numpy(trim(pf%outdir)//c_null_char, trim(fnamep)//c_null_char, '<f8'//c_null_char//c_null_char, &
         qend%dim, qend%shape, size(qend%pflatarray), qend%pflatarray)
  end subroutine ndarray_oc_dump_hook

  
  subroutine ndarray_oc_dump_all_hook(pf, lev, state)
    type(pf_pfasst_t),   intent(inout) :: pf
    class(pf_level_t),    intent(inout) :: lev
    type(pf_state_t),    intent(in)    :: state

    character(len=256)     :: fnamey, fnamep
    integer                :: m
    
    type(ndarray_oc), pointer :: qend
   

    do m=1, lev%nnodes
      qend => cast_as_ndarray_oc(lev%Q(m))

      write(fnamey, "('y_s',i0.2,'l',i0.2,'m',i0.2,'.npy')") &
           state%step, lev%index, m

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fnamey)//c_null_char, '<f8'//c_null_char//c_null_char, &
           qend%dim,qend%shape, size(qend%yflatarray), qend%yflatarray)

      write(fnamep, "('p_s',i0.2,'l',i0.2,'m',i0.2,'.npy')") &
           state%step, lev%index, m

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fnamep)//c_null_char, '<f8'//c_null_char//c_null_char, &
           qend%dim, qend%shape, size(qend%pflatarray), qend%pflatarray)
   end do

  end subroutine ndarray_oc_dump_all_hook



  
    !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine ndarray_oc_eprint(this,flags)
    class(ndarray_oc), intent(inout) :: this
    integer,     intent(in   ), optional :: flags    
    !  Just print the first few values
    print *, this%yflatarray(1:10)
    print *, this%pflatarray(1:10)
  end subroutine ndarray_oc_eprint

!   function array2_oc(solptr, flags) result(r)
!     type(c_ptr), intent(in   ), value :: solptr
!     integer,     intent(in   ), optional :: flags
!     real(pfdp), pointer :: r(:,:)
! 
!     integer                :: shp(2)
!     integer                :: which
!     type(ndarray_oc), pointer :: sol
!     call c_f_pointer(solptr, sol)
! 
!     if (sol%dim == 2) then
!        shp = sol%shape
!        
!        which = 0
!        if (present(flags)) which = flags
! 
!        select case (which)
!        case (0)
!           stop "ERROR in array2_oc: only 1, 2 allowed as flags"
!        case (1)
!           call c_f_pointer(sol%yptr, r, shp)
!        case (2)
!           call c_f_pointer(sol%pptr, r, shp)
!        case default
!           stop "ERROR in array2_oc: only 1, 2 allowed as flags"
!        end select
!     else
!        stop "ERROR: array2 dimension mismatch."
!     end if
!   end function array2_oc
! 
!   function array3_oc(solptr, flags) result(r)
!     type(c_ptr), intent(in   ), value :: solptr
!     integer,     intent(in   ), optional :: flags
!     real(pfdp), pointer :: r(:,:,:)
! 
!     integer                :: shp(3)
!     integer                :: which
!     type(ndarray_oc), pointer :: sol
!     call c_f_pointer(solptr, sol)
! 
!     if (sol%dim == 3) then
!        shp = sol%shape
! 
!        which = 0
!        if (present(flags)) which = flags
! 
!        select case (which)
!        case (0)
!           stop "ERROR in array3_oc: only 1, 2 allowed as flags"
!        case (1)
!           call c_f_pointer(sol%yptr, r, shp)
!        case (2)
!           call c_f_pointer(sol%pptr, r, shp)
!        case default
!           stop "ERROR in array3_oc: only 1, 2 allowed as flags"
!        end select
!     else
!        stop "ERROR: array3 dimension mismatch."
!     end if
!   end function array3_oc
! 
!   function array4_oc(solptr, flags) result(r)
!     type(c_ptr), intent(in   ), value :: solptr
!     integer,     intent(in   ), optional :: flags
!     real(pfdp), pointer :: r(:,:,:,:)
! 
!     integer                :: shp(4)
!     integer                :: which
!     type(ndarray_oc), pointer :: sol
!     call c_f_pointer(solptr, sol)
! 
!     if (sol%dim == 4) then
!        shp = sol%shape
!        
!        which = 0
!        if (present(flags)) which = flags
! 
!        select case (which)
!        case (0)
!           stop "ERROR in array4_oc: only 1, 2 allowed as flags"
!        case (1)
!           call c_f_pointer(sol%yptr, r, shp)
!        case (2)
!           call c_f_pointer(sol%pptr, r, shp)
!        case default
!           stop "ERROR in array4_oc: only 1, 2 allowed as flags"
!        end select
!     else
!        stop "ERROR: array4 dimension mismatch."
!     end if
!   end function array4_oc
! 
!   function array5_oc(solptr, flags) result(r)
!     type(c_ptr), intent(in   ), value :: solptr
!     integer,     intent(in   ), optional :: flags
!     real(pfdp), pointer :: r(:,:,:,:,:)
! 
!     integer                :: shp(5)
!     integer                :: which
!     type(ndarray_oc), pointer :: sol
!     call c_f_pointer(solptr, sol)
! 
!     if (sol%dim == 5) then
!        shp = sol%shape
! 
!        which = 0
!        if (present(flags)) which = flags
! 
!        select case (which)
!        case (0)
!           stop "ERROR in array5_oc: only 1, 2 allowed as flags"
!        case (1)
!           call c_f_pointer(sol%yptr, r, shp)
!        case (2)
!           call c_f_pointer(sol%pptr, r, shp)
!        case default
!           stop "ERROR in array5_oc: only 1, 2 allowed as flags"
!        end select
!     else
!        stop "ERROR: array5 dimension mismatch."
!     end if
!   end function array5_oc
! 
!    function get_y(ptr) result(y)
!      type(c_ptr), intent(in), value :: ptr
!      real(pfdp),    pointer :: y(:)
!      type(ndarray_oc), pointer :: sol
! 
!      call c_f_pointer(ptr, sol)
!      y => sol%yflatarray
!    end function get_y
! 
!    function get_p(ptr) result(p)
!      type(c_ptr), intent(in), value :: ptr
!      real(pfdp),    pointer :: p(:)
!      type(ndarray_oc), pointer :: sol
! 
!      call c_f_pointer(ptr, sol)
!      p => sol%pflatarray
!    end function get_p
! 
!   subroutine ndarray_oc_encap_create(encap)
!     type(pf_encap_t), intent(out) :: encap
! 
!     encap%create  => ndarray_oc_create
!     encap%destroy => ndarray_oc_destroy
!     encap%setval  => ndarray_oc_setval
!     encap%copy    => ndarray_oc_copy
!     encap%norm    => ndarray_oc_norm
!     encap%pack    => ndarray_oc_pack
!     encap%unpack  => ndarray_oc_unpack
!     encap%axpy    => ndarray_oc_saxpy
!   end subroutine ndarray_oc_encap_create

end module pf_mod_ndarray_oc
