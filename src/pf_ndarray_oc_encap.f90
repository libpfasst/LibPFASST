!! N-dimensional array encapsulation for optimal control problems.
!
! This file is part of LIBPFASST.
!
!!
!> When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'lev_shape' attribute to create a new array with that
!! shape.  Thus, the 'lev_shape' attributes of the PFASST levels should be
!! set appropriately.  For example, before calling pf_pfasst_run we can
!! set the shape of the coarsest level by doing:
!!
!!   allocate(pf%levels(1)%lev_shape(2))
!!   pf%levels(1)%lev_shape = [ 3, 10 ]
!!
!! The helper routines get_array1d_oc, get_array2d_oc, get_array3d_oc, etc can be used to
!! extract pointers to the encapsulated array from a C pointer without
!! performing any copies.
module pf_mod_ndarray_oc
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_utils
!  use fnpy
  implicit none

  !>  Type to create and destroy N-dimenstional arrays for optimal control
  type, extends(pf_factory_t) :: pf_ndarray_oc_factory_t
   contains
     procedure :: create_single  => ndarray_oc_create_single
     procedure :: create_array   => ndarray_oc_create_array
     procedure :: destroy_single => ndarray_oc_destroy_single
     procedure :: destroy_array  => ndarray_oc_destroy_array
  end type pf_ndarray_oc_factory_t

  !>  N-dimensional array type for optimal control,  extends the abstract encap type  
  type, extends(pf_encap_t) :: pf_ndarray_oc_t
     integer                 :: ndim
     integer,    allocatable :: arr_shape(:)    
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
  end type pf_ndarray_oc_t

  ! interfaces to routines in pf_numpy.c
  interface
     subroutine ndarray_dump_numpy(dname, fname, endian, ndim, shape_in, nvars, array) bind(c)
       use iso_c_binding
       use pf_mod_dtype
       character(c_char), intent(in   )        :: dname, fname, endian(5)
       integer(c_int),    intent(in   ), value :: ndim, nvars
       integer(c_int),    intent(in   )        :: shape_in(ndim)
       real(pfdp),    intent(in   )        :: array(nvars)
     end subroutine ndarray_dump_numpy
  end interface

contains
  !>  Subroutine to allocate the array and set the size parameters
  subroutine ndarray_oc_build(q, shape_in)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape_in(:)
    integer :: ierr
    select type (q)
    class is (pf_ndarray_oc_t)
       allocate(q%arr_shape(SIZE(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
       allocate(q%yflatarray(product(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
       allocate(q%pflatarray(product(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
       q%ndim   = SIZE(shape_in)
       q%arr_shape = shape_in
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine ndarray_oc_build

  !> Subroutine to  create a single array
  subroutine ndarray_oc_create_single(this, x, level_index, lev_shape)
    class(pf_ndarray_oc_factory_t), intent(inout)           :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)

    integer :: ierr
    allocate(pf_ndarray_oc_t::x,stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
    
    call ndarray_oc_build(x, lev_shape)
  end subroutine ndarray_oc_create_single

  !> Subroutine to create an array of arrays
  subroutine ndarray_oc_create_array(this, x, n, level_index, lev_shape)
    class(pf_ndarray_oc_factory_t), intent(inout)              :: this
    class(pf_encap_t),         intent(inout), allocatable :: x(:)
    integer,                   intent(in   )              :: n
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: i,ierr
    allocate(pf_ndarray_oc_t::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
    
    do i = 1, n
       call ndarray_oc_build(x(i), lev_shape)
    end do
  end subroutine ndarray_oc_create_array
  
  !>  Subroutine to destroy array
  subroutine ndarray_oc_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_ndarray_oc_t), pointer :: ndarray_oc_obj

    ndarray_oc_obj => cast_as_ndarray_oc(encap) !??

    deallocate(ndarray_oc_obj%pflatarray)
    deallocate(ndarray_oc_obj%yflatarray)
    deallocate(ndarray_oc_obj%arr_shape)
    nullify(ndarray_oc_obj)
  end subroutine ndarray_oc_destroy

    !> Subroutine to destroy an single array
  subroutine ndarray_oc_destroy_single(this, x)
    class(pf_ndarray_oc_factory_t), intent(inout)              :: this
    class(pf_encap_t),         intent(inout), allocatable :: x
    
    select type (x)
    class is (pf_ndarray_oc_t)
       deallocate(x%pflatarray)
       deallocate(x%yflatarray)
       deallocate(x%arr_shape)
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
    deallocate(x)
  end subroutine ndarray_oc_destroy_single
  
    !> Subroutine to destroy an array of arrays
  subroutine ndarray_oc_destroy_array(this, x)
    class(pf_ndarray_oc_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer                                            :: i
    
    select type(x)
    class is (pf_ndarray_oc_t)
       do i = 1,SIZE(x)
          deallocate(x(i)%pflatarray)
          deallocate(x(i)%yflatarray)
          deallocate(x(i)%arr_shape)
       end do
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
    deallocate(x)
  end subroutine ndarray_oc_destroy_array
  
  
  !> Subroutine to set array to a scalar  value.
  subroutine ndarray_oc_setval(this, val, flags)
    class(pf_ndarray_oc_t), intent(inout)     :: this
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
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Select case error',which)
    end select

  end subroutine ndarray_oc_setval

  !> Subroutine to copy an array
  subroutine ndarray_oc_copy(this, src, flags)
    class(pf_ndarray_oc_t),    intent(inout)  :: this
    class(pf_encap_t),    intent(in   )  :: src
    integer,     intent(in   ), optional :: flags
    integer :: which

    which = 0
    if (present(flags)) which = flags
!     if(.not.present(flags)) print *, "copy without flags"
    
    select type(src)
    type is (pf_ndarray_oc_t)
      select case (which)
      case (0)
        this%yflatarray = src%yflatarray
        this%pflatarray = src%pflatarray
      case (1)
        this%yflatarray = src%yflatarray
      case (2)
        this%pflatarray = src%pflatarray
      case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',which)
      end select
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine ndarray_oc_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine ndarray_oc_pack(this, z, flags)
    class(pf_ndarray_oc_t), intent(in   ) :: this
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
       call pf_stop(__FILE__,__LINE__,'ndarray_oc_pack: only 1, 2 allowed as flags, which=',which)
    case (1)
       z = this%yflatarray
    case (2)
       z = this%pflatarray
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',which)
    end select
  end subroutine ndarray_oc_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine ndarray_oc_unpack(this, z, flags)
    class(pf_ndarray_oc_t), intent(inout) :: this
    real(pfdp),  intent(in   )        :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: which

    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (0)
       call pf_stop(__FILE__,__LINE__,'ndarray_oc_unpack: only 1, 2 allowed as flags, which=',which)
    case (1)
       this%yflatarray = z
    case (2)
       this%pflatarray = z 
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',which)
    end select
  end subroutine ndarray_oc_unpack
  

  !> Subroutine to define the norm of the array (here the max norm)
  function ndarray_oc_norm(this, flags) result (norm)
    class(pf_ndarray_oc_t), intent(in   ) :: this
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
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',which)
    end select
  end function ndarray_oc_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine ndarray_oc_axpy(this, a, x, flags)
    class(pf_ndarray_oc_t), intent(inout) :: this    
    class(pf_encap_t), intent(in   )     :: x
    real(pfdp),  intent(in   )           :: a
    integer,     intent(in   ), optional :: flags
    integer :: which

    if (a .eq. 0.0_pfdp) return
    which = 0
    if (present(flags)) which = flags
!     if (.not.present(flags)) stop "axpy without flags" 
    
    select type(x)
    type is (pf_ndarray_oc_t)
      select case (which)
      case (0)
        this%yflatarray = a * x%yflatarray + this%yflatarray
        this%pflatarray = a * x%pflatarray + this%pflatarray
      case (1)
        this%yflatarray = a * x%yflatarray + this%yflatarray
      case (2)
        this%pflatarray = a * x%pflatarray + this%pflatarray
      case DEFAULT
         call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',which)
      end select
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select  
  end subroutine ndarray_oc_axpy

  

  function cast_as_ndarray_oc(encap_polymorph) result(ndarray_oc_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pf_ndarray_oc_t), pointer :: ndarray_oc_obj
    
    select type(encap_polymorph)
    type is (pf_ndarray_oc_t)
       ndarray_oc_obj => encap_polymorph
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end function cast_as_ndarray_oc

  
  function get_array1d_oc(x, flags) result(r)
    class(pf_encap_t), target, intent(in) :: x
    integer,     intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:)
    integer                :: which

    which = 0
    if (present(flags)) which = flags
!     if(.not.present(flags)) print *, "array1d_oc without flags"
    
    select type (x)
    type is (pf_ndarray_oc_t)
      select case (which)
        case (1)
          r => x%yflatarray
        case (2)
          r => x%pflatarray
       case DEFAULT
          call pf_stop(__FILE__,__LINE__,'gerarray1d: only 1, 2 allowed as flags, which=',which)
       end select
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end function get_array1d_oc

  
  function get_array2d_oc(x, flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,     intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:)
    integer                :: which

    which = 0
    if (present(flags)) which = flags
!     if(.not.present(flags)) print *, "array2d_oc without flags"
    
    select type (x)
    type is (pf_ndarray_oc_t)
      select case (which)
        case (1)
          r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%yflatarray
        case (2)
          r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%pflatarray
       case DEFAULT
          call pf_stop(__FILE__,__LINE__,'gerarray2d: only 1, 2 allowed as flags, which=',which)          
       end select
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end function get_array2d_oc
  
  function get_array3d_oc(x, flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,     intent(in   ), optional :: flags
    real(pfdp), pointer :: r(:,:,:)
    integer                :: which

    which = 0
    if (present(flags)) which = flags
    select type (x)
    type is (pf_ndarray_oc_t)
      select case (which)
        case (1)
          r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%yflatarray
        case (2)
          r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%pflatarray
        case DEFAULT
          call pf_stop(__FILE__,__LINE__,'gerarray3d: only 1, 2 allowed as flags, which=',which)          
       end select
    class DEFAULT
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end function get_array3d_oc
  
  
  subroutine ndarray_oc_dump_hook(pf, level_index)
    type(pf_pfasst_t),   intent(inout) :: pf
    integer,             intent(in)    :: level_index

    character(len=256)     :: fnamey, fnamep
    type(pf_ndarray_oc_t), pointer :: qend

    qend => cast_as_ndarray_oc(pf%levels(level_index)%qend)
    
    write(fnamey, "('y_s',i0.2,'i',i0.3,'l',i0.2,'.npy')") &
         pf%state%step, pf%state%iter, level_index

!    call ndarray_dump_numpy(trim(pf%outdir)//c_null_char, trim(fnamey)//c_null_char, '<f8'//c_null_char//c_null_char, &
!         qend%ndim, qend%arr_shape, SIZE(qend%yflatarray), qend%yflatarray)


    write(fnamep, "('p_s',i0.2,'i',i0.3,'l',i0.2,'.npy')") &
         pf%state%step, pf%state%iter, level_index


!    call ndarray_dump_numpy(trim(pf%outdir)//c_null_char, trim(fnamep)//c_null_char, '<f8'//c_null_char//c_null_char, &
!         qend%ndim, qend%arr_shape, SIZE(qend%pflatarray), qend%pflatarray)
  end subroutine ndarray_oc_dump_hook

  
  subroutine ndarray_oc_dump_all_hook(pf, level_index)
    type(pf_pfasst_t),   intent(inout) :: pf
    integer,             intent(in)    :: level_index

    character(len=256)     :: fnamey, fnamep
    integer                :: m
    
    type(pf_ndarray_oc_t), pointer :: qend

    do m=1, pf%levels(level_index)%nnodes
      qend => cast_as_ndarray_oc(pf%levels(level_index)%Q(m))

      write(fnamey, "('y_s',i0.2,'l',i0.2,'m',i0.2,'.npy')") &
           pf%state%step, level_index, m



      write(fnamep, "('p_s',i0.2,'l',i0.2,'m',i0.2,'.npy')") &
           pf%state%step, level_index, m

   end do

  end subroutine ndarray_oc_dump_all_hook

  
    !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine ndarray_oc_eprint(this,flags)
    class(pf_ndarray_oc_t), intent(inout) :: this
    integer,     intent(in   ), optional :: flags    
    !  Just print the first few values
    print *, this%yflatarray(1:10)
    print *, this%pflatarray(1:10)
  end subroutine ndarray_oc_eprint

end module pf_mod_ndarray_oc
