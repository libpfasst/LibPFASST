!!  Hamiltonian system  encapsulation.
!
! This file is part of LIBPFASST.
!
!> Hamiltonian system encapsulation.
!!
!! This encapsulation uses the levels 'lev_shape' attribute to create define 2 new arrays
!! of size (ndim,nparticles) corresponding to the two components (p,q) in the Hamiltonian particles system
!!
!! The helper routines get_array can be used to
!! extract pointers to the different components of the system array 
!! performing any copies.
!!
module pf_mod_ham_sys
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  !>  Type to create and destroy N-dimenstional arrays
  type, extends(pf_factory_t) :: pf_ham_sys_factory
   contains
     procedure :: create_single  => ham_sys_create_single
     procedure :: create_array  => ham_sys_create_array
     procedure :: destroy_single => ham_sys_destroy_single
     procedure :: destroy_array => ham_sys_destroy_array
  end type pf_ham_sys_factory

  !>  N-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: pf_ham_sys_t
     integer             :: ndim        !  spatial dimension
     integer             :: nparticles  !  number of grid points or grid points per dimension
     integer             :: ndof   !  ndim*nparticles
     real(pfdp), pointer :: p(:,:),q(:,:) !  pointers to (p,q) solution     
     real(pfdp), allocatable :: flatarray(:)  !  of size 2*ndim*nparticles or 2*ndof
   contains
     procedure :: setval => ham_sys_setval
     procedure :: copy => ham_sys_copy
     procedure :: norm => ham_sys_norm
     procedure :: pack => ham_sys_pack
     procedure :: unpack => ham_sys_unpack
     procedure :: axpy => ham_sys_axpy
     procedure :: eprint => ham_sys_eprint
  end type pf_ham_sys_t

  !> Interfaces to output routines in pf_numpy.c
  interface
     !>  Subroutine to make a directory for output
     subroutine ham_sys_mkdir(dname, dlen) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname
       integer,    intent(in   ), value :: dlen
     end subroutine ham_sys_mkdir
  end interface

contains
  !>  Subroutine to allocate the array and set the size parameters
  subroutine ham_sys_build(q, lev_shape)
    class(pf_encap_t), intent(inout),target :: q
    integer,           intent(in   ) :: lev_shape(2)

    select type (q)
    class is (pf_ham_sys_t)
       q%nparticles  =lev_shape(2)
       q%ndim   = lev_shape(1)
       q%ndof   = q%nparticles*q%ndim

       allocate(q%flatarray(2*q%ndof))
       q%flatarray=0.0_pfdp
       
       q%p(1:q%ndim,1:q%nparticles) => q%flatarray(1:q%ndof)
       q%q(1:q%ndim,1:q%nparticles) => q%flatarray(q%ndof+1:2*q%ndof)              
    end select
  end subroutine ham_sys_build

  !> Subroutine to  create a single array
  subroutine ham_sys_create_single(this, x, level_index, lev_shape)
    class(pf_ham_sys_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index, lev_shape(:)
    integer :: i
    allocate(pf_ham_sys_t::x)
    call ham_sys_build(x,lev_shape)
  end subroutine ham_sys_create_single

  !> Subroutine to create an array of arrays
  subroutine ham_sys_create_array(this, x, n, level_index, lev_shape)
    class(pf_ham_sys_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n, level_index, lev_shape(:)
    integer :: i
    allocate(pf_ham_sys_t::x(n))
    do i = 1, n
       call ham_sys_build(x(i), lev_shape)
    end do
  end subroutine ham_sys_create_array

  !>  Subroutine to destroy array
  subroutine ham_sys_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_ham_sys_t), pointer :: ham_sys_obj

    ham_sys_obj => cast_as_ham_sys(encap)

    deallocate(ham_sys_obj%flatarray)
    nullify(ham_sys_obj)

  end subroutine ham_sys_destroy

  !> Subroutine to destroy an single array
  subroutine ham_sys_destroy_single(this, x)
    class(pf_ham_sys_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x

    select type (x)
    class is (pf_ham_sys_t)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine ham_sys_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine ham_sys_destroy_array(this, x)
    class(pf_ham_sys_factory), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer                                            :: i

    select type(x)
    class is (pf_ham_sys_t)
       do i = 1,SIZE(x)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
  end subroutine ham_sys_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  !!  The integer flag "flags" determines which value is set
  !!  Usually it is either 0,1,2  corresponding to both, 1st, or 2nd
  !> Subroutine to set array to a scalar  value.
  subroutine ham_sys_setval(this, val, flags)
    class(pf_ham_sys_t), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags

    integer :: icomp
    if (present(flags)) then
       icomp = flags
    else
       icomp = 0
    end if
    
    select case (icomp)
    case(0)
       this%flatarray = val
    case(1)
       this%p =  val
    case(2)
       this%q =  val
    case default
       stop "bad flags"
    end select
    
  end subroutine ham_sys_setval

  !> Subroutine to copy an array
  subroutine ham_sys_copy(this, src, flags)
    class(pf_ham_sys_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags

    integer :: icomp

    if (present(flags)) then
       icomp = flags
    else
       icomp = 0
    end if

    select type(src)
    type is (pf_ham_sys_t)
       select case (icomp)
       case(0)
          this%flatarray = src%flatarray
       case(1)
          this%p =  src%p
       case(2)
          this%q =  src%q          
       case default
          stop "bad flags"
       end select
    class default
       stop "TYPE ERROR"
    end select

  end subroutine ham_sys_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine ham_sys_pack(this, z, flags)
    class(pf_ham_sys_t), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    z = this%flatarray
  end subroutine ham_sys_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine ham_sys_unpack(this, z, flags)
    class(pf_ham_sys_t), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags
    this%flatarray = z
  end subroutine ham_sys_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function ham_sys_norm(this, flags) result (norm)
    class(pf_ham_sys_t), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm = maxval(abs(this%flatarray))
  end function ham_sys_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine ham_sys_axpy(this, a, x, flags)
    class(pf_ham_sys_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags

    integer :: icomp

     if (present(flags)) then
       icomp = flags
    else
       icomp = 0
    end if
    select type (x)
    type is (pf_ham_sys_t)
       select case (icomp)
       case(0)
          this%flatarray = a * x%flatarray + this%flatarray
       case(1)
          this%p = a * x%p + this%p
       case(2)
          this%q = a * x%q + this%q          
       case(12)
          this%q = a * x%p + this%q
       case default
          stop "bad flags"
       end select
    class default       
       stop "TYPE ERROR"
    end select
  end subroutine ham_sys_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine ham_sys_eprint(this,flags)
    class(pf_ham_sys_t), intent(inout) :: this
    integer, optional, intent(in) :: flags

    integer :: n

    !  Just print the first few values
    n=min(10,this%nparticles)
    print *,'p=', this%p(1:this%ndim,1:n)
    print *,'q=', this%q(1:this%ndim,1:n)    

  end subroutine ham_sys_eprint

  function cast_as_ham_sys(encap_polymorph) result(ham_sys_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pf_ham_sys_t), pointer :: ham_sys_obj
    
    select type(encap_polymorph)
    type is (pf_ham_sys_t)
       ham_sys_obj => encap_polymorph
    end select
  end function cast_as_ham_sys

end module pf_mod_ham_sys
