!!  N-dimensional array encapsulation using petsc.
!
! This file is part of LIBPFASST.
!

!> Module to define and encapsulation of a Petsc Vector.
!!
!! When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'shape_lev' attribute to create a new array with that
!! shape.  Thus, the 'shape'_lev attributes of the PFASST levels should be
!! set appropriately.  For example, before calling pf_pfasst_run we can
!! set the shape of the coarsest level by doing:
!!
!!   allocate(pf%levels_lev(1)%shape(2))
!!   pf%levels(1)%shape_lev = [ 3, 10 ]
!!
module pf_mod_petscVec
#include <petsc/finclude/petscvec.h>
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_utils
  use petscvec
  implicit none

  !>  Type to create and destroy N-dimenstional arrays 
  type, extends(pf_factory_t) :: pf_petscVec_factory_t
   contains
     procedure :: create_single  => pf_petscVec_create_single
     procedure :: create_array  => pf_petscVec_create_array
     procedure :: destroy_single => pf_petscVec_destroy_single
     procedure :: destroy_array => pf_petscVec_destroy_array
  end type pf_petscVec_factory_t
  
  !>  1-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: pf_petscVec_t
     integer             :: ndim
     integer,    allocatable :: arr_shape(:)
     type(tVec) ::  petscVec
     integer :: ierr
   contains
     procedure :: setval => pf_petscVec_setval
     procedure :: copy => pf_petscVec_copy
     procedure :: norm => pf_petscVec_norm
     procedure :: pack => pf_petscVec_pack
     procedure :: unpack => pf_petscVec_unpack
     procedure :: axpy => pf_petscVec_axpy
     procedure :: eprint => pf_petscVec_eprint
  end type pf_petscVec_t

contains
  function cast_as_pf_petscVec(encap_polymorph) result(pf_petscVec_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pf_petscVec_t), pointer :: pf_petscVec_obj
    
    select type(encap_polymorph)
    type is (pf_petscVec_t)
       pf_petscVec_obj => encap_polymorph
    end select
  end function cast_as_pf_petscVec

  !>  Subroutine to allocate the array and set the size parameters
  subroutine pf_petscVec_build(q, shape_in)
    use pf_mod_comm_mpi
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: shape_in(:)

    integer nn,psize,rank,ierr
    select type (q)
    class is (pf_petscVec_t)

       call VecCreate(PETSC_COMM_WORLD,q%petscVec,ierr);CHKERRQ(ierr)
       call VecSetSizes(q%petscVec,PETSC_DECIDE,shape_in(1),ierr);CHKERRQ(ierr)
       
       call VecSetFromOptions(q%petscVec,ierr);CHKERRQ(ierr)
       call VecGetLocalSize(q%petscVec,psize,ierr);CHKERRQ(ierr)
       call mpi_comm_rank(PETSC_COMM_WORLD, rank,ierr);CHKERRQ(ierr)

       allocate(q%arr_shape(SIZE(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                                
       q%ndim   = SIZE(shape_in)
       q%arr_shape = shape_in
    end select
  end subroutine pf_petscVec_build

  !> Subroutine to  create a single array
  subroutine pf_petscVec_create_single(this, x, level_index, lev_shape)
    class(pf_petscVec_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: ierr
    allocate(pf_petscVec_t::x,stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    call pf_petscVec_build(x, lev_shape)
  end subroutine pf_petscVec_create_single

  !> Subroutine to create an array of arrays
  subroutine pf_petscVec_create_array(this, x, n, level_index,  lev_shape)
    class(pf_petscVec_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: i,ierr
    allocate(pf_petscVec_t::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    do i = 1, n
       call pf_petscVec_build(x(i), lev_shape)
    end do
  end subroutine pf_petscVec_create_array

  !>  Subroutine to destroy array
  subroutine pf_petscVec_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_petscVec_t), pointer :: pf_petscVec_obj
    integer ::  ierr
    pf_petscVec_obj => cast_as_pf_petscVec(encap)
!    print *,'destroying petscVec'
    call VecDestroy(pf_petscVec_obj%petscVec,ierr);CHKERRQ(ierr)
    deallocate(pf_petscVec_obj%arr_shape)

    nullify(pf_petscVec_obj)

  end subroutine pf_petscVec_destroy

  !> Subroutine to destroy an single array
  subroutine pf_petscVec_destroy_single(this, x)
    class(pf_petscVec_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer ::  ierr

    select type (x)
    class is (pf_petscVec_t)
       call VecDestroy(x%petscVec,ierr);CHKERRQ(ierr)
       deallocate(x%arr_shape)
    end select
    deallocate(x)
  end subroutine pf_petscVec_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine pf_petscVec_destroy_array(this, x)
    class(pf_petscVec_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout),allocatable :: x(:)
    integer                                            :: i,ierr
    select type(x)
    class is (pf_petscVec_t)
       do i = 1,SIZE(x)
          call VecDestroy(x(i)%petscVec,ierr);CHKERRQ(ierr)
          deallocate(x(i)%arr_shape)
       end do
    end select
    deallocate(x)
  end subroutine pf_petscVec_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalar  value.
  subroutine pf_petscVec_setval(this, val, flags)
    class(pf_petscVec_t), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    call VecSet(this%petscVec,val,this%ierr);CHKERRQ(this%ierr)
  end subroutine pf_petscVec_setval

  !> Subroutine to copy an array
  subroutine pf_petscVec_copy(this, src, flags)
    class(pf_petscVec_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags

    integer m,n ! for debug

    select type(src)
    type is (pf_petscVec_t)

        call VecGetSize(src%petscVec,m, this%ierr);CHKERRQ(this%ierr)
        call VecGetSize(this%petscVec,n, this%ierr);CHKERRQ(this%ierr)

        call VecCopy(src%petscVec,this%petscVec,this%ierr);CHKERRQ(this%ierr)

    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine pf_petscVec_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine pf_petscVec_pack(this, z, flags)
    class(pf_petscVec_t), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags

    real(pfdp),  pointer :: p_petscVec(:)
    integer :: psize

    call VecGetArrayReadF90(this%petscVec,p_petscVec,this%ierr);CHKERRQ(this%ierr)
    call VecGetLocalSize(this%petscVec,psize,this%ierr);CHKERRQ(this%ierr)
    z=p_petscVec(1:psize)
    
    call VecRestoreArrayReadF90(this%petscVec,p_petscVec,this%ierr);CHKERRQ(this%ierr)
    
  end subroutine pf_petscVec_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine pf_petscVec_unpack(this, z, flags)
    class(pf_petscVec_t), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags

    real(pfdp),  pointer :: p_petscVec(:)
    integer :: psize

    call VecGetArrayF90(this%petscVec,p_petscVec,this%ierr);CHKERRQ(this%ierr)
    call VecGetLocalSize(this%petscVec,psize,this%ierr);CHKERRQ(this%ierr)

    p_petscVec=z(1:psize)
    call VecRestoreArrayF90(this%petscVec,p_petscVec,this%ierr);CHKERRQ(this%ierr)
    
  end subroutine pf_petscVec_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function pf_petscVec_norm(this, flags) result (norm)
    class(pf_petscVec_t), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm

    call VecNorm(this%petscVec,NORM_INFINITY,norm,this%ierr);CHKERRQ(this%ierr)
  end function pf_petscVec_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine pf_petscVec_axpy(this, a, x, flags)
    class(pf_petscVec_t),    intent(inout)       :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags

    select type(x)
    type is (pf_petscVec_t)


       call VecAXPY(this%petscVec,a,x%petscVec,this%ierr);CHKERRQ(this%ierr)

    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine pf_petscVec_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine pf_petscVec_eprint(this,flags)
    class(pf_petscVec_t), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
    if (product(this%arr_shape) < 10) then

    else

    endif
  end subroutine pf_petscVec_eprint


end module pf_mod_petscVec
