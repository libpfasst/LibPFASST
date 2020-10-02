!!  LibPFASST array encapsulation using amrex multifab.
!
! This file is part of LIBPFASST.
!

!> Module to define and encapsulation of a Petsc Vector.
!!
!! When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'shape_lev' attribute to create a new array with that
!! shape.  Thus, the 'shape'_lev attributes of the PFASST levels should be
!! set appropriately.

module pf_mod_AMReX_mfab
  use amrex_base_module
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  !>  Type to create and destroy N-dimenstional arrays 
  type, extends(pf_factory_t) :: pf_AMReX_mfab_factory_t
   contains
     procedure :: create_single  => AMReX_mfab_create_single
     procedure :: create_array  => AMReX_mfab_create_array
     procedure :: destroy_single => AMReX_mfab_destroy_single
     procedure :: destroy_array => AMReX_mfab_destroy_array
  end type pf_AMReX_mfab_factory_t
  
  !>  1-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: pf_amrex_mfab_t
     integer             :: ndim
     integer   :: arr_shape(4)
     integer   :: pack_size(4)
     integer :: n_cell, max_grid_size, nsteps, plot_int
     integer :: ncomp     ! Number of componentns
     integer ::  nghost   ! number of  ghost cells
     integer :: istep
     type(amrex_parmparse) :: pp
     type(amrex_box) :: domain
     type(amrex_geometry)  :: geom
     type(amrex_multifab)  :: mfab
     
   contains
     procedure :: setval => AMReX_mfab_setval
     procedure :: copy => AMReX_mfab_copy
     procedure :: norm => AMReX_mfab_norm
     procedure :: pack => AMReX_mfab_pack
     procedure :: unpack => AMReX_mfab_unpack
     procedure :: axpy => AMReX_mfab_axpy
     procedure :: eprint => AMReX_mfab_eprint
  end type pf_amrex_mfab_t
  
  

contains
  function cast_as_AMReX_mfab(encap_polymorph) result(pf_AMReX_mfab_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pf_AMReX_mfab_t), pointer :: pf_AMReX_mfab_obj
    
    select type(encap_polymorph)
    type is (pf_AMReX_mfab_t)
       pf_AMReX_mfab_obj => encap_polymorph
    end select
  end function cast_as_AMReX_mfab

  !>  Subroutine to allocate the array and set the size parameters
  subroutine AMReX_mfab_build(this, shape_in)
    use pf_mod_comm_mpi
    class(pf_encap_t), intent(inout) :: this
    integer,           intent(in   ) :: shape_in(:)

    integer nn,psize,rank,ierr,k
    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_parmparse) :: pp
    integer :: n_cell, max_grid_size

    
    select type (this)
    class is (pf_AMReX_mfab_t)
       this%ndim   = SIZE(shape_in)-1
       this%ncomp=shape_in(1)
       this%nghost=0
       this%arr_shape = shape_in(2:this%ndim+1)
       print *,'ndim=',this%ndim,' ncomp=',this%ncomp,' nghost=',this%nghost, ' array_shape=',this%arr_shape(1:this%ndim)
       
       ! Define a single box covering the domain
       this%domain = amrex_box((/0,0,0/), (/shape_in(2),shape_in(3),shape_in(4)/))

       ! Initialize the boxarray "ba" from the single box "bx"
       call amrex_boxarray_build(ba, this%domain)

       this%max_grid_size=64
       ! Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
       call ba%maxSize(this%max_grid_size)
       
       ! Build a DistributionMapping for the boxarray

       call amrex_distromap_build(dm, ba)

       ! This defines a amrex_geometry object.
!       print *,'set geometry coord'
       call amrex_geometry_set_coord_sys(0)
       call amrex_geometry_set_prob_domain((/0.0d0,0.0d0,0.0d0/), (/1.0d0,1.0d0,1.0d0/))
       call amrex_geometry_set_periodic ((/ .true.,.true.,.true./))
       call amrex_geometry_build(this%geom, this%domain)
       
       ! Build data multifabs
       call amrex_multifab_build(this%mfab, ba, dm, this%ncomp, this%nghost)
       
       call amrex_distromap_destroy(dm)
       call amrex_boxarray_destroy(ba)
       
       !  Make a shape the size of the grid with ghost cells
       this%pack_size(1:size(shape_in)) = shape_in
       do k=2,this%ndim+1
          if (this%pack_size(k)>1) this%pack_size(k)=this%pack_size(k)+2*this%nghost
       end do
       
       
    end select
  end subroutine AMReX_mfab_build

  !> Subroutine to  create a single array
  subroutine AMReX_mfab_create_single(this, x, level_index, lev_shape)
    class(pf_AMReX_mfab_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: ierr
    allocate(pf_AMReX_mfab_t::x,stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    call AMReX_mfab_build(x, lev_shape)
  end subroutine AMReX_mfab_create_single

  !> Subroutine to create an array of arrays
  subroutine AMReX_mfab_create_array(this, x, n, level_index,  lev_shape)
    class(pf_AMReX_mfab_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: i,ierr
    allocate(pf_AMReX_mfab_t::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    do i = 1, n
       call AMReX_mfab_build(x(i), lev_shape)
    end do
  end subroutine AMReX_mfab_create_array

  !>  Subroutine to destroy array
  subroutine AMReX_mfab_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_AMReX_mfab_t), pointer :: pf_AMReX_mfab_obj
    integer ::  ierr
    pf_AMReX_mfab_obj => cast_as_AMReX_mfab(encap)
!    print *,'destroying AMReX_mfab'

    call amrex_multifab_destroy(pf_AMReX_mfab_obj%mfab)
    call amrex_geometry_destroy(pf_AMReX_mfab_obj%geom) 
    nullify(pf_AMReX_mfab_obj)


  end subroutine AMReX_mfab_destroy

  !> Subroutine to destroy an single array
  subroutine AMReX_mfab_destroy_single(this, x)
    class(pf_AMReX_mfab_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer ::  ierr

    select type (x)
    class is (pf_AMReX_mfab_t)
       call amrex_multifab_destroy(x%mfab)
       call amrex_geometry_destroy(x%geom)
    end select
    deallocate(x)
  end subroutine AMReX_mfab_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine AMReX_mfab_destroy_array(this, x)
    class(pf_AMReX_mfab_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout),allocatable :: x(:)
    integer                                            :: i,ierr
    select type(x)
    class is (pf_AMReX_mfab_t)
       do i = 1,SIZE(x)
          call amrex_multifab_destroy(x(i)%mfab)
          call amrex_geometry_destroy(x(i)%geom)
       end do
    end select
    deallocate(x)
  end subroutine AMReX_mfab_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalar  value.
  subroutine AMReX_mfab_setval(this, val, flags)
    class(pf_AMReX_mfab_t), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags

    call this%mfab%setval(val,1,this%ncomp,this%nghost)       
  end subroutine AMReX_mfab_setval

  !> Subroutine to copy an array
  subroutine AMReX_mfab_copy(this, src, flags)
    class(pf_AMReX_mfab_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags

    integer ng,nc ! for debug

    
    select type(src)
    type is (pf_AMReX_mfab_t)
       ng=this%nghost
       nc=this%ncomp
       !call this%mfab%amrex_multifab_copy(src%mfab,1,1,nc,ng)
       !       call this%mfab%copy(src%mfab,1,1,nc,ng)
       call this%mfab%parallel_copy(src%mfab,src%geom)              
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine AMReX_mfab_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine AMReX_mfab_pack(this, z, flags)
    class(pf_AMReX_mfab_t), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: psize
    real(pfdp),  pointer :: mfab_data(:,:,:,:)

    mfab_data=>this%mfab%dataPtr(0)

    z=reshape(mfab_data,[product(this%pack_size)])

    
  end subroutine AMReX_mfab_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine AMReX_mfab_unpack(this, z, flags)
    class(pf_AMReX_mfab_t), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags

    real(pfdp),  pointer :: mfab_data(:,:,:,:)

    mfab_data=>this%mfab%dataPtr(0)

    mfab_data=reshape(z,this%pack_size)

    
  end subroutine AMReX_mfab_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function AMReX_mfab_norm(this, flags) result (norm)
    class(pf_AMReX_mfab_t), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    norm=1.0
!    if (present(flags)) then
       norm = this%mfab%norm2()
!    else
!       norm = this%mfab%norm2(0)
!    end if
    
!    call VecNorm(this%AMReX_mfab,NORM_INFINITY,norm,this%ierr);CHKERRQ(this%ierr)
  end function AMReX_mfab_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine AMReX_mfab_axpy(this, a, x, flags)
    class(pf_AMReX_mfab_t),    intent(inout)       :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags

    select type(x)
    type is (pf_AMReX_mfab_t)
       call this%mfab%saxpy(a,x%mfab,1,1,this%ncomp,this%nghost)       
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine AMReX_mfab_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine AMReX_mfab_eprint(this,flags)
    class(pf_AMReX_mfab_t), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    real(pfdp),  pointer :: mfab_data(:,:,:,:)

    mfab_data=>this%mfab%dataPtr(0)

 
    !  Just print the first few values
    print *,mfab_data(1,1:10,1,1)

  end subroutine AMReX_mfab_eprint


end module pf_mod_AMReX_mfab
