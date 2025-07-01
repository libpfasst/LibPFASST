!!  LibPFASST array encapsulation using amrex multifab.
!
! This file is part of LIBPFASST.
!
!!
!! When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'shape_lev' attribute to create a new array with that
!! shape.  Thus, the 'shape'_lev attributes of the PFASST levels should be
!! set appropriately.
module pf_mod_AMReX_mfab
  use amrex_base_module
  use amrex_fort_module
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  !>  Type to create and destroy N-dimenstional arrays 
  type, extends(pf_factory_t) :: pf_AMReX_mfab_factory_t
     type(amrex_boxarray),  pointer :: ba   => null()
     type(amrex_distromap), pointer :: dm   => null()
     type(amrex_geometry),  pointer :: geom => null()
   contains
     procedure :: create_single  => AMReX_mfab_create_single
     procedure :: create_array  => AMReX_mfab_create_array
     procedure :: destroy_single => AMReX_mfab_destroy_single
     procedure :: destroy_array => AMReX_mfab_destroy_array
  end type pf_AMReX_mfab_factory_t
  
  !>  AMReX array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: pf_amrex_mfab_t
     integer                        :: ndim                 ! Number of spatial dimensions
     integer                        :: ncomp                ! Number of solution components
     integer                        :: nghost               ! Number of ghost cells
     integer                        :: ndof                 ! Number of DOFs
     integer                        :: max_grid_size        ! AMReX BoxArray Parameter
     integer                        :: arr_shape(4)         ! (nx,ny,nz,ncomp)
     integer                        :: pack_size(4)         ! arr_shape plus ghost cells (nx+2*nghost,ny+2*nghost,nz+2*nghost,ncomp)
     integer, allocatable           :: am_to_flat(:,:,:,:)  ! index for the array
     type(amrex_geometry)           :: geom                 ! pointer to global geometry - needed for fill_boundary in AMReX mfab
     type(amrex_multifab)  :: mfab  !  The actual multifab
     
   contains
     procedure :: setval => AMReX_mfab_setval
     procedure :: copy => AMReX_mfab_copy
     procedure :: norm => AMReX_mfab_norm
     procedure :: pack => AMReX_mfab_pack
     procedure :: unpack => AMReX_mfab_unpack
     procedure :: axpy => AMReX_mfab_axpy
     procedure :: eprint => AMReX_mfab_eprint
     !procedure, private  :: get_array_func
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

  function cast_as_AMReX_factory(factory_polymorph) result(pf_AMReX_factory_obj)
    class(pf_factory_t), intent(in), target :: factory_polymorph
    type(pf_AMReX_mfab_factory_t), pointer :: pf_AMReX_factory_obj
    
    select type(factory_polymorph)
    type is (pf_AMReX_mfab_factory_t)
       pf_AMReX_factory_obj => factory_polymorph
    end select
  end function cast_as_AMReX_factory

  !>  Subroutine to allocate the mfab and set its pararmeters
  subroutine AMReX_mfab_build(this, shape_in, ba, dm, geom)
    use pf_mod_comm_mpi
    class(pf_encap_t),              intent(inout)   :: this
    integer,                        intent(in   )   :: shape_in(:)  ! new def: [nx,(ny,nz,)ncomp,nghost] 
    type(amrex_boxarray),           intent(in   )   :: ba
    type(amrex_distromap),          intent(in   )   :: dm
    type(amrex_geometry),  target,  intent(in   )   :: geom

    ! local 
    integer :: k      ! iterator
    integer :: ierr
    !
    select type (this)
    class is (pf_AMReX_mfab_t)
       ! dimension check
       this%ndim   = SIZE(shape_in)-2   ! -2 to account for ncomp and nghost
       if (this%ndim .ne. amrex_spacedim) then
          call pf_stop(__FILE__,__LINE__,'bad dimension in amrex encap, ndim=',this%ndim)
       end if
       
       !
       this%arr_shape               = 1.0_pfdp
       this%arr_shape(1:this%ndim)  = shape_in(1:this%ndim)
       this%arr_shape(4)            = shape_in(this%ndim+1)   ! kinda double here
       this%ncomp                   = shape_in(this%ndim+1)   ! kinda double here
       this%nghost                  = shape_in(this%ndim+2)
       this%ndof                    = product(shape_in(1:this%ndim+1))
       
       !  Make a shape the size of the grid with ghost cells
       this%pack_size = 1
       this%pack_size(1:this%ndim+1) = shape_in(1:this%ndim+1)
       do k=1,this%ndim
           this%pack_size(k)=this%pack_size(k)+2*this%nghost 
       end do

       ! Only point to amrex_geometry object here, needed in sweeper for fill_boundary
       this%geom = geom
       ! DistributionMapping "dm" & Boxarray "ba" build in main and used for multifab build are available via the mfab:
       !  dm -> this%mfab%dm & ba -> this%mfab%ba

       !> Build data multifabs
       call amrex_multifab_build(this%mfab, ba, dm, this%ncomp, this%nghost)

       !> Build the index for the flat array
       ! we use C++ indexing convention for am_to_flat since bx%lo,bx%hi use that convention 
       allocate(this%am_to_flat(this%arr_shape(1), this%arr_shape(2), this%arr_shape(3), this%arr_shape(4)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
       this%am_to_flat = reshape( (/ (k, k=1,product(this%arr_shape)) /), this%arr_shape)

      end select

  end subroutine AMReX_mfab_build

  !> Subroutine to  create a single array
  subroutine AMReX_mfab_create_single(this, x, level_index, lev_shape)
    class(pf_AMReX_mfab_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    ! local
    integer :: ierr
    
    !> allocate
    allocate(pf_AMReX_mfab_t::x,stat=ierr)
    
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    
    call AMReX_mfab_build(x, lev_shape, this%ba, this%dm, this%geom)
    
  end subroutine AMReX_mfab_create_single

  !> Subroutine to create an array of arrays
  subroutine AMReX_mfab_create_array(this, x, n, level_index,  lev_shape)
    class(pf_AMReX_mfab_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    ! local
    integer :: i,ierr
    
    !> allocate
    allocate(pf_AMReX_mfab_t::x(n),stat=ierr)
    
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    do i = 1, n
       call AMReX_mfab_build(x(i), lev_shape, this%ba, this%dm, this%geom)
    end do
  end subroutine AMReX_mfab_create_array

  !>  Subroutine to destroy array
  subroutine AMReX_mfab_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pf_AMReX_mfab_t), pointer :: pf_AMReX_mfab_obj
    
    pf_AMReX_mfab_obj => cast_as_AMReX_mfab(encap)
    
    call amrex_multifab_destroy(pf_AMReX_mfab_obj%mfab)
    call amrex_geometry_destroy(pf_AMReX_mfab_obj%geom) 
    nullify(pf_AMReX_mfab_obj)

  end subroutine AMReX_mfab_destroy

  !> Subroutine to destroy an single array
  subroutine AMReX_mfab_destroy_single(this, x)
    class(pf_AMReX_mfab_factory_t), intent(inout)              :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    !
    select type (x)
    class is (pf_AMReX_mfab_t)
       call amrex_multifab_destroy(x%mfab)
       call amrex_geometry_destroy(x%geom)
    end select
    deallocate(x)
  end subroutine AMReX_mfab_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine AMReX_mfab_destroy_array(this, x)
    class(pf_AMReX_mfab_factory_t), intent(inout)      :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer                                            :: i
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

    ! 1 defines icomp in AMReX - icomp defines which components are set (here 1:this%ncomp)
    call this%mfab%setval(val,1,this%ncomp,this%nghost)       
  end subroutine AMReX_mfab_setval

  !> Subroutine to copy an array
  subroutine AMReX_mfab_copy(this, src, flags)
    class(pf_AMReX_mfab_t),    intent(inout)   :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags

    integer ng, nc
    integer this_comp, src_comp

    
    select type(src)
    type is (pf_AMReX_mfab_t)
       ng = this%nghost              ! ghost cells
       nc = this%ncomp               ! number of components
       ! the 1,1 defines the component to copy from src
       this_comp = 1
       src_comp  = 1
       !call this%mfab%copy(src%mfab, src_comp, this_comp, nc, ng)
       call this%mfab%parallel_copy(src%mfab, src%geom)              
    class default
       call pf_stop(__FILE__,__LINE__,'Type error')
    end select
  end subroutine AMReX_mfab_copy

  !> Subroutine to pack mfab into a flat array for sending
  subroutine AMReX_mfab_pack(this, z, flags)
    class(pf_AMReX_mfab_t), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    ! local
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: mfab_data
    type(amrex_box) :: bx    
    type(amrex_mfiter) :: mfi
    integer :: lo(3), hi(3), n(3), i
    integer, allocatable :: flat_idx(:)

    !> loop over boxes
    if (this%mfab%ba%nboxes() == 1) then    
      ! single box -> single reshape is sufficient
      mfab_data=>this%mfab%dataPtr(0)       ! get c-pointer to data - 0 index due to c++ indexing
      bx = this%mfab%ba%get_box(0)          ! get box to get rid of ghost cells
      z = reshape(mfab_data(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),1:this%ncomp),[product(this%arr_shape)])
    else    ! multiple boxes -> use mfiter to loop over boxes
      call amrex_mfiter_build(mfi, this%mfab, tiling=.true.)
      z = 0.0_pfdp
      do while (mfi%next())
        ! get current box and data pointer
        bx = mfi%tilebox()
        mfab_data => this%mfab%dataptr(mfi)
        ! get index ranges
        lo = bx%lo          ! min and max index per dimension in c++ style
        hi = bx%hi          ! min and max index per dimension in c++ style
        lo(1:amrex_spacedim) = lo(1:amrex_spacedim) + 1     ! +1 to go  from c++ indexing to fortran indexing   
        hi(1:amrex_spacedim) = hi(1:amrex_spacedim) + 1     ! +1 to go  from c++ indexing to fortran indexing
        n = hi - lo + 1     ! dofs per dim inside current box 
        ! generate the index for the flat array
        if (allocated(flat_idx)) deallocate(flat_idx)
        allocate(flat_idx(product(n)))  
        flat_idx = reshape(this%am_to_flat(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:this%ncomp) , [product(n)*this%ncomp])
        ! assign the data to the flat array
        z(flat_idx) = reshape(mfab_data(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),1:this%ncomp),[product(n)*this%ncomp])
        deallocate(flat_idx)
      end do   
      call amrex_mfiter_destroy(mfi)
    end if
    
  end subroutine AMReX_mfab_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine AMReX_mfab_unpack(this, z, flags)
    class(pf_AMReX_mfab_t), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags

    real(pfdp),  pointer :: mfab_data(:,:,:,:)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi 
    integer :: lo(3), hi(3), n(3), i
    integer, allocatable :: flat_idx(:)

    if (this%mfab%ba%nboxes() == 1) then
      ! single box -> single reshape is sufficient
      mfab_data=>this%mfab%dataPtr(0)         ! get c-pointer to data - 0 index due to c++ indexing
      bx = this%mfab%ba%get_box(0)            ! get box to get rid of ghost cells
      mfab_data(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),1:this%ncomp) = reshape(z,this%arr_shape)
    else
      ! multiple boxes -> use mfiter to loop over boxes
      call amrex_mfiter_build(mfi, this%mfab, tiling=.true.)
      do while (mfi%next())
        ! get box and data pointer
        bx = mfi%tilebox()
        mfab_data => this%mfab%dataptr(mfi)
        ! get index ranges
        lo = bx%lo          ! min and max index per dimension in c++ style
        hi = bx%hi          ! min and max index per dimension in c++ style
        lo(1:amrex_spacedim) = lo(1:amrex_spacedim) + 1     ! +1 to go  from c++ indexing to fortran indexing   
        hi(1:amrex_spacedim) = hi(1:amrex_spacedim) + 1     ! +1 to go  from c++ indexing to fortran indexing
        n = hi - lo + 1     ! dofs per dim inside current box
        ! generate the index for the flat array
        if (allocated(flat_idx)) deallocate(flat_idx)
        allocate(flat_idx(product(n)))  
        flat_idx = reshape(this%am_to_flat(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:this%ncomp) , [product(n)*this%ncomp])
        ! assign the data from the flat array
        mfab_data(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),1:this%ncomp) = reshape(z(flat_idx),[n(1), n(2), n(3), this%ncomp])
        deallocate(flat_idx)
      end do   
      call amrex_mfiter_destroy(mfi)
    end if

  end subroutine AMReX_mfab_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function AMReX_mfab_norm(this, flags) result (norm)
    class(pf_AMReX_mfab_t), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    
    norm=1.0
    !if (present(flags)) then
      norm = this%mfab%norm0()
    !else
       !norm = this%mfab%norm2()
    !end if
  end function AMReX_mfab_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine AMReX_mfab_axpy(this, a, x, flags)
    class(pf_AMReX_mfab_t),    intent(inout)   :: this
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
    ! local
    real(pfdp), allocatable :: mfab_data_flat(:)
    
    allocate(mfab_data_flat(this%ndof))
    call this%pack(mfab_data_flat)
    !print *, mfab_data_flat         ! prints without ghost cells
    print *, mfab_data_flat(1:10)         ! prints without ghost cells
    deallocate(mfab_data_flat)  
  
  end subroutine AMReX_mfab_eprint

  !>  Helper functions to return the array part
  function get_array_func(x) result(r)
    class(pf_encap_t), target, intent(in) :: x
    ! local
    real(pfdp), pointer :: r(:,:,:,:)
    integer :: i, lo(3), hi(3)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    !
    select type (x)
    type is (pf_amrex_mfab_t)
      call amrex_mfiter_build(mfi, x%mfab, tiling=.true.)
      ! amrex MultiFabIterator(mfi) 
      do while (mfi%next())
        bx = mfi%tilebox()  
        lo = bx%lo
        hi = bx%hi
        do i = 1,3
          if (lo(i) .ne. hi(i)) then
            lo(i) = lo(i) + 1
            hi(i) = hi(i) + 1
          end if
        end do
        r(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:x%ncomp) => x%mfab%dataPtr(mfi)
      end do
      call amrex_mfiter_destroy(mfi)
    end select
  end function get_array_func

end module pf_mod_AMReX_mfab