! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

module mod_ctx

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_ctx': Define degrees of freedom (independent and dependent)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! External
  use petsc,only : DM,Vec
  use pfasst,only : pfdp

  ! ================================================================
  ! Safety
  ! ================================================================

  implicit none

  ! ================================================================ 
  ! Parameters
  ! ================================================================ 

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! Number of spatial steps
  integer,dimension(:,:),allocatable,save,public :: num_step_space
  ! Spatial step sizes (nonpositives values are ignored)
  real(pfdp),dimension(:,:),allocatable,save,public :: size_step_space
  ! Space domain sizes (positive `size_step_space' lead to these arguments being inferred)
  real(pfdp),dimension(:),allocatable,save,public :: size_dom_space
  ! Space domain origins
  real(pfdp),dimension(:),allocatable,save,public :: orig_dom_space

  ! Number of sol components
  integer,save,public :: num_var_dep
  
  ! List of arguments
  namelist /arg_core_ctx/ num_step_space,size_step_space,size_dom_space,orig_dom_space,num_var_dep

  ! ================================================================
  ! Data
  ! ================================================================

  ! Solution norm
  real(pfdp),dimension(:),allocatable,public :: norm_

  ! Stencil width from derivative and interpolation order
  integer,save,public :: stenc_wid_

  ! PFASST level context
  type,public :: type_ctx
    integer :: level
    type(DM) :: mesh
    type(Vec) :: pad,ref
    real(pfdp) :: size_step_time
    real(pfdp),dimension(:),pointer :: hdf
    integer,dimension(:),pointer :: shape,count,offset
  end type type_ctx

  ! All Fortran pointed-to PFASST level contexts
  type,private :: f_type_ctx
    type(type_ctx),pointer :: f_ctx
  end type f_type_ctx
  type(f_type_ctx),dimension(:),pointer,public :: all_f_ctx_

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_ctx,destroy_ctx,ctx_write
  private :: ctx_read_arg_core_ctx

contains

  subroutine create_ctx(io_unit,name_file_input_pfasst,pfasst_pf)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_ctx:create_ctx(...)': Prepare memory etc. for spatial DoFs and dependent sols
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_loc

    ! External
    use mpi,only : MPI_INTEGER
    use petsc,only : DM_BOUNDARY_PERIODIC,DMDA_STENCIL_BOX,PETSC_NULL_INTEGER
    use pfasst,only : pf_pfasst_t,pfdp
    use omp_lib

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : max_len_char_,num_dim_space,num_rank_space,mpi_space_,mpi_exit_gracefully
    use mod_deriv,only : deriv_wid,deriv_set_coeffs
    use mod_interp,only : interp_wid,interp_set_weights

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: io_unit
    character(len = *),intent(in) :: name_file_input_pfasst
    type(pf_pfasst_t),intent(in) :: pfasst_pf

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: num_arg_read,io_stat,mpi_stat,alloc_stat,i_level,petsc_stat,i_dim_space
    character(len = max_len_char_) :: cmd_line_arg,namelist_arg,io_msg

    real(pfdp),dimension(:),pointer :: f_flat_arr

    ! ================================================================
    ! Work
    ! ================================================================

    ! Get input arguments for DoFs definitions
    call ctx_read_arg_core_ctx(io_unit,name_file_input_pfasst)

    ! Sanity check for number of steps
    do i_level = 1,num_level
      do i_dim_space = 1,num_dim_space
        if(num_step_space(i_dim_space,i_level) .lt. num_rank_space(i_dim_space)) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"There are not enough spatial steps for the number of PEs assigned.")
        end if ! `num_step_space(i_dim_space,i_level) .lt. num_rank_space(i_dim_space)'
      end do ! `i_dim_space'
    end do ! `i_level'

    ! Space step size (a nonpositive value is overwritten)
    do i_level = 1,num_level
      do i_dim_space = 1,num_dim_space
        if(size_step_space(i_dim_space,i_level) .le. 0) then
          size_step_space(i_dim_space,i_level) = size_dom_space(i_dim_space)/num_step_space(i_dim_space,i_level)
        end if ! `size_step_space(i_dim_space,i_level) .le. 0'
      end do ! `i_dim_space'
    end do ! `i_level'

    ! Space domain size (if possible then derived from spatial step sizes)
    do i_level = 1,num_level
      do i_dim_space = 1,num_dim_space
        if(size_step_space(i_dim_space,i_level) .gt. 0) then
          size_dom_space(i_dim_space) = num_step_space(i_dim_space,i_level)*size_step_space(i_dim_space,i_level)
        end if ! `size_step_space(i_dim_space,i_level) .gt. 0'
      end do ! `i_dim_space'
    end do ! `i_level'

    ! Reserve memory for vector norm
    allocate(norm_(num_var_dep),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `norm_(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Initialize derivative coefficients
    call deriv_set_coeffs(size_step_space)

    ! Assign values to interpolation weights
    call interp_set_weights()

    ! Stencil width from derivative and interpolation order
    stenc_wid_ = max(deriv_wid,interp_wid)

    ! One `c_loc'-able Fortran pointer to each PFASST level context
    allocate(all_f_ctx_(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_f_ctx_(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! One `c_loc'-able Fortran pointer to each PFASST level context
    do i_level = 1,num_level
      allocate(all_f_ctx_(i_level)%f_ctx,stat = alloc_stat)
      if(alloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_f_ctx_(...)%f_ctx' failed.")
      end if ! `alloc_stat .ne. 0'
    end do ! `i_level'

    ! Assign associated PFASST level
    do i_level = 1,num_level
      all_f_ctx_(i_level)%f_ctx%level = i_level
    end do ! `i_level'

    ! Create spatial mesh for space of PFASST level
    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      do i_level = 1,num_level
        call DMDACreate2d(mpi_space_%comm,DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,DMDA_STENCIL_BOX,num_step_space(1,i_level),num_step_space(2,i_level),num_rank_space(1),num_rank_space(2),num_var_dep,stenc_wid_,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,all_f_ctx_(i_level)%f_ctx%mesh,petsc_stat)
      end do ! `i_level'
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! Allocate memory for local vector padded by ghosts
    do i_level = 1,num_level
      call DMCreateLocalVector(all_f_ctx_(i_level)%f_ctx%mesh,all_f_ctx_(i_level)%f_ctx%pad,petsc_stat)
    end do ! `i_level'

    ! Allocate memory for reference solution (no ghosts)
    do i_level = 1,num_level
      call DMCreateGlobalVector(all_f_ctx_(i_level)%f_ctx%mesh,all_f_ctx_(i_level)%f_ctx%ref,petsc_stat)
    end do ! `i_level'

    ! Global (w.r.t. MPI) file dimensions
    do i_level = 1,num_level
      ! Reserve memory
      allocate(all_f_ctx_(i_level)%f_ctx%shape(1),stat = alloc_stat)
      if(alloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_f_ctx_(...)%f_ctx%shape(...)'' failed.")
      end if ! `alloc_stat .ne. 0'
    end do ! `i_level'

    ! Assign values
    do i_level = 1,num_level
      call VecGetSize(all_f_ctx_(i_level)%f_ctx%ref,all_f_ctx_(i_level)%f_ctx%shape(1),petsc_stat)
    end do ! `i_level'

    ! Local (w.r.t. MPI) memory dimensions
    do i_level = 1,num_level
      ! Reserve memory
      allocate(all_f_ctx_(i_level)%f_ctx%count(1),stat = alloc_stat)
      if(alloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_f_ctx_(...)%f_ctx%count(...)'' failed.")
      end if ! `alloc_stat .ne. 0'
    end do ! `i_level'

    ! Assign values
    do i_level = 1,num_level
      call VecGetLocalSize(all_f_ctx_(i_level)%f_ctx%ref,all_f_ctx_(i_level)%f_ctx%count(1),petsc_stat)
    end do ! `i_level'

    ! Local (w.r.t. MPI) memory offsets
    do i_level = 1,num_level
      ! Reserve memory for local (w.r.t. MPI) memory offsets
      allocate(all_f_ctx_(i_level)%f_ctx%offset(1),stat = alloc_stat)
      if(alloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_f_ctx_(...)%f_ctx%offset(...)'' failed.")
      end if ! `alloc_stat .ne. 0'
    end do ! `i_level'

    ! Assign data
    do i_level = 1,num_level
      call VecGetOwnershipRange(all_f_ctx_(i_level)%f_ctx%ref,all_f_ctx_(i_level)%f_ctx%offset(1),PETSC_NULL_INTEGER,petsc_stat)
    end do ! `i_level'

    do i_level = 1,num_level
      allocate(all_f_ctx_(i_level)%f_ctx%hdf(num_var_dep*all_f_ctx_(i_level)%f_ctx%count(1)),stat = alloc_stat)
      if(alloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_f_ctx_(...)%f_ctx%hdf(...)' failed.")
      end if ! `alloc_stat .ne. 0'
    end do ! `i_level'

    ! Assign location of Fortran pointer to C-pointer
    do i_level = 1,num_level
      pfasst_pf%levels(i_level)%ctx = c_loc(all_f_ctx_(i_level)%f_ctx) ! Define C-pointer
    end do ! `i_level'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_ctx

  subroutine destroy_ctx(pfasst_pf)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_ctx:destroy_ctx(...)': Call PETSc routine to free memory reserved by spatial meshes
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================ 
    ! Modules
    ! ================================================================ 

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use pfasst,only : pf_pfasst_t

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully

    ! ================================================================ 
    ! Safety
    ! ================================================================ 

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    type(pf_pfasst_t),target,intent(in) :: pfasst_pf

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: dealloc_stat,i_level,petsc_stat,i_dim_space

    ! ================================================================
    ! Work
    ! ================================================================

    ! Get location of C-pointer
    do i_level = 1,num_level
      ! Access PFASST level context
      call c_f_pointer(pfasst_pf%levels(i_level)%ctx,all_f_ctx_(i_level)%f_ctx)
    end do ! `i_level'

    do i_level = 1,num_level
      deallocate(all_f_ctx_(i_level)%f_ctx%hdf,stat = dealloc_stat)
      if(dealloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_f_ctx_(...)%f_ctx%hdf' failed.")
      end if ! `dealloc_stat .ne. 0'
    end do ! `i_level'

    ! Local (w.r.t. MPI) memory offsets
    do i_level = 1,num_level
      deallocate(all_f_ctx_(i_level)%f_ctx%offset,stat = dealloc_stat)
      if(dealloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_f_ctx_(...)%f_ctx%offset'' failed.")
      end if ! `dealloc_stat .ne. 0'
    end do ! `i_level'

    ! Local (w.r.t. MPI) memory dimensions
    do i_level = 1,num_level
      deallocate(all_f_ctx_(i_level)%f_ctx%count,stat = dealloc_stat)
      if(dealloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_f_ctx_(...)%f_ctx%count'' failed.")
      end if ! `dealloc_stat .ne. 0'
    end do ! `i_level'

    ! Local (w.r.t. MPI) file dimensions
    do i_level = 1,num_level
      deallocate(all_f_ctx_(i_level)%f_ctx%shape,stat = dealloc_stat)
      if(dealloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_f_ctx_(...)%f_ctx%shape'' failed.")
      end if ! `dealloc_stat .ne. 0'
    end do ! `i_level'

    ! Remove reference solution
    do i_level = 1,num_level
      call VecDestroy(all_f_ctx_(i_level)%f_ctx%ref,petsc_stat)
    end do ! `i_level'

    ! Remove reference solution
    do i_level = 1,num_level
      call VecDestroy(all_f_ctx_(i_level)%f_ctx%pad,petsc_stat)
    end do ! `i_level'

    ! Deallocate spatial meshes
    do i_level = 1,num_level
      call DMDestroy(all_f_ctx_(i_level)%f_ctx%mesh,petsc_stat)
      if(petsc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_f_ctx_(...)%mesh' failed.")
      end if ! `petsc_stat .ne. 0'
    end do ! `i_level'

    ! One `c_loc'-able Fortran pointer to each PFASST level context
    do i_level = 1,num_level
      deallocate(all_f_ctx_(i_level)%f_ctx,stat = dealloc_stat)
      if(dealloc_stat .ne. 0) then
        call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_ctxs_(...)%f_ctx' failed.")
      end if ! `dealloc_stat .ne. 0'
    end do ! `i_level'

    ! One `c_loc'-able Fortran pointer to each PFASST level context
    deallocate(all_f_ctx_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_f_ctx_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Destruction of memory for vector norm
    deallocate(norm_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `norm_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Space domain origins
    deallocate(orig_dom_space,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `orig_dom_space' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Space domain sizes (positive `size_step_space' lead to these arguments being inferred)
    deallocate(size_dom_space,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `size_dom_space' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Spatial step sizes (nonpositives values are ignored)
    deallocate(size_step_space,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `size_step_space' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Number of spatial steps
    deallocate(num_step_space,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `num_step_space' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_ctx

  subroutine ctx_read_arg_core_ctx(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_ctx:ctx_read_arg_core_ctx(...)': Get DoF arguments from file or command line input
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_INTEGER,MPI_DOUBLE_PRECISION

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : max_len_char_,num_dim_space,num_rank_space,mpi_world_,mpi_exit_gracefully

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: io_unit
    character(len = *),intent(in) :: name_file_input_pfasst

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: num_arg_read,io_stat,mpi_stat,alloc_stat,i_level,petsc_stat,i_dim_space
    character(len = max_len_char_) :: cmd_line_arg,namelist_arg,io_msg

    ! ================================================================
    ! Work
    ! ================================================================

    ! Number of spatial steps
    allocate(num_step_space(num_dim_space,num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `num_step_space(...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Default value
    num_step_space = 32

    ! Spatial step sizes (nonpositives values are ignored)
    allocate(size_step_space(num_dim_space,num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `size_step_space(...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Default value
    size_step_space = -1d0

    ! Space domain sizes (positive `size_step_space' lead to these arguments being inferred)
    allocate(size_dom_space(num_dim_space),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `size_dom_space(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Default value
    size_dom_space = 1d0

    ! Space domain origins
    allocate(orig_dom_space(num_dim_space),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `orig_dom_space(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Default value
    orig_dom_space = 0d0

    ! Default value for number of dependent sol components
    num_var_dep = 1

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_ctx)
      close(unit = io_unit)

      ! Read command line arguments to possibly overwrite content from input file
      num_arg_read = 0
      do
        call get_command_argument(num_arg_read,cmd_line_arg)
        if(len_trim(adjustL(cmd_line_arg)) .eq. 0) then
          exit
        end if ! `len_trim(adjustL(cmd_line_arg)) .eq. 0'
        if(num_arg_read .gt. 0) then
          ! Read namelist
          namelist_arg = "&arg_core_ctx "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_core_ctx,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(num_step_space,size(num_step_space),MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(size_step_space,size(size_step_space),MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(size_dom_space,size(size_dom_space),MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(orig_dom_space,size(orig_dom_space),MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(num_var_dep,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine ctx_read_arg_core_ctx

  subroutine ctx_write(level,name_file)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:ctx_write(...)': Print PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr

    ! External
    use hdf5,only : hSize_t,hid_t,h5open_f,h5pCreate_f,H5P_FILE_ACCESS_F,h5pSet_fapl_mpio_f,h5fCreate_f,H5F_ACC_TRUNC_F,h5pClose_f,h5sCreate_simple_f,h5dCreate_f,H5T_NATIVE_DOUBLE,h5dGet_space_f,h5sSelect_hyperslab_f,H5S_SELECT_SET_F,H5P_DATASET_XFER_F,h5pSet_dxpl_mpio_f,H5FD_MPIO_COLLECTIVE_F,h5dWrite_f,h5sClose_f,h5dClose_f,h5pClose_f,h5fClose_f,h5close_f
    use pfasst,only : pfdp

    ! Project
    use mod_mpi,only : mpi_space_

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: level
    character(len = *),intent(in) :: name_file

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: hdf_stat
    integer(hid_t) :: props_id,file_id,file_space,data_id,mem_space

    ! ================================================================
    ! Work
    ! ================================================================

    ! Initialize HDF library (see [h5open_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5.html#Library-Open ))
    call h5open_f(hdf_stat)

    ! Create property list for file access (see [h5pCreate_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-Create ))
    call h5pCreate_f(H5P_FILE_ACCESS_F,props_id,hdf_stat)
    ! Store MPI output comm. info. in property list (see [h5pSet_fapl_mpio_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetFaplMpio ))
    call h5pSet_fapl_mpio_f(props_id,mpi_space_%comm,mpi_space_%info,hdf_stat)
    ! Create HDF file (see [h5Create_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5F.html#File-Create ))
    call h5fCreate_f(trim(adjustL(name_file)),H5F_ACC_TRUNC_F,file_id,hdf_stat,access_prp = props_id)
    ! Terminate access to property list (see [h5pClose_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-Close ))
    call h5pClose_f(props_id,hdf_stat)

    ! Create simple file space and open it for access (see [h5sCreate_simple_f]( https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-CreateSimple ))
    call h5sCreate_simple_f(1,int(all_f_ctx_(level)%f_ctx%shape,hSize_t),file_space,hdf_stat)
    ! Create data set and link it to a location in the file (see [h5dCreate_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Create ))
    call h5dCreate_f(file_id,"hdf",H5T_NATIVE_DOUBLE,file_space,data_id,hdf_stat)
    ! Release and terminate access to file space (see [h5sClose_f]( https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-Close ))
    call h5sClose_f(file_space,hdf_stat)

    ! Make copy of file space that is associated with data set belonging to `data_id' and then return `file_space' as identifier for this new data space (see [h5dGet_space_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-GetSpace ))
    call h5dGet_space_f(data_id,file_space,hdf_stat)
    ! Select hyperslab region from data space defined by `file_space' and add it to currently selected region (see [h5sSelect_hyperslab_f]( https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab ))
    call h5sSelect_hyperslab_f(file_space,H5S_SELECT_SET_F,int(all_f_ctx_(level)%f_ctx%offset,hSize_t),int(all_f_ctx_(level)%f_ctx%count,hSize_t),hdf_stat)

    ! Create simple memory space and open it for access (see [h5sCreate_simple_f]( https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-CreateSimple ))
    call h5sCreate_simple_f(1,int(all_f_ctx_(level)%f_ctx%count,hSize_t),mem_space,hdf_stat)
    ! Create property list for data-to-file transfer (see [h5pCreate_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-Create ))
    call h5pCreate_f(H5P_DATASET_XFER_F,props_id,hdf_stat)
    ! Set mode for data-to-file transfer to 'collective' (see [h5pSet_dxpl_mpio_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetDxplMpio ))
    call h5pSet_dxpl_mpio_f(props_id,H5FD_MPIO_COLLECTIVE_F,hdf_stat)
    ! Write data set collectively (see [h5dWrite_f]( https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Write ))
    call h5dWrite_f(data_id,H5T_NATIVE_DOUBLE,all_f_ctx_(level)%f_ctx%hdf,int(all_f_ctx_(level)%f_ctx%shape,hSize_t),hdf_stat,mem_space_id = mem_space,file_space_id = file_space,xfer_prp = props_id)

    ! Close memory space for flat data set
    call h5sClose_f(mem_space,hdf_stat)
    ! Close file space for flat data set
    call h5sClose_f(file_space,hdf_stat)
    ! Close flat data set
    call h5dClose_f(data_id,hdf_stat)
    ! Close property list for data-to-file transfer
    call h5pClose_f(props_id,hdf_stat)
    ! Close file that has properties
    call h5fClose_f(file_id,hdf_stat)

    ! Close Fortran inferface
    call h5close_f(hdf_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine ctx_write

end module mod_ctx
