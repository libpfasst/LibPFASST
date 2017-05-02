! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

module mod_transfer

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_transfer': Routines for inter-PFASST level data transfer
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

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

  ! ================================================================
  ! Data
  ! ================================================================

  integer,dimension(:),allocatable,private :: ref_rats

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: transfer_interpolate,transfer_restrict

contains

  subroutine create_transfer()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_transfer:create_transfer(...)': Create refinement ratios
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: alloc_stat

    ! ================================================================
    ! Work
    ! ================================================================

    allocate(ref_rats(num_dim_space),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `ref_rats(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_transfer

  subroutine destroy_transfer()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_transfer:destroy_transfer(...)': Destroy refinement ratios
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : mpi_exit_gracefully

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: dealloc_stat

    ! ================================================================
    ! Work
    ! ================================================================

    deallocate(ref_rats,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `ref_rats' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_transfer

  subroutine transfer_interpolate(c_encap_fn,c_encap_cr,level_fn,f_ctxn_c,level_cr,c_ctxr_c,time)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_transfer:transfer_interpolate(...)': Given coarse data, get fine data
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use petsc,only : DMDAVecGetArrayF90,DMDAVecRestoreArrayF90,INSERT_VALUES
    use pfasst,only : pfdp
    use omp_lib

    ! Project
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully
    use mod_ctx,only : num_step_space,num_var_dep,all_f_ctx_
    use mod_encap,only : type_encap
    use mod_interp,only : interp_wid,weights

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    type(c_ptr),intent(in),value :: c_ctxr_c,f_ctxn_c,c_encap_cr,c_encap_fn
    integer,intent(in) :: level_cr,level_fn
    real(pfdp),intent(in) :: time

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: i_var_dep,petsc_stat,step_space_x_cr,step_space_y_cr,step_space_x_cr1,step_space_y_cr1,step_space_x_fn,step_space_y_fn
    type(type_encap),pointer :: f_encap_cr,f_encap_fn
    type(DMDALocalInfo),dimension(DMDA_LOCAL_INFO_SIZE) :: info_cr,info_fn
    real(pfdp),dimension(:,:,:),pointer :: f_grid_arr_2d_cr,f_grid_arr_2d_fn

    ! ================================================================
    ! Work
    ! ================================================================
    
    ! Access PFASST encapsulations
    call c_f_pointer(c_encap_cr,f_encap_cr)
    call c_f_pointer(c_encap_fn,f_encap_fn)

    ! Define refinement ratios
    ref_rats = num_step_space(:,f_encap_fn%level)/num_step_space(:,f_encap_cr%level)

    ! Sanity check for `x'-refinement ratio
    if((ref_rats(1) .ne. 1) .and. (ref_rats(1) .ne. 2)) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Fine to coarse `x'-refinement ratio must be less than 3.")
    end if ! `(ref_rats(1) .ne. 1) .and. (ref_rats(1) .ne. 2)'

    ! Sanity check for `y'-refinement ratio
    if((ref_rats(2) .ne. 1) .and. (ref_rats(2) .ne. 2)) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Fine to coarse `y'-refinement ratio must be less than 3.")
    end if ! `(ref_rats(2) .ne. 1) .and. (ref_rats(2) .ne. 2)'

    call DMDAGetLocalInfoF90(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,info_cr,petsc_stat)
    call DMDAGetLocalInfoF90(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,info_fn,petsc_stat)

    ! Assign ghost values to local (w.r.t. MPI) coarse grid data
    call DMGlobalToLocalBegin(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,f_encap_cr%u,INSERT_VALUES,all_f_ctx_(f_encap_cr%level)%f_ctx%pad,petsc_stat)
    call DMGlobalToLocalEnd(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,f_encap_cr%u,INSERT_VALUES,all_f_ctx_(f_encap_cr%level)%f_ctx%pad,petsc_stat)

    ! Dimensionality dependent interpolation routine
    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      ! Access local (w.r.t. MPI) *ghosted* PETSc vector as multi-dimensional array
      call DMDAVecGetArrayF90(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,all_f_ctx_(f_encap_cr%level)%f_ctx%pad,f_grid_arr_2d_cr,petsc_stat)

      ! Access global (w.r.t. MPI) ghost-free PETSc vector as multi-dimensional array
      call DMDAVecGetArrayF90(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,f_encap_fn%u,f_grid_arr_2d_fn,petsc_stat)

      ! Point injection in `(x,y)'
!      !$OMP parallel do
      do step_space_y_fn = info_fn(DMDA_LOCAL_INFO_YS),info_fn(DMDA_LOCAL_INFO_YS)+info_fn(DMDA_LOCAL_INFO_YM)-ref_rats(2),ref_rats(2)
        step_space_y_cr = step_space_y_fn/ref_rats(2)
        ! Get value at `(x,y)'
        do step_space_x_fn = info_fn(DMDA_LOCAL_INFO_XS),info_fn(DMDA_LOCAL_INFO_XS)+info_fn(DMDA_LOCAL_INFO_XM)-ref_rats(1),ref_rats(1)
          step_space_x_cr = step_space_x_fn/ref_rats(1)
          do i_var_dep = 0,num_var_dep-1
            f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn) = f_grid_arr_2d_cr(i_var_dep,step_space_x_cr,step_space_y_cr)
          end do ! `i_var_dep'
        end do ! `step_space_x_fn'
      end do ! `step_space_y_fn'
!      !$OMP end parallel do

      ! Straight lines in `x'
      if(ref_rats(1) .eq. 2) then
!        !$OMP parallel do
        do step_space_y_fn = info_fn(DMDA_LOCAL_INFO_YS),info_fn(DMDA_LOCAL_INFO_YS)+info_fn(DMDA_LOCAL_INFO_YM)-ref_rats(2),ref_rats(2)
          step_space_y_cr = step_space_y_fn/ref_rats(2)
          do step_space_x_fn = info_fn(DMDA_LOCAL_INFO_XS)+1,info_fn(DMDA_LOCAL_INFO_XS)+info_fn(DMDA_LOCAL_INFO_XM)-1,ref_rats(1)
            step_space_x_cr = (step_space_x_fn-1)/ref_rats(1)
            do i_var_dep = 0,num_var_dep-1
              f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn) = 0d0
              ! Get new value at intermediate `x'
              do step_space_x_cr1 = step_space_x_cr-interp_wid+1,step_space_x_cr+interp_wid
                f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn) = f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn)+weights(step_space_x_cr1-step_space_x_cr+interp_wid,1)*f_grid_arr_2d_cr(i_var_dep,step_space_x_cr1,step_space_y_cr)
              end do ! `step_space_x_cr1'
            end do ! `i_var_dep'
          end do ! `step_space_x_fn'
        end do ! `step_space_y_fn'
!        !$OMP end parallel do
      end if ! `ref_rats(1) .eq. 2'

      ! Straight lines in `y'
      if(ref_rats(2) .eq. 2) then
!        !$OMP parallel do
        do step_space_y_fn = info_fn(DMDA_LOCAL_INFO_YS)+1,info_fn(DMDA_LOCAL_INFO_YS)+info_fn(DMDA_LOCAL_INFO_YM)-1,ref_rats(2)
          step_space_y_cr = (step_space_y_fn-1)/ref_rats(2)
          do step_space_x_fn = info_fn(DMDA_LOCAL_INFO_XS),info_fn(DMDA_LOCAL_INFO_XS)+info_fn(DMDA_LOCAL_INFO_XM)-ref_rats(1),ref_rats(1)
            step_space_x_cr = step_space_x_fn/ref_rats(1)
            do i_var_dep = 0,num_var_dep-1
              f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn) = 0d0
              ! Get new value at intermediate `y'
              do step_space_y_cr1 = step_space_y_cr-interp_wid+1,step_space_y_cr+interp_wid
                f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn) = f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn)+weights(step_space_y_cr1-step_space_y_cr+interp_wid,2)*f_grid_arr_2d_cr(i_var_dep,step_space_x_cr,step_space_y_cr1)
              end do ! `step_space_y_cr1'
            end do ! `i_var_dep'
          end do ! `step_space_x_fn'
        end do ! `step_space_y_fn'
!        !$OMP end parallel do
      end if ! `ref_rats(2) .eq. 2'

      ! Planes in `(x,y)'
      if((ref_rats(1) .eq. 2) .and. (ref_rats(2) .eq. 2)) then
!        !$OMP parallel do
        do step_space_y_fn = info_fn(DMDA_LOCAL_INFO_YS)+1,info_fn(DMDA_LOCAL_INFO_YS)+info_fn(DMDA_LOCAL_INFO_YM)-1,ref_rats(2)
          step_space_y_cr = (step_space_y_fn-1)/ref_rats(2)
          do step_space_x_fn = info_fn(DMDA_LOCAL_INFO_XS)+1,info_fn(DMDA_LOCAL_INFO_XS)+info_fn(DMDA_LOCAL_INFO_XM)-1,ref_rats(1)
            step_space_x_cr = (step_space_x_fn-1)/ref_rats(1)
            do i_var_dep = 0,num_var_dep-1
              f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn) = 0d0
              ! Get new value at intermediate `(x,y)'
              do step_space_y_cr1 = step_space_y_cr-interp_wid+1,step_space_y_cr+interp_wid
                do step_space_x_cr1 = step_space_x_cr-interp_wid+1,step_space_x_cr+interp_wid
                  f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn) = f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn)+weights(step_space_x_cr1-step_space_x_cr+interp_wid,1)*weights(step_space_y_cr1-step_space_y_cr+interp_wid,2)*f_grid_arr_2d_cr(i_var_dep,step_space_x_cr1,step_space_y_cr1)
                end do ! `step_space_x_cr1'
              end do ! `step_space_y_cr1'
            end do ! `i_var_dep'
          end do ! `step_space_x_fn'
        end do ! `step_space_y_fn'
!        !$OMP end parallel do
      end if ! `(ref_rats(1) .eq. 2) .and. (ref_rats(2) .eq. 2)'

      ! Return pointer to global (w.r.t. MPI) ghost-free PETSc vector as multi-dimensional array
      call DMDAVecRestoreArrayF90(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,f_encap_fn%u,f_grid_arr_2d_fn,petsc_stat)

      ! Return pointer to local (w.r.t. MPI) *ghosted* PETSc vector as multi-dimensional array
      call DMDAVecRestoreArrayF90(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,all_f_ctx_(f_encap_cr%level)%f_ctx%pad,f_grid_arr_2d_cr,petsc_stat)
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine transfer_interpolate

  subroutine transfer_restrict(c_encap_fn,c_encap_cr,level_fn,f_ctxn_c,level_cr,c_ctxr_c,time)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_transfer:transfer_restrict(...)': Given fine data, get coarse data
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use petsc,only : DMDAVecGetArrayF90,DMDAVecRestoreArrayF90,INSERT_VALUES
    use pfasst,only : pfdp

    ! Project
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully
    use mod_ctx,only : num_step_space,num_var_dep,all_f_ctx_
    use mod_encap,only : type_encap
    use mod_interp,only : interp_wid

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    type(c_ptr),intent(in),value :: f_ctxn_c,c_ctxr_c,c_encap_fn,c_encap_cr
    integer,intent(in) :: level_fn,level_cr
    real(pfdp),intent(in) :: time

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: i_var_dep,petsc_stat,step_space_x_fn,step_space_y_fn,step_space_x_cr,step_space_y_cr
    type(type_encap),pointer :: f_encap_fn,f_encap_cr
    type(DMDALocalInfo),dimension(DMDA_LOCAL_INFO_SIZE) :: info_fn,info_cr
    real(pfdp),dimension(:,:,:),pointer :: f_grid_arr_2d_cr,f_grid_arr_2d_fn

    ! ================================================================
    ! Work
    ! ================================================================
    
    ! Access PFASST encapsulations
    call c_f_pointer(c_encap_fn,f_encap_fn)
    call c_f_pointer(c_encap_cr,f_encap_cr)

    ! Define refinement ratios
    ref_rats = num_step_space(:,f_encap_fn%level)/num_step_space(:,f_encap_cr%level)

    ! Sanity check for `x'-refinement ratio
    if((ref_rats(1) .ne. 1) .and. (ref_rats(1) .ne. 2)) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Fine to coarse `x'-refinement ratio must be less than 3.")
    end if ! `(ref_rats(1) .ne. 1) .and. (ref_rats(1) .ne. 2)'

    ! Sanity check for `y'-refinement ratio
    if((ref_rats(2) .ne. 1) .and. (ref_rats(2) .ne. 2)) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Fine to coarse `y'-refinement ratio must be less than 3.")
    end if ! `(ref_rats(2) .ne. 1) .and. (ref_rats(2) .ne. 2)'

    call DMDAGetLocalInfoF90(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,info_fn,petsc_stat)
    call DMDAGetLocalInfoF90(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,info_cr,petsc_stat)

    ! Assign ghost values to local (w.r.t. MPI) coarse grid data
    call DMGlobalToLocalBegin(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,f_encap_fn%u,INSERT_VALUES,all_f_ctx_(f_encap_fn%level)%f_ctx%pad,petsc_stat)
    call DMGlobalToLocalEnd(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,f_encap_fn%u,INSERT_VALUES,all_f_ctx_(f_encap_fn%level)%f_ctx%pad,petsc_stat)

    ! Dimensionality dependent transfer_restriction routine
    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      ! Access local (w.r.t. MPI) *ghosted* PETSc vector as multi-dimensional array
      call DMDAVecGetArrayF90(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,all_f_ctx_(f_encap_fn%level)%f_ctx%pad,f_grid_arr_2d_fn,petsc_stat)

      ! Access global (w.r.t. MPI) ghost-free PETSc vector as multi-dimensional array
      call DMDAVecGetArrayF90(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,f_encap_cr%u,f_grid_arr_2d_cr,petsc_stat)

      ! Point injection in `(x,y)'
!      !$OMP parallel do
      do step_space_y_cr = info_cr(DMDA_LOCAL_INFO_YS),info_cr(DMDA_LOCAL_INFO_YS)+info_cr(DMDA_LOCAL_INFO_YM)-1
        step_space_y_fn = step_space_y_cr*ref_rats(2)
        ! Get value at `(x,y)'
        do step_space_x_cr = info_cr(DMDA_LOCAL_INFO_XS),info_cr(DMDA_LOCAL_INFO_XS)+info_cr(DMDA_LOCAL_INFO_XM)-1
          step_space_x_fn = step_space_x_cr*ref_rats(1)
          do i_var_dep = 0,num_var_dep-1
            f_grid_arr_2d_cr(i_var_dep,step_space_x_cr,step_space_y_cr) = f_grid_arr_2d_fn(i_var_dep,step_space_x_fn,step_space_y_fn)
          end do ! `i_var_dep'
        end do ! `step_space_x_cr'
      end do ! `step_space_y_cr'
!      !$OMP end parallel do

      ! Return pointer to global (w.r.t. MPI) ghost-free PETSc vector as multi-dimensional array
      call DMDAVecRestoreArrayF90(all_f_ctx_(f_encap_cr%level)%f_ctx%mesh,f_encap_cr%u,f_grid_arr_2d_cr,petsc_stat)

      ! Return pointer to local (w.r.t. MPI) *ghosted* PETSc vector as multi-dimensional array
      call DMDAVecRestoreArrayF90(all_f_ctx_(f_encap_fn%level)%f_ctx%mesh,all_f_ctx_(f_encap_fn%level)%f_ctx%pad,f_grid_arr_2d_fn,petsc_stat)
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine transfer_restrict

end module mod_transfer
