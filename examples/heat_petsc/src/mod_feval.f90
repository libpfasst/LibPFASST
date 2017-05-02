! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

module mod_feval

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_feval': Define time steppers
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! External
  use petsc,only : Mat
  use pfasst,only : pfdp

  ! ================================================================
  ! Safety
  ! ================================================================

  implicit none

  ! ================================================================ 
  ! Parameters
  ! ================================================================ 

  ! Sweeper `1' for explicit, `2' for implicit, and `3' for IMEX
  integer,dimension(3),parameter,private :: avail_sweepers_ = (/1,2,3/)

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! `1' for explicit, `2' for implicit, and `3' for IMEX
  integer,save,public :: sweeper
  ! LU factorization in implicit sweeper
  logical,save,public :: use_LU

  ! List of arguments
  namelist /arg_core_feval/ sweeper,use_LU

  ! ================================================================
  ! Data
  ! ================================================================

  ! Number of stencil rows
  integer,save,public :: num_stenc_rows_
  ! Number of stencil columns
  integer,save,public :: num_stenc_cols_

  ! Linear operator type
  type,public :: lops_type
    type(MatStencil),dimension(:,:),allocatable :: rows
    type(MatStencil),dimension(:,:),allocatable :: cols
    real(pfdp),dimension(:,:),allocatable :: vals
  end type lops_type

  ! Matrix definitions
  type(lops_type),dimension(:,:),pointer,public :: mats_

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_feval,destroy_feval,feval_f1eval,feval_f2eval,feval_f2comp
  private :: feval_set_2d_res,feval_set_2d_jac

contains

  subroutine create_feval(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_feval:create_feval(...)': Create solver infrastructure
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : num_dim_space,mpi_sanity_check,mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,all_f_ctx_
    use mod_deriv,only : deriv_acc

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

    integer :: i_level,i_var_dep,alloc_stat,petsc_stat

		! ================================================================
		! Work
		! ================================================================

    call feval_read_arg_core_feval(io_unit,name_file_input_pfasst)

    ! Sanity check for sweeper to use
    call mpi_sanity_check(avail_sweepers_,sweeper)

    ! Number of stencil rows
    num_stenc_rows_ = 1

    ! Number of stencil columns
    num_stenc_cols_ = num_dim_space*deriv_acc+1

    ! Reserve space for matrices for linear equations
    allocate(mats_(num_var_dep,num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `mats_(...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Create one Jacobian for each PFASST level
    do i_level = 1,num_level
      do i_var_dep = 1,num_var_dep
        ! Create rows (fourth row for dependent variable index)
        allocate(mats_(i_var_dep,i_level)%rows(4,num_stenc_rows_),stat = alloc_stat)
        if(alloc_stat .ne. 0) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `mats_(...,...)%rows(...,...)' failed.")
        end if ! `alloc_stat .ne. 0'

        ! Create columns (fourth row for dependent variable index)
        allocate(mats_(i_var_dep,i_level)%cols(4,num_stenc_cols_),stat = alloc_stat)
        if(alloc_stat .ne. 0) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `mats_(...,...)%cols(...,...)' failed.")
        end if ! `alloc_stat .ne. 0'

        ! Create values (fourth row for dependent variable index)
        allocate(mats_(i_var_dep,i_level)%vals(num_stenc_rows_,num_stenc_cols_),stat = alloc_stat)
        if(alloc_stat .ne. 0) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `mats_(...,...)%vals(...,...)' failed.")
        end if ! `alloc_stat .ne. 0'
      end do ! `i_var_dep'
    end do ! `i_level'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_feval

  subroutine destroy_feval()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_feval:destroy_feval(...)': Destroy solver infrastructure
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : mpi_exit_gracefully
    use mod_pf,only : num_level
    use mod_ctx,only : num_var_dep,type_ctx

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

    integer :: i_level,i_var_dep,petsc_stat,dealloc_stat

		! ================================================================
		! Work
		! ================================================================

    do i_level = 1,num_level
      do i_var_dep = 1,num_var_dep
        ! Remove values
        deallocate(mats_(i_var_dep,i_level)%vals,stat = dealloc_stat)
        if(dealloc_stat .ne. 0) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `mats_(...,...)%vals' failed.")
        end if

        ! Remove columns
        deallocate(mats_(i_var_dep,i_level)%cols,stat = dealloc_stat)
        if(dealloc_stat .ne. 0) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `mats_(...,...)%cols' failed.")
        end if

        ! Remove rows
        deallocate(mats_(i_var_dep,i_level)%rows,stat = dealloc_stat)
        if(dealloc_stat .ne. 0) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `mats_(...,...)%rows' failed.")
        end if
      end do ! `i_var_dep'
    end do ! `i_level'

    ! Free space occupied by matrices for linear equations
    deallocate(mats_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `mats_' failed.")
    end if

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_feval

  subroutine feval_read_arg_core_feval(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_prob:feval_read_arg_core_feval(...)': Read input arguments for problem and solver definitions
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHAR

    ! Project
    use mod_mpi,only : max_len_char_,mpi_world_,mpi_time_,mpi_sanity_check,mpi_exit_gracefully
    use mod_pf,only : num_level

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

    integer :: num_arg_read,io_stat,mpi_stat,alloc_stat
    character(len = max_len_char_) :: cmd_line_arg,namelist_arg,io_msg

    ! ================================================================
    ! Work
    ! ================================================================

    ! Define sweeper
    sweeper = 2 ! `1' for explicit, `2' for implicit, and `3' for IMEX
    use_LU = .true. ! LU factorization in implicit sweeper

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_feval)
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
          namelist_arg = "&arg_core_feval "//trim(adjustL(cmd_line_arg))//" /"
          read(namelist_arg,nml = arg_core_feval,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(sweeper,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(use_LU,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine feval_read_arg_core_feval

  subroutine feval_f1eval(c_encap_sol,time,level,c_ctx,c_encap_f1)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_feval:feval_f1eval(...)': Define explicit time stepper
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr

    ! External
    use pfasst,only : pfdp

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: level
    real(pfdp),intent(in) :: time
    type(c_ptr),value,intent(in) :: c_ctx,c_encap_sol,c_encap_f1

    ! ================================================================
    ! Locals
    ! ================================================================

    ! ================================================================
    ! Work
    ! ================================================================

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine feval_f1eval

  subroutine feval_f2eval(c_encap_sol,time,level,c_ctx,c_encap_f2)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_feval:feval_f2eval(...)': Define implicit time stepper
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
    use mod_ctx,only : num_var_dep,type_ctx,all_f_ctx_
    use mod_prob,only : prob,therm_diffu
    use mod_deriv,only : deriv_acc,deriv_wid,ord2_coeffs
    use mod_encap,only : type_encap

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: level
    real(pfdp),intent(in) :: time
    type(c_ptr),value,intent(in) :: c_ctx,c_encap_sol,c_encap_f2

    ! ================================================================
    ! Locals
    ! ================================================================

    type(type_encap),pointer :: f_encap_sol,f_encap_f2
    type(DMDALocalInfo),dimension(DMDA_LOCAL_INFO_SIZE) :: info
    integer :: step_space_x,step_space_y,stenc_col,petsc_stat,i_var_dep
    real(pfdp),dimension(:,:,:),pointer :: f_grid_arr_2d_sol,f_grid_arr_2d_f2

		! ================================================================
		! Work
		! ================================================================

    ! Access gird information
    call DMDAGetLocalInfoF90(all_f_ctx_(level)%f_ctx%mesh,info,petsc_stat)

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap_f2,f_encap_f2)

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap_sol,f_encap_sol)

    call DMGlobalToLocalBegin(all_f_ctx_(level)%f_ctx%mesh,f_encap_sol%u,INSERT_VALUES,all_f_ctx_(level)%f_ctx%pad,petsc_stat)
    call DMGlobalToLocalEnd(all_f_ctx_(level)%f_ctx%mesh,f_encap_sol%u,INSERT_VALUES,all_f_ctx_(level)%f_ctx%pad,petsc_stat)

    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      call DMDAVecGetArrayF90(all_f_ctx_(level)%f_ctx%mesh,f_encap_f2%u,f_grid_arr_2d_f2,petsc_stat)

      call DMDAVecGetArrayF90(all_f_ctx_(level)%f_ctx%mesh,all_f_ctx_(level)%f_ctx%pad,f_grid_arr_2d_sol,petsc_stat)

      select case(prob)
      case(1)
  !      !$OMP parallel do
        do step_space_y = info(DMDA_LOCAL_INFO_YS),info(DMDA_LOCAL_INFO_YS)+info(DMDA_LOCAL_INFO_YM)-1
          do step_space_x = info(DMDA_LOCAL_INFO_XS),info(DMDA_LOCAL_INFO_XS)+info(DMDA_LOCAL_INFO_XM)-1
            do i_var_dep = 0,num_var_dep-1
              f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y) = 0d0

              do stenc_col = 1,deriv_wid
                f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y) = f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y)+therm_diffu*ord2_coeffs(deriv_acc+1-stenc_col+1,2,level)*f_grid_arr_2d_sol(i_var_dep,step_space_x,step_space_y+deriv_wid-stenc_col+1)
              end do ! `stenc_col'

              do stenc_col = deriv_wid+1,deriv_acc
                f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y) = f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y)+therm_diffu*ord2_coeffs(3*deriv_wid+1-stenc_col+1,1,level)*f_grid_arr_2d_sol(i_var_dep,step_space_x+deriv_acc-stenc_col+1,step_space_y)
              end do ! `stenc_col'

              f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y) = f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y)+therm_diffu*sum(ord2_coeffs(deriv_wid+1,:,level))*f_grid_arr_2d_sol(i_var_dep,step_space_x,step_space_y)

              do stenc_col = deriv_acc+2,3*deriv_wid+1
                f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y) = f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y)+therm_diffu*ord2_coeffs(3*deriv_wid+1-stenc_col+1,1,level)*f_grid_arr_2d_sol(i_var_dep,step_space_x+deriv_acc-stenc_col+1,step_space_y)
              end do ! `stenc_col'

              do stenc_col = 3*deriv_wid+2,num_stenc_cols_
                f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y) = f_grid_arr_2d_f2(i_var_dep,step_space_x,step_space_y)+therm_diffu*ord2_coeffs(num_stenc_cols_-stenc_col+1,2,level)*f_grid_arr_2d_sol(i_var_dep,step_space_x,step_space_y+3*deriv_wid-stenc_col+1)
              end do ! `stenc_col'
            end do ! `i_var_dep'
          end do ! `step_space_x'
        end do ! `step_space_y'
  !      !$OMP end parallel do
      case(2)
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      case default
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      end select ! `prob'

      call DMDAVecRestoreArrayF90(all_f_ctx_(level)%f_ctx%mesh,all_f_ctx_(level)%f_ctx%pad,f_grid_arr_2d_sol,petsc_stat)

      call DMDAVecRestoreArrayF90(all_f_ctx_(level)%f_ctx%mesh,f_encap_f2%u,f_grid_arr_2d_f2,petsc_stat)
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine feval_f2eval

  subroutine feval_set_2d_res(info,grid_arr_2d_sol,grid_arr_2d_res,ctx,petsc_stat)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_feval:feval_set_2d_res(...)': Define local (w.r.t. MPI) residual function
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use omp_lib

    ! Project
    use mod_mpi,only : mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,type_ctx
    use mod_prob,only : prob,therm_diffu
    use mod_deriv,only : deriv_acc,deriv_wid,ord2_coeffs

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

		! ================================================================
		! Dummies
		! ================================================================
    
    integer,intent(inout) :: petsc_stat
    type(type_ctx),intent(in) :: ctx
    type(DMDALocalInfo),dimension(DMDA_LOCAL_INFO_SIZE),intent(inout) :: info
    real(pfdp),dimension(0:info(DMDA_LOCAL_INFO_DOF)-1,info(DMDA_LOCAL_INFO_GXS):info(DMDA_LOCAL_INFO_GXS)+info(DMDA_LOCAL_INFO_GXM)-1,info(DMDA_LOCAL_INFO_GYS):info(DMDA_LOCAL_INFO_GYS)+info(DMDA_LOCAL_INFO_GYM)-1),intent(in) :: grid_arr_2d_sol
    real(pfdp),dimension(0:info(DMDA_LOCAL_INFO_DOF)-1,info(DMDA_LOCAL_INFO_XS):info(DMDA_LOCAL_INFO_XS)+info(DMDA_LOCAL_INFO_XM)-1,info(DMDA_LOCAL_INFO_YS):info(DMDA_LOCAL_INFO_YS)+info(DMDA_LOCAL_INFO_YM)-1),intent(inout) :: grid_arr_2d_res

		! ================================================================
		! Locals
		! ================================================================

    integer :: step_space_x,step_space_y,stenc_col,i_var_dep

		! ================================================================
		! Work
		! ================================================================

    select case(prob)
    case(1)
!      !$OMP parallel do
      do step_space_y = info(DMDA_LOCAL_INFO_YS),info(DMDA_LOCAL_INFO_YS)+info(DMDA_LOCAL_INFO_YM)-1
        do step_space_x = info(DMDA_LOCAL_INFO_XS),info(DMDA_LOCAL_INFO_XS)+info(DMDA_LOCAL_INFO_XM)-1
          do i_var_dep = 0,num_var_dep-1
            grid_arr_2d_res(i_var_dep,step_space_x,step_space_y) = 0d0

            do stenc_col = 1,deriv_wid
              grid_arr_2d_res(i_var_dep,step_space_x,step_space_y) = grid_arr_2d_res(i_var_dep,step_space_x,step_space_y)-ctx%size_step_time*therm_diffu*ord2_coeffs(deriv_acc+1-stenc_col+1,2,ctx%level)*grid_arr_2d_sol(i_var_dep,step_space_x,step_space_y+deriv_wid-stenc_col+1)
            end do ! `stenc_col'

            do stenc_col = deriv_wid+1,deriv_acc
              grid_arr_2d_res(i_var_dep,step_space_x,step_space_y) = grid_arr_2d_res(i_var_dep,step_space_x,step_space_y)-ctx%size_step_time*therm_diffu*ord2_coeffs(3*deriv_wid+1-stenc_col+1,1,ctx%level)*grid_arr_2d_sol(i_var_dep,step_space_x+deriv_acc-stenc_col+1,step_space_y)
            end do ! `stenc_col'

            grid_arr_2d_res(i_var_dep,step_space_x,step_space_y) = grid_arr_2d_res(i_var_dep,step_space_x,step_space_y)+grid_arr_2d_sol(i_var_dep,step_space_x,step_space_y)-ctx%size_step_time*therm_diffu*sum(ord2_coeffs(deriv_wid+1,:,ctx%level))*grid_arr_2d_sol(i_var_dep,step_space_x,step_space_y)

            do stenc_col = deriv_acc+2,3*deriv_wid+1
              grid_arr_2d_res(i_var_dep,step_space_x,step_space_y) = grid_arr_2d_res(i_var_dep,step_space_x,step_space_y)-ctx%size_step_time*therm_diffu*ord2_coeffs(3*deriv_wid+1-stenc_col+1,1,ctx%level)*grid_arr_2d_sol(i_var_dep,step_space_x+deriv_acc-stenc_col+1,step_space_y)
            end do ! `stenc_col'

            do stenc_col = 3*deriv_wid+2,num_stenc_cols_
              grid_arr_2d_res(i_var_dep,step_space_x,step_space_y) = grid_arr_2d_res(i_var_dep,step_space_x,step_space_y)-ctx%size_step_time*therm_diffu*ord2_coeffs(num_stenc_cols_-stenc_col+1,2,ctx%level)*grid_arr_2d_sol(i_var_dep,step_space_x,step_space_y+3*deriv_wid-stenc_col+1)
            end do ! `stenc_col'
          end do ! `i_var_dep'
        end do ! `step_space_x`
      end do ! `step_space_y`
!      !$OMP end parallel do
    case(2)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
    end select ! `prob'

		! ================================================================
		! Exit
		! ================================================================

		return

  end subroutine feval_set_2d_res

  subroutine feval_set_2d_jac(info,grid_arr_2d_sol,jac,precon,ctx,petsc_stat)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_feval:feval_set_2d_jac(...)': Define local (w.r.t. MPI) Jacobian
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use petsc,only : Vec,SNES,INSERT_VALUES,MAT_FINAL_ASSEMBLY
    use omp_lib

    ! Project
    use mod_mpi,only : mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,type_ctx,all_f_ctx_
    use mod_prob,only : prob,therm_diffu
    use mod_deriv,only : deriv_acc,deriv_wid,ord2_coeffs

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

		! ================================================================
		! Dummies
		! ================================================================

    integer,intent(inout) :: petsc_stat
    type(type_ctx),intent(in) :: ctx
    type(Mat),intent(inout) :: jac,precon
    type(DMDALocalInfo),dimension(DMDA_LOCAL_INFO_SIZE),intent(inout) :: info
    real(pfdp),dimension(0:info(DMDA_LOCAL_INFO_DOF)-1,info(DMDA_LOCAL_INFO_GXS):info(DMDA_LOCAL_INFO_GXS)+info(DMDA_LOCAL_INFO_GXM)-1,info(DMDA_LOCAL_INFO_GYS):info(DMDA_LOCAL_INFO_GYS)+info(DMDA_LOCAL_INFO_GYM)-1),intent(in) :: grid_arr_2d_sol

		! ================================================================
		! Locals
		! ================================================================

    integer :: step_space_y,step_space_x,stenc_row,stenc_col,i_var_dep

		! ================================================================
		! Work
		! ================================================================

    select case(prob)
    case(1) ! 2D heat equation
!      !$OMP parallel do
      do step_space_y = info(DMDA_LOCAL_INFO_YS),info(DMDA_LOCAL_INFO_YS)+info(DMDA_LOCAL_INFO_YM)-1
        do step_space_x = info(DMDA_LOCAL_INFO_XS),info(DMDA_LOCAL_INFO_XS)+info(DMDA_LOCAL_INFO_XM)-1
          do i_var_dep = 0,num_var_dep-1
            do stenc_row = 1,num_stenc_rows_
              mats_(i_var_dep+1,ctx%level)%rows(MatStencil_i,stenc_row) = step_space_x
              mats_(i_var_dep+1,ctx%level)%rows(MatStencil_j,stenc_row) = step_space_y
              mats_(i_var_dep+1,ctx%level)%rows(MatStencil_c,stenc_row) = i_var_dep
            end do

            ! Is `{1,...,4}` for `deriv_acc = 8`
            do stenc_col = 1,deriv_wid
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_i,stenc_col) = step_space_x ! Is `step_space_x`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_j,stenc_col) = step_space_y+deriv_wid-stenc_col+1 ! Is `step_space_y+{4,...,1}` for `deriv_acc = 8`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_c,stenc_col) = i_var_dep

              ! Is `{1,...,4}` for `deriv_acc = 8`
              do stenc_row = 1,num_stenc_rows_
                mats_(i_var_dep+1,ctx%level)%vals(stenc_row,stenc_col) = -ctx%size_step_time*therm_diffu*ord2_coeffs(deriv_acc+1-stenc_col+1,2,ctx%level) ! Is `{9,...,6}` for `deriv_acc = 8`
              end do
            end do

            ! Is `{5,...,8}` for `deriv_acc = 8`
            do stenc_col = deriv_wid+1,deriv_acc
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_i,stenc_col) = step_space_x+deriv_acc-stenc_col+1 ! Is `step_space_x+{4,...,1}` for `deriv_acc = 8`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_j,stenc_col) = step_space_y ! Is `step_space_y`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_c,stenc_col) = i_var_dep

              ! Is `{5,...,8}` for `deriv_acc = 8`
              do stenc_row = 1,num_stenc_rows_
                mats_(i_var_dep+1,ctx%level)%vals(stenc_row,stenc_col) = -ctx%size_step_time*therm_diffu*ord2_coeffs(3*deriv_wid+1-stenc_col+1,1,ctx%level) ! Is `{9,...,6}` for `deriv_acc = 8`
              end do
            end do

            mats_(i_var_dep+1,ctx%level)%cols(MatStencil_i,deriv_acc+1) = step_space_x ! Is `step_space_x`
            mats_(i_var_dep+1,ctx%level)%cols(MatStencil_j,deriv_acc+1) = step_space_y ! Is `step_space_y`
            mats_(i_var_dep+1,ctx%level)%cols(MatStencil_c,deriv_acc+1) = i_var_dep

            ! Is `9` for `deriv_acc = 8`
            do stenc_row = 1,num_stenc_rows_
              mats_(i_var_dep+1,ctx%level)%vals(stenc_row,deriv_acc+1) = 1d0-ctx%size_step_time*therm_diffu*sum(ord2_coeffs(deriv_wid+1,:,ctx%level))
            end do

            ! Is `{10,...,13}` for `deriv_acc = 8`
            do stenc_col = deriv_acc+2,3*deriv_wid+1
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_i,stenc_col) = step_space_x+deriv_acc-stenc_col+1 ! Is `step_space_x+{-1,...,-4}` for `deriv_acc = 8`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_j,stenc_col) = step_space_y ! Is `step_space_y`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_c,stenc_col) = i_var_dep

              ! Is `{10,...,13}` for `deriv_acc = 8`
              do stenc_row = 1,num_stenc_rows_
                mats_(i_var_dep+1,ctx%level)%vals(stenc_row,stenc_col) = -ctx%size_step_time*therm_diffu*ord2_coeffs(3*deriv_wid+1-stenc_col+1,1,ctx%level) ! Is `{4,...,1}` for `deriv_acc = 8`
              end do
            end do

            ! Is `{14,...,17}` for `deriv_acc = 8`
            do stenc_col = 3*deriv_wid+2,num_stenc_cols_
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_i,stenc_col) = step_space_x ! Is `step_space_x`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_j,stenc_col) = step_space_y+3*deriv_wid-stenc_col+1 ! Is `step_space_y+{-1,...,-4}` for `deriv_acc = 8`
              mats_(i_var_dep+1,ctx%level)%cols(MatStencil_c,stenc_col) = i_var_dep

              ! Is `{14,...,17}` for `deriv_acc = 8`
              do stenc_row = 1,num_stenc_rows_
                mats_(i_var_dep+1,ctx%level)%vals(stenc_row,stenc_col) = -ctx%size_step_time*therm_diffu*ord2_coeffs(num_stenc_cols_-stenc_col+1,2,ctx%level) ! Is `{4,...,1}` for `deriv_acc = 8`
              end do
            end do

            call MatSetValuesStencil(jac,num_stenc_rows_,mats_(i_var_dep+1,ctx%level)%rows,num_stenc_cols_,mats_(i_var_dep+1,ctx%level)%cols,mats_(i_var_dep+1,ctx%level)%vals,INSERT_VALUES,petsc_stat)
          end do ! `dof'
        end do ! `step_space_x'
      end do ! `step_space_y'
!      !$OMP end parallel do

      call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,petsc_stat)
      call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,petsc_stat)
    case(2) ! 2D nonlinear diffusion
      call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
    end select ! `prob'

		! ================================================================
		! Exit
		! ================================================================

		return

  end subroutine feval_set_2d_jac

  subroutine feval_f2comp(c_encap_sol,time,size_step_time,c_encap_rhs,level,c_ctx,c_encap_f2)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_feval:feval_f2comp(...)': Define equation solver for implicit time stepper
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use petsc,only : DMDAVecGetArrayF90,DMDAVecRestoreArrayF90,INSERT_VALUES,PETSC_NULL_OBJECT
    use pfasst,only : pfdp

    ! Project
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,all_f_ctx_
    use mod_encap,only : type_encap
    use mod_prob,only : prob
    use mod_petsc,only : all_f_outer_

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: level
    real(pfdp),intent(in) :: size_step_time,time
    type(c_ptr),value,intent(in) :: c_ctx,c_encap_sol,c_encap_f2,c_encap_rhs

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: petsc_stat
    type(type_encap),pointer :: f_encap_sol,f_encap_rhs

    ! ================================================================
    ! Work
    ! ================================================================

    all_f_ctx_(level)%f_ctx%size_step_time = size_step_time

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap_sol,f_encap_sol)

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap_rhs,f_encap_rhs)

    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      select case(prob)
      case(1)
        call DMDASNESSetFunctionLocal(all_f_ctx_(level)%f_ctx%mesh,INSERT_VALUES,feval_set_2d_res,all_f_ctx_(level)%f_ctx,petsc_stat)
        call DMDASNESSetJacobianLocal(all_f_ctx_(level)%f_ctx%mesh,feval_set_2d_jac,all_f_ctx_(level)%f_ctx,petsc_stat)
        call SNESSolve(all_f_outer_(level),f_encap_rhs%u,f_encap_sol%u,petsc_stat)
      case(2)
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      case default
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      end select ! `prob'
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! Get "implicit" source evaluated at advanced time
    call feval_f2eval(c_encap_sol,time,level,c_ctx,c_encap_f2)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine feval_f2comp

end module mod_feval
