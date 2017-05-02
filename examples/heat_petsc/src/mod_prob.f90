! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

module mod_prob

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_prob:': Define problem and solver specifics (e.g. for [linear] KSP or [scalable] nonlinear [equation] solver, i.e. SNES)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! Fortran
  use iso_c_binding,only : c_ptr

  ! External
  use petsc,only : Mat,KSP,SNES
  use pfasst,only : pfdp

  ! Project
  use mod_mpi,only : max_len_char_

  ! ================================================================
  ! Safety
  ! ================================================================

  implicit none

  ! ================================================================ 
  ! Parameters
  ! ================================================================ 

  ! Numerical constants
  real(pfdp),parameter,private :: pi_ = 4d0*atan(1d0) ! Constant `pi'
  real(pfdp),parameter,private :: ti_ = 8d0*atan(1d0) ! Two times `pi'

  ! Available problems to solve (`1' for 2D heat equation, `2' for 1D nonlinear diffusion)
  integer,dimension(2),parameter,private :: avail_probs_ = (/1,2/)

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! Problem to solve (`1' for 2D heat equation, `2' for 1D nonlinear diffusion)
  integer,save,public :: prob 
  ! Thermal diffussivity for `prob = 1'
  real(pfdp),save,public :: therm_diffu 

  ! List of arguments
  namelist /arg_core_prob/ prob,therm_diffu

  ! ================================================================
  ! Data
  ! ================================================================

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_prob,destroy_prob,prob_set_ini,prob_get_err,prob_set_exa
  private :: prob_read_arg_core_prob

contains

  subroutine create_prob(io_unit,name_file_input_pfasst,name_file_input_petsc,pfasst_pf)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_pf:create_pf(...)': Allocate memory and define objects for problems and solver objects
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION
    use pfasst,only : pf_pfasst_t

    ! Project
    use mod_mpi,only : max_len_char_,mpi_exit_gracefully,mpi_world_
    use mod_pf,only : num_level
    use mod_ctx,only : all_f_ctx_
    use mod_encap,only : create_encap

    ! ================================================================ 
    ! Safety
    ! ================================================================ 

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: io_unit
    character(len = *),intent(in) :: name_file_input_pfasst,name_file_input_petsc
    type(pf_pfasst_t),intent(in) :: pfasst_pf

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: num_arg_read,io_stat,mpi_stat,alloc_stat,i_level,petsc_stat
    character(len = max_len_char_) :: cmd_line_arg,namelist_arg,io_msg

    ! ================================================================
    ! Work
    ! ================================================================

    ! Get input arguments
    call prob_read_arg_core_prob(io_unit,name_file_input_pfasst)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_prob

  subroutine destroy_prob()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_prob:destroy_prob(...)`: Free occupied memory
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================ 
    ! Modules
    ! ================================================================ 

    ! Project
    use mod_mpi,only : mpi_exit_gracefully
    use mod_pf,only : num_level
    use mod_ctx,only : all_f_ctx_
    use mod_encap,only : destroy_encap

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

    integer :: dealloc_stat,i_level,petsc_stat

    ! ================================================================
    ! Work
    ! ================================================================

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_prob

  subroutine prob_read_arg_core_prob(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_prob:prob_read_arg_core_prob(...)': Read input arguments for problem and solver definitions
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

    ! Problem choice
    prob = 1 ! Problem to solve (`1' for 2D heat equation, `2' for 1D nonlinear diffusion)

    ! Problem parameters
    therm_diffu = 6.25d-2 ! Thermal diffussivity for `prob = 1'

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_prob)
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
          namelist_arg = "&arg_core_prob "//trim(adjustL(cmd_line_arg))//" /"
          read(namelist_arg,nml = arg_core_prob,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(prob,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(therm_diffu,1,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! Sanity check for problem to solve
    call mpi_sanity_check(avail_probs_,prob)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine prob_read_arg_core_prob

  subroutine prob_set_ini(c_encap)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_prob:prob_set_ini(...)': Assign initial data to PFASST encapsulation
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================ 
    ! Modules
    ! ================================================================ 

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! Project
    use mod_pf,only : orig_dom_time
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully
    use mod_ctx,only : all_f_ctx_
    use mod_encap,only : type_encap

    ! ================================================================ 
    ! Safety
    ! ================================================================ 

    implicit none

    ! ================================================================ 
    ! Dummies
    ! ================================================================ 

    type(c_ptr),value,intent(in) :: c_encap

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: petsc_stat
    type(type_encap),pointer :: f_encap

    ! ================================================================
    ! Work
    ! ================================================================

    call c_f_pointer(c_encap,f_encap)
    
    ! Solution initialization depends on spatial dimensionality
    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      select case(prob)
      case(1)
        call prob_set_exa(orig_dom_time,f_encap%level,f_encap%u)
      case(2)
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      case default
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      end select ! `prob'
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine prob_set_ini

  function prob_get_err(time,c_encap) result(err)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_prob:prob_get_err(...)': Get the err of PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use petsc,only : NORM_INFINITY
    use pfasst,only : pfdp

    ! Project
    use mod_ctx,only : all_f_ctx_,norm_
    use mod_encap,only : type_encap

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    real(pfdp),intent(in) :: time
    type(c_ptr),value,intent(in) :: c_encap

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: petsc_stat
    real(pfdp) :: err
    type(type_encap),pointer :: f_encap

    ! ================================================================
    ! Work
    ! ================================================================

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap,f_encap)

    ! Get err of PFASST solution data for each sol component
    call prob_set_exa(time,f_encap%level,all_f_ctx_(f_encap%level)%f_ctx%ref)

    ! Get difference between solution and exact counterpart
    call VecAXPY(all_f_ctx_(f_encap%level)%f_ctx%ref,-1d0,f_encap%u,petsc_stat)
    
    ! Get err of PFASST solution data over all sol components
    call VecNorm(all_f_ctx_(f_encap%level)%f_ctx%ref,NORM_INFINITY,norm_,petsc_stat)

    ! Define error
    err = maxVal(norm_)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end function prob_get_err

  subroutine prob_set_exa(time,level,ref)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_prob:prob_set_exa(...)': Assign exact solution values
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================ 
    ! Modules
    ! ================================================================ 

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use omp_lib
    use petsc,only : Vec,DMDAVecGetArrayF90,DMDAVecRestoreArrayF90


    ! Project
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully
    use mod_ctx,only : size_dom_space,size_step_space,orig_dom_space,num_var_dep,all_f_ctx_
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
    type(Vec),intent(inout) :: ref

    ! ================================================================
    ! Locals
    ! ================================================================

    type(type_encap),pointer :: f_encap
    real(pfdp),dimension(:,:,:),pointer :: f_grid_arr_2d
    integer :: step_space_x,step_space_y,petsc_stat,var_dep
    type(DMDALocalInfo),dimension(DMDA_LOCAL_INFO_SIZE) :: info

    ! ================================================================
    ! Work
    ! ================================================================

    call DMDAGetLocalInfoF90(all_f_ctx_(level)%f_ctx%mesh,info,petsc_stat)

    ! Solution initialization depends on spatial dimensionality
    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      select case(prob)
      case(1)
        ! Associate with PETSc vector
        call DMDAVecGetArrayF90(all_f_ctx_(level)%f_ctx%mesh,ref,f_grid_arr_2d,petsc_stat)

        !$OMP parallel do
        do step_space_y = info(DMDA_LOCAL_INFO_YS),info(DMDA_LOCAL_INFO_YS)+info(DMDA_LOCAL_INFO_YM)-1
          do step_space_x = info(DMDA_LOCAL_INFO_XS),info(DMDA_LOCAL_INFO_XS)+info(DMDA_LOCAL_INFO_XM)-1
            do var_dep = 0,num_var_dep-1
              ! Define exact value
              f_grid_arr_2d(var_dep,step_space_x,step_space_y) = sin(ti_*(orig_dom_space(1)+step_space_x*size_step_space(1,level))/size_dom_space(1))*exp(-therm_diffu*time*(ti_/size_dom_space(1))**2d0)*sin(ti_*(orig_dom_space(2)+step_space_y*size_step_space(2,level))/size_dom_space(2))*exp(-therm_diffu*time*(ti_/size_dom_space(2))**2d0)
            end do ! `var_dep'
          end do ! `step_space_x'
        end do ! `step_space_y'
        !$OMP end parallel do

        ! Give to PETSc vector
        call DMDAVecRestoreArrayF90(all_f_ctx_(level)%f_ctx%mesh,ref,f_grid_arr_2d,petsc_stat)
      case(2)
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      case default
        call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
      end select ! `prob'
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine prob_set_exa

end module mod_prob
