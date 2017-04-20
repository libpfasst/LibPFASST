! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

module mod_encap

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_encap': PFASST solution encapsulation
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! External
  use petsc,only : Vec
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

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! ================================================================
  ! Data
  ! ================================================================

  ! Encapsulation data
  type,public :: type_encap
    integer :: level
    type(Vec) :: u
  end type type_encap

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_encap,destroy_encap,encap_set_val,encap_copy,encap_get_norm,encap_set_norm,encap_pack,encap_unpack,encap_get_aXpY,encap_print,encap_create_pf,encap_destroy_pf

contains

  subroutine create_encap(c_encap,level,kind,nVars,shape,c_ctx)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:create_encap(...)': Initialize PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer,c_loc

    ! Project
    use mod_mpi,only : mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,all_f_ctx_

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: level,kind,nVars
    integer,dimension(:),intent(in) :: shape
    type(c_ptr),value,intent(in) :: c_ctx
    type(c_ptr),intent(inout) :: c_encap

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: alloc_stat,petsc_stat
    type(type_encap),pointer :: f_encap

    ! ================================================================
    ! Work
    ! ================================================================

    ! Encapsulation data
    allocate(f_encap,stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `f_encap' failed.")
    end if ! `alloc_stat .ne. 0'

    ! PFASST level identifier
    f_encap%level = level

    ! Create global (w.r.t. MPI) ghost-free PETSc vector
    call DMCreateGlobalVector(all_f_ctx_(f_encap%level)%f_ctx%mesh,f_encap%u,petsc_stat)

    ! Initialize C-pointer to PFASST encapsulation
    c_encap = c_loc(f_encap)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_encap

  subroutine destroy_encap(c_encap)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:destroy_encap(...)': Destroy PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! Project
    use mod_mpi,only : mpi_exit_gracefully

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

    integer :: dealloc_stat,petsc_stat
    type(type_encap),pointer :: f_encap

    ! ================================================================
    ! Work
    ! ================================================================

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap,f_encap)

    ! Create global (w.r.t. MPI) ghost-free PETSc vector
    call VecDestroy(f_encap%u,petsc_stat)

    ! Encapsulation data
    deallocate(f_encap,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `f_encap' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_encap

  subroutine encap_set_val(c_encap,val,flags)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_set_val(...)': Assign values to PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use pfasst,only : pfdp

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in),optional :: flags
    real(pfdp),intent(in) :: val
    type(c_ptr),value,intent(in) :: c_encap

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: petsc_stat
    type(type_encap),pointer :: f_encap

    ! ================================================================
    ! Work
    ! ================================================================

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap,f_encap)

    ! Assign value to PETSc vector
    call VecSet(f_encap%u,val,petsc_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_set_val

  subroutine encap_copy(c_encap_dst,c_encap_src,flags)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_copy(...)': Copy PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in),optional :: flags
    type(c_ptr),value,intent(in) :: c_encap_src,c_encap_dst

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: petsc_stat
    type(type_encap),pointer :: f_encap_src,f_encap_dst

    ! ================================================================
    ! Work
    ! ================================================================

    ! Access PFASST encapsulations
    call c_f_pointer(c_encap_src,f_encap_src)
    call c_f_pointer(c_encap_dst,f_encap_dst)

    ! Copy PETSc vectors
    call VecCopy(f_encap_src%u,f_encap_dst%u,petsc_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_copy

  function encap_get_norm(c_encap) result(norm)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_get_norm(...)': Get the norm of PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use pfasst,only : pfdp

    ! Project
    use mod_ctx,only : norm_

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

    real(pfdp) :: norm
    type(type_encap),pointer :: f_encap

    ! ================================================================
    ! Work
    ! ================================================================
    
    ! Get norm of PFASST solution data for each sol component
    call encap_set_norm(c_encap)

    ! Get norm of PFASST solution data over all sol components
    norm = maxVal(norm_)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end function encap_get_norm

  subroutine encap_set_norm(c_encap)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_set_norm(...)': Get the norm of all sol components of PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use petsc,only : NORM_INFINITY

    ! Project
    use mod_ctx,only : norm_

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

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap,f_encap)
    
    ! Get norm of PETSc vector
    call VecNorm(f_encap%u,NORM_INFINITY,norm_,petsc_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_set_norm

  subroutine encap_pack(flat_arr,c_encap)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_pack(...)': Pack PFASST encapsulations to flat array
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use petsc,only : DMDAVecGetArrayF90,DMDAVecRestoreArrayF90
    use pfasst,only : pfdp

    ! Project
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,all_f_ctx_

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    real(pfdp),dimension(:),intent(out) :: flat_arr
    type(c_ptr),value,intent(in) :: c_encap

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: i_var_dep,petsc_stat
    type(type_encap),pointer :: f_encap
    real(pfdp),dimension(:,:,:),pointer :: f_grid_arr_2d

    ! ================================================================
    ! Work
    ! ================================================================

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap,f_encap)

    ! Dimension-dependent pack
    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      ! Read from PETSc vector
      call DMDAVecGetArrayF90(all_f_ctx_(f_encap%level)%f_ctx%mesh,f_encap%u,f_grid_arr_2d,petsc_stat)

      ! Assign to flat PFASST container `flat_arr'
      flat_arr = reshape(f_grid_arr_2d,(/num_var_dep*all_f_ctx_(f_encap%level)%f_ctx%count(1)/))

      ! Disassociate from PETSc vector
      call DMDAVecRestoreArrayF90(all_f_ctx_(f_encap%level)%f_ctx%mesh,f_encap%u,f_grid_arr_2d,petsc_stat)
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_pack

  subroutine encap_unpack(c_encap,flat_arr)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_unpack(...)': Revert packing of PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use petsc,only : DMDAVecGetArrayF90,DMDAVecRestoreArrayF90
    use pfasst,only : pfdp

    ! Project
    use mod_mpi,only : num_dim_space,mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,all_f_ctx_

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    real(pfdp),dimension(:),intent(in) :: flat_arr
    type(c_ptr),value,intent(in) :: c_encap

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: i_var_dep,petsc_stat
    type(type_encap),pointer :: f_encap
    real(pfdp),dimension(:,:,:),pointer :: f_grid_arr_2d

    ! ================================================================
    ! Work
    ! ================================================================

    ! Access PFASST encapsulation
    call c_f_pointer(c_encap,f_encap)

    ! Dimension-dependent unpack
    select case(num_dim_space)
    case(1)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    case(2)
      ! Associate with PETSc vector
      call DMDAVecGetArrayF90(all_f_ctx_(f_encap%level)%f_ctx%mesh,f_encap%u,f_grid_arr_2d,petsc_stat)

      ! Get from flat PFASST container `flat_arr'
      f_grid_arr_2d = reshape(flat_arr,shape(f_grid_arr_2d))

      ! Give to PETSc vector
      call DMDAVecRestoreArrayF90(all_f_ctx_(f_encap%level)%f_ctx%mesh,f_encap%u,f_grid_arr_2d,petsc_stat)
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
    end select ! `num_dim_space'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_unpack

  subroutine encap_get_aXpY(c_encap_y,a,c_encap_x,flags)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_get_aXpY(...)': Calculate the `aX+Y' of PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr,c_f_pointer

    ! External
    use pfasst,only : pfdp

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in),optional :: flags
    real(pfdp),intent(in) :: a
    type(c_ptr),value,intent(in) :: c_encap_x,c_encap_y

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: petsc_stat
    type(type_encap),pointer :: f_encap_x,f_encap_y

    ! ================================================================
    ! Work
    ! ================================================================

    ! Access PFASST encapsulations
    call c_f_pointer(c_encap_x,f_encap_x)
    call c_f_pointer(c_encap_y,f_encap_y)

    ! Assign `aXpY' by replacement
    call VecAXPY(f_encap_y%u,a,f_encap_x%u,petsc_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_get_aXpY

  subroutine encap_print(c_encap)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_print(...)': Print PFASST encapsulations
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr

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

    ! ================================================================
    ! Work
    ! ================================================================

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_print

  subroutine encap_create_pf(encap_pf)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_create_pf(...)': Create PFASST encapsulation containers
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use pfasst,only : pf_encap_t

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    type(pf_encap_t),intent(out) :: encap_pf

    ! ================================================================
    ! Locals
    ! ================================================================

    ! ================================================================
    ! Work
    ! ================================================================

    encap_pf%create => create_encap
    encap_pf%destroy => destroy_encap
    encap_pf%setVal => encap_set_val
    encap_pf%copy => encap_copy
    encap_pf%norm => encap_get_norm
    encap_pf%pack => encap_pack
    encap_pf%unpack => encap_unpack
    encap_pf%aXpY => encap_get_aXpY
    encap_pf%ePrint => encap_print

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine encap_create_pf

  subroutine encap_destroy_pf(pf_encap)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_encap:encap_destroy_pf(...)': Destroy PFASST encapsulation containers
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use pfasst,only : pf_encap_t

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    type(pf_encap_t),intent(inout) :: pf_encap

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

  end subroutine encap_destroy_pf

end module mod_encap
