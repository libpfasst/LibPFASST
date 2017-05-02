! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

module mod_petsc

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_petsc:': Define problem and solver specifics (e.g. for [linear] KSP or [scalable] nonlinear [equation] solver, i.e. SNES)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! Fortran
  use iso_c_binding,only : c_ptr

  ! External
  use petsc,only : Mat,KSP,SNES
  use pfasst,only : pfdp,pf_state_t

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

  ! Define specifics for nonlinear PETSc SNES problem solvers
  real(pfdp),dimension(:),allocatable,save,public :: tol_abs_outer ! Absolute tolerances
  real(pfdp),dimension(:),allocatable,save,public :: tol_rel_outer ! Relative tolerances
  real(pfdp),dimension(:),allocatable,save,public :: tol_step_outer ! Correction step tolerances, i.e. `||dx|| < tol*||x||'
  integer,dimension(:),allocatable,save,public :: max_iter_outer ! Maxiumum SNES (e.g. Newton) iterations
  integer,dimension(:),allocatable,save,public :: max_func_eval_outer ! Maximum residual function evaluations

  ! Define specifics for nonlinear PETSc SNES problem solvers
  logical,save,public :: monitor_outer ! Print or save solver statistics
  logical,save,public :: print_outer ! Print solver statistics
  logical,save,public :: save_outer ! Save solver statistics
  integer,save,public :: print_freq_outer ! Steps between output
  integer,save,public :: save_freq_outer ! Steps between output

  ! Define specifics for linear PETSc KPS problem solvers
  real(pfdp),dimension(:),allocatable,save,public :: tol_rel_inner ! Relative tolerances
  real(pfdp),dimension(:),allocatable,save,public :: tol_abs_inner ! Absolute tolerances
  real(pfdp),dimension(:),allocatable,save,public :: tol_div_inner ! Divergence tolerances
  integer,dimension(:),allocatable,save,public :: max_iter_inner ! Maximum iterations

  ! Define specifics for linear PETSc KPS problem solvers
  logical,save,public :: monitor_inner ! Print or save solver statistics
  logical,save,public :: print_inner ! Print solver statistics
  logical,save,public :: save_inner ! Save solver statistics
  integer,save,public :: print_freq_inner ! Steps between output
  integer,save,public :: save_freq_inner ! Steps between output

  ! List of arguments
  namelist /arg_core_petsc/ monitor_outer,print_outer,print_freq_outer,save_outer,save_freq_outer,tol_abs_outer,tol_rel_outer,tol_step_outer,max_iter_outer,max_func_eval_outer,monitor_inner,print_inner,print_freq_inner,save_inner,save_freq_inner,tol_rel_inner,tol_abs_inner,tol_div_inner,max_iter_inner

  ! ================================================================
  ! Data
  ! ================================================================

  ! PFASST state
  type(pf_state_t),pointer :: state_pf_

  ! Problem solvers (a `num_level' array)
  type(SNES),dimension(:),pointer,public :: all_f_outer_ ! SNES solvers

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_petsc,destroy_petsc
  private :: petsc_read_arg_core_petsc,outer_viewer,inner_viewer

contains

  subroutine create_petsc(io_unit,name_file_input_pfasst,name_file_input_petsc,pfasst_pf)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_pf:create_pf(...)': Allocate memory and define objects for problems and solver objects
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION
    use petsc,only : PETSC_COMM_WORLD,PETSC_NULL_OBJECT,PETSC_TRUE,PETSC_NULL_FUNCTION,SNESMonitorDefault,PETSC_NULL_CHARACTER
    use pfasst,only : pf_pfasst_t

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : max_len_char_,mpi_exit_gracefully,mpi_world_,mpi_space_
    use mod_ctx,only : num_var_dep,all_f_ctx_
    use mod_prob,only : prob
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
    type(pf_pfasst_t),target,intent(in) :: pfasst_pf

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: num_arg_read,io_stat,mpi_stat,alloc_stat,dealloc_stat,i_level,petsc_stat
    character(len = max_len_char_) :: cmd_line_arg,namelist_arg,io_msg
    type(KSP),dimension(:),pointer :: all_inner_

    ! ================================================================
    ! Work
    ! ================================================================

    state_pf_ => pfasst_pf%state

    ! Get input arguments
    call petsc_read_arg_core_petsc(io_unit,name_file_input_pfasst)

    ! Prepare array of SNES problem solvers
    allocate(all_f_outer_(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_f_outer_(...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Prepare array of KSP problem solvers
    allocate(all_inner_(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_inner_(...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Create SNES problem solvers
    do i_level = 1,num_level
      ! Get input parameters
      call petscOptionsInsertFile(PETSC_COMM_WORLD,PETSC_NULL_OBJECT,trim(adjustL(name_file_input_petsc)),PETSC_TRUE,petsc_stat)

      ! Create SNES and set outer tolerances
      call SNESCreate(PETSC_COMM_WORLD,all_f_outer_(i_level),petsc_stat)
      call SNESSetDM(all_f_outer_(i_level),all_f_ctx_(i_level)%f_ctx%mesh,petsc_stat)
      call SNESSetTolerances(all_f_outer_(i_level),tol_abs_outer(i_level),tol_rel_outer(i_level),tol_step_outer(i_level),max_iter_outer(i_level),max_func_eval_outer(i_level),petsc_stat)

      ! Set inner tolerances
      call SNESGetKSP(all_f_outer_(i_level),all_inner_(i_level),petsc_stat)
      call KSPSetTolerances(all_inner_(i_level),tol_rel_inner(i_level),tol_abs_inner(i_level),tol_div_inner(i_level),max_iter_inner(i_level),petsc_stat)

      ! Define monitoring routines
      if(monitor_outer) then
        call SNESMonitorSet(all_f_outer_(i_level),outer_viewer,pfasst_pf%levels(i_level),PETSC_NULL_FUNCTION,petsc_stat)
      end if ! `monitor_outer'

      ! Define monitoring routines
      if(monitor_inner) then
        call KSPMonitorSet(all_inner_(i_level),inner_viewer,pfasst_pf%levels(i_level),PETSC_NULL_FUNCTION,petsc_stat)
      end if ! `monitor_inner'

      ! Finalize solver setup
      call KSPSetFromOptions(all_inner_(i_level),petsc_stat)
      call SNESSetFromOptions(all_f_outer_(i_level),petsc_stat)

      ! Ignore outer vector operations and use inner KSP only (e.g. for a linear problem)
      if(prob .eq. 1) then
        call SNESSetType(all_f_outer_(i_level),SNESKSPONLY,petsc_stat)
      end if ! `prob .eq. 1'
    end do ! `i_level'

    ! If KSP was created only for inner tolerances then destroy it again
    deallocate(all_inner_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_inner_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_petsc

  subroutine destroy_petsc()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_petsc:destroy_petsc(...)`: Free occupied memory
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================ 
    ! Modules
    ! ================================================================ 

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : mpi_exit_gracefully
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

    ! Call PETSc's SNES destructor
    do i_level = 1,num_level
      call SNESDestroy(all_f_outer_(i_level),petsc_stat)
    end do ! `i_level'

    ! Destroy array of SNES problem solvers
    deallocate(all_f_outer_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_f_outer_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Maximum iterations
    deallocate(max_iter_inner,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `max_iter_inner' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Divergence tolerances
    deallocate(tol_div_inner,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `tol_div_inner' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Absolute tolerances
    deallocate(tol_abs_inner,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `tol_abs_inner' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Relative tolerances
    deallocate(tol_rel_inner,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `tol_rel_inner' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Maximum residual function evaluations
    deallocate(max_func_eval_outer,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `max_func_eval_outer' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Maxiumum SNES (e.g. Newton) iterations
    deallocate(max_iter_outer,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `max_iter_outer' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Correction step tolerances, i.e. `||dx|| < tol*||x||'
    deallocate(tol_step_outer,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `tol_step_outer' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Relative tolerances
    deallocate(tol_rel_outer,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `tol_rel_outer' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Absolute tolerances
    deallocate(tol_abs_outer,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `tol_abs_outer' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_petsc

  subroutine petsc_read_arg_core_petsc(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_petsc:petsc_read_arg_core_petsc(...)': Read input arguments for problem and solver definitions
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHAR

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : max_len_char_,mpi_world_,mpi_time_,mpi_sanity_check,mpi_exit_gracefully

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

    ! Define specifics for nonlinear PETSc SNES problem solvers
    monitor_outer = .true. ! Print or save solver statistics
    print_outer = .true. ! Print solver statistics
    save_outer = .true. ! Save solver statistics
    print_freq_outer = 1 ! Steps between output
    save_freq_outer = 1 ! Steps between output

    ! Absolute tolerances
    allocate(tol_abs_outer(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `tol_abs_outer(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for nonlinear PETSc SNES problem solvers
    tol_abs_outer = 1e-15

    ! Relative tolerances
    allocate(tol_rel_outer(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `tol_rel_outer(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for nonlinear PETSc SNES problem solvers
    tol_rel_outer = 1e-15

    ! Correction step tolerances, i.e. `||dx|| < tol*||x||'
    allocate(tol_step_outer(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `tol_step_outer(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for nonlinear PETSc SNES problem solvers
    tol_step_outer = 1e-15

    ! Maxiumum SNES (e.g. Newton) iterations
    allocate(max_iter_outer(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `max_iter_outer(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for nonlinear PETSc SNES problem solvers
    max_iter_outer = 10000

    ! Maximum residual function evaluations
    allocate(max_func_eval_outer(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `max_func_eval_outer(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for nonlinear PETSc SNES problem solvers
    max_func_eval_outer = 10000

    ! Define specifics for linear PETSc KPS problem solvers
    monitor_inner = .true. ! Print or save solver statistics
    print_inner = .true. ! Print solver statistics
    save_inner = .true. ! Save solver statistics
    print_freq_inner = 1 ! Steps between output
    save_freq_inner = 1 ! Steps between output

    ! Relative tolerances
    allocate(tol_rel_inner(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `tol_rel_inner(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for linear PETSc KPS problem solvers
    tol_rel_inner = 1e-15

    ! Absolute tolerances
    allocate(tol_abs_inner(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `tol_abs_inner(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for linear PETSc KPS problem solvers
    tol_abs_inner = 1e-15

    ! Divergence tolerances
    allocate(tol_div_inner(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `tol_div_inner(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for linear PETSc KPS problem solvers
    tol_div_inner = 1e+4

    ! Maximum iterations
    allocate(max_iter_inner(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `max_iter_inner(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Define specifics for linear PETSc KPS problem solvers
    max_iter_inner = 10000

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_petsc)
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
          namelist_arg = "&arg_core_petsc "//trim(adjustL(cmd_line_arg))//" /"
          read(namelist_arg,nml = arg_core_petsc,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(monitor_outer,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_outer,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_outer,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_freq_outer,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_freq_outer,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(tol_abs_outer,num_level,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(tol_rel_outer,num_level,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(tol_step_outer,num_level,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(max_iter_outer,num_level,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(max_func_eval_outer,num_level,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(monitor_inner,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_inner,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_inner,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_freq_inner,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_freq_inner,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(tol_rel_inner,num_level,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(tol_abs_inner,num_level,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(tol_div_inner,num_level,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(max_iter_inner,num_level,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine petsc_read_arg_core_petsc

  subroutine outer_viewer(in_outer,outer_iter,outer_res,level_pf,petsc_stat)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_petsc:outer_viewer(...)': Define KSP solvers
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use petsc,only : PETSC_NULL_OBJECT,PETSC_TRUE,PETSC_COMM_WORLD,PETSC_NULL_CHARACTER
    use pfasst,only : pf_level_t

    ! Project
    use mod_pf,only : num_level,num_step_time
    use mod_mpi,only : max_len_char_,mpi_space_,mpi_time_,mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,type_ctx
    use mod_hook,only : num_outer_,name_outer_,time_outer_,data_outer_,base_io_unit_outer_,hook_print,hook_save

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: outer_iter,petsc_stat
    real(pfdp),intent(in) :: outer_res
    type(pf_level_t),intent(in) :: level_pf
    type(SNES),intent(in) :: in_outer

    ! ================================================================
    ! Locals
    ! ================================================================

    logical :: is_open,does_exist
    integer :: io_unit,io_stat,i_outer
    character(len = max_len_char_) :: name_file

    ! ================================================================
    ! Work
    ! ================================================================

    if((print_outer .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,print_freq_outer) .eq. 0))) .or. (save_outer .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,save_freq_outer) .eq. 0)))) then
      ! Time values
      write(time_outer_(1),"(1es24.16e3)") state_pf_%t0
      write(time_outer_(2),"(1es24.16e3)") state_pf_%t0

      ! Field values
      write(data_outer_(1),"(1i24)") outer_iter
      write(data_outer_(2),"(1es24.16e3)") outer_res

      ! Print data
      if(print_outer) then
        call hook_print(num_outer_,name_outer_,time_outer_,data_outer_,state_pf_%hook,state_pf_%step,level_pf%level-1,state_pf_%sweep-1,state_pf_%iter-1)
      end if ! `print_outer'

      ! Save data
      if(save_outer) then
        call hook_save(num_outer_,name_outer_,time_outer_,data_outer_,base_io_unit_outer_,save_freq_outer,state_pf_%hook,state_pf_%step,level_pf%level-1,state_pf_%sweep-1,state_pf_%iter-1)
      end if ! `save_outer'
    end if ! `(print_outer .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,print_freq_outer) .eq. 0))) .or. (save_outer .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,save_freq_outer) .eq. 0)))'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine outer_viewer

  subroutine inner_viewer(in_inner,inner_iter,inner_res,level_pf,petsc_stat)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_petsc:inner_viewer(...)': Define KSP solvers
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use petsc,only : PETSC_NULL_OBJECT,PETSC_TRUE,PETSC_COMM_WORLD,PETSC_NULL_CHARACTER
    use pfasst,only : pf_level_t

    ! Project
    use mod_pf,only : num_level,num_step_time
    use mod_mpi,only : max_len_char_,mpi_space_,mpi_time_,mpi_exit_gracefully
    use mod_ctx,only : num_var_dep,type_ctx
    use mod_hook,only : num_inner_,name_inner_,time_inner_,data_inner_,base_io_unit_inner_,hook_print,hook_save

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: inner_iter,petsc_stat
    real(pfdp),intent(in) :: inner_res
    type(pf_level_t),intent(in) :: level_pf
    type(KSP),intent(in) :: in_inner

    ! ================================================================
    ! Locals
    ! ================================================================

    logical :: is_open,does_exist
    integer :: io_unit,io_stat,i_inner
    character(len = max_len_char_) :: name_file

    ! ================================================================
    ! Work
    ! ================================================================

    if((print_inner .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,print_freq_inner) .eq. 0))) .or. (save_inner .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,save_freq_inner) .eq. 0)))) then
      ! Time values
      write(time_inner_(1),"(1es24.16e3)") state_pf_%t0
      write(time_inner_(2),"(1es24.16e3)") state_pf_%t0

      ! Field values
      write(data_inner_(1),"(1i24)") inner_iter
      write(data_inner_(2),"(1es24.16e3)") inner_res

      ! Print data
      if(print_inner) then
        call hook_print(num_inner_,name_inner_,time_inner_,data_inner_,state_pf_%hook,state_pf_%step,level_pf%level-1,state_pf_%sweep-1,state_pf_%iter-1)
      end if ! `print_inner'

      ! Save data
      if(save_inner) then
        call hook_save(num_inner_,name_inner_,time_inner_,data_inner_,base_io_unit_inner_,save_freq_inner,state_pf_%hook,state_pf_%step,level_pf%level-1,state_pf_%sweep-1,state_pf_%iter-1)
      end if ! `save_inner'
    end if ! `(print_inner .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,print_freq_inner) .eq. 0))) .or. (save_inner .and. ((state_pf_%step .eq. 0) .or. (mod(state_pf_%step+1,save_freq_inner) .eq. 0)))'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine inner_viewer

end module mod_petsc
