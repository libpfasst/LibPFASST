! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

module mod_hook

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_hook': Output routines
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! Project
  use mod_mpi,only : max_len_char_

  ! ================================================================
  ! Safety
  ! ================================================================

  implicit none

  ! ================================================================ 
  ! Parameters
  ! ================================================================ 

  ! Hook names
  integer,parameter,public :: num_hook_ = 18 ! Number of hooks
  ! Outer sols
  integer,parameter,public :: num_outer_ = 2 ! Number of SNES sols
  ! Inner sols
  integer,parameter,public :: num_inner_ = 2 ! Number of KSP sols
  ! Residual sols
  integer,parameter,public :: num_res_ = 3 ! Number of residual sols
  ! Error sols
  integer,parameter,public :: num_err_ = 2 ! Number of error sols

  ! Field names
  integer,parameter,public :: max_len_name_ = 10 ! Length of longest sol name
  ! Time values
  integer,parameter,public :: max_len_time_ = 24 ! Length of time strings
  ! Data values
  integer,parameter,public :: max_len_data_ = 24 ! Length of data strings

  ! Lower bound on output units
  integer,parameter,public :: base_io_unit_ = 111 ! Base output unit

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! Hooks to attach for output
  logical,save,public :: ante_pred ! Before predictor
  logical,save,public :: post_pred ! After predictor
  logical,save,public :: ante_iter ! Before iteration
  logical,save,public :: post_iter ! After iteration
  logical,save,public :: ante_sweep ! Before sweep
  logical,save,public :: post_sweep ! After sweep
  logical,save,public :: ante_step ! Before step
  logical,save,public :: post_step ! After step
  logical,save,public :: ante_prolo ! Before prolongation
  logical,save,public :: post_prolo ! After prolongation
  logical,save,public :: ante_ipq0 ! Before interpolation `Q0'
  logical,save,public :: post_ipq0 ! After interpolation `Q0'
  logical,save,public :: ante_rest ! Before restriction
  logical,save,public :: post_rest ! After restriction

  ! Residuals
  logical,save,public :: print_res ! Print maximum residual
  logical,save,public :: save_res ! Save maximum residual
  integer,save,public :: print_freq_res ! Steps between output
  integer,save,public :: save_freq_res ! Steps between output

  ! Errors if there is an exact solution
  logical,save,public :: print_err ! Shows time step tick error
  logical,save,public :: save_err ! Save time step tick error
  integer,save,public :: print_freq_err ! Steps between output
  integer,save,public :: save_freq_err ! Steps between output

  ! Numerical solution if there is no exact solution
  logical,save,public :: write_tic_sol ! Write time step tick sols using HDF
  logical,save,public :: write_toc_sol ! Write time step tock sols using HDF
  integer,save,public :: write_freq_sol ! Steps between output

  ! List of arguments
  namelist /arg_core_hook/ ante_pred,post_pred,ante_iter,post_iter,ante_sweep,post_sweep,ante_step,post_step,ante_prolo,post_prolo,ante_ipq0,post_ipq0,ante_rest,post_rest,ante_pred,post_pred,ante_iter,post_iter,ante_sweep,post_sweep,ante_step,post_step,ante_prolo,post_prolo,ante_ipq0,post_ipq0,ante_rest,post_rest,print_res,save_res,print_freq_res,save_freq_res,print_err,save_err,print_freq_err,save_freq_err,write_tic_sol,write_toc_sol,write_freq_sol

  ! ================================================================
  ! Data
  ! ================================================================

  ! Hook names
  character(len = max_len_name_),dimension(:),allocatable,save,public :: name_hook_
  ! Outer sol names, times, and values
  character(len = max_len_name_),dimension(:),allocatable,save,public :: name_outer_
  character(len = max_len_time_),dimension(:),allocatable,save,public :: time_outer_
  character(len = max_len_data_),dimension(:),allocatable,save,public :: data_outer_
  ! Inner sol names, times, and values
  character(len = max_len_name_),dimension(:),allocatable,save,public :: name_inner_
  character(len = max_len_time_),dimension(:),allocatable,save,public :: time_inner_
  character(len = max_len_data_),dimension(:),allocatable,save,public :: data_inner_
  ! Residual sol names, times, and values
  character(len = max_len_name_),dimension(:),allocatable,save,public :: name_res_
  character(len = max_len_time_),dimension(:),allocatable,save,public :: time_res_
  character(len = max_len_data_),dimension(:),allocatable,save,public :: data_res_
  ! Error sol names, times, and values
  character(len = max_len_name_),dimension(:),allocatable,save,public :: name_err_
  character(len = max_len_time_),dimension(:),allocatable,save,public :: time_err_
  character(len = max_len_data_),dimension(:),allocatable,save,public :: data_err_

  ! Base output unit for SNES sols
  integer,save,public :: base_io_unit_outer_
  ! Base output unit for KSP sols
  integer,save,public :: base_io_unit_inner_
  ! Base output unit for residual sols
  integer,save,public :: base_io_unit_res_
  ! Base output unit for error sols
  integer,save,public :: base_io_unit_err_

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_hook,destroy_hook,hook_print,hook_save,hook_log,hook_output
  private :: hook_read_arg_core_hook

contains

  subroutine create_hook(io_unit,name_file_input_pfasst,pfasst_pf)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_hook:create_hook(...)': Create output environment
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_INTEGER
    use pfasst,only : pf_pfasst_t,pf_add_hook,PF_PRE_PREDICTOR,PF_POST_PREDICTOR,PF_PRE_ITERATION,PF_POST_ITERATION,PF_PRE_SWEEP,PF_POST_SWEEP,PF_PRE_STEP,PF_POST_STEP,PF_PRE_INTERP_ALL,PF_POST_INTERP_ALL,PF_PRE_INTERP_Q0,PF_POST_INTERP_Q0,PF_PRE_RESTRICT_ALL,PF_POST_RESTRICT_ALL

    ! Project
    use mod_mpi,only : mpi_space_,mpi_time_,mpi_exit_gracefully
    use mod_pf,only : num_level
    use mod_ctx,only : num_step_space,num_var_dep,all_f_ctx_

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: io_unit
    character(len = *),intent(in) :: name_file_input_pfasst
    type(pf_pfasst_t),intent(inout) :: pfasst_pf

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: level,alloc_stat,dealloc_stat,mpi_stat
    integer,dimension(:),pointer :: all_mem_counts

    ! ================================================================
    ! Work
    ! ================================================================

    ! Hook names
    allocate(character(len = max_len_name_) :: name_hook_(num_hook_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `name_hook_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! List of hook names
    write(name_hook_(1),"(a)") "ante_block"
    write(name_hook_(2),"(a)") "post_block"
    write(name_hook_(3),"(a)") "ante_pred "
    write(name_hook_(4),"(a)") "post_pred "
    write(name_hook_(5),"(a)") "ante_iter "
    write(name_hook_(6),"(a)") "post_iter "
    write(name_hook_(7),"(a)") "ante_sweep"
    write(name_hook_(8),"(a)") "post_sweep"
    write(name_hook_(9),"(a)") "ante_step "
    write(name_hook_(10),"(a)") "post_step "
    write(name_hook_(11),"(a)") "ante_prolo"
    write(name_hook_(12),"(a)") "post_prolo"
    write(name_hook_(13),"(a)") "ante_ipq0 "
    write(name_hook_(14),"(a)") "post_ipq0 "
    write(name_hook_(15),"(a)") "ante_rest "
    write(name_hook_(16),"(a)") "post_rest "
    write(name_hook_(17),"(a)") "ante_conv "
    write(name_hook_(18),"(a)") "post_conv "

    ! Outer sol names
    allocate(character(len = max_len_name_) :: name_outer_(num_outer_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `name_outer_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! List of SNES sol names
    write(name_outer_(1),"(a)") "outer_iter"
    write(name_outer_(2),"(a)") "outer_res "

    ! Outer times
    allocate(character(len = max_len_time_) :: time_outer_(num_outer_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `time_outer_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Outer values
    allocate(character(len = max_len_data_) :: data_outer_(num_outer_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `data_outer_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Inner sol names
    allocate(character(len = max_len_name_) :: name_inner_(num_inner_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `name_inner_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! List of KSP sol names
    write(name_inner_(1),"(a)") "inner_iter"
    write(name_inner_(2),"(a)") "inner_res "

    ! Inner times
    allocate(character(len = max_len_time_) :: time_inner_(num_inner_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `time_inner_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Inner values
    allocate(character(len = max_len_data_) :: data_inner_(num_inner_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `data_inner_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Residual sol names
    allocate(character(len = max_len_name_) :: name_res_(num_res_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `name_res_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! List of residual sol names
    write(name_res_(1),"(a)") "res_max"
    write(name_res_(2),"(a)") "res_tic"
    write(name_res_(3),"(a)") "res_toc"

    ! Residual times
    allocate(character(len = max_len_time_) :: time_res_(num_res_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `time_res_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Residual values
    allocate(character(len = max_len_data_) :: data_res_(num_res_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `data_res_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Error sol names
    allocate(character(len = max_len_name_) :: name_err_(num_err_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `name_err_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! List of error sol names
    write(name_err_(1),"(a)") "err_tic"
    write(name_err_(2),"(a)") "err_toc"

    ! Error times
    allocate(character(len = max_len_time_) :: time_err_(num_err_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `time_err_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Error values
    allocate(character(len = max_len_data_) :: data_err_(num_err_),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `data_err_' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Base output unit for SNES sols
    base_io_unit_outer_ = 0
    ! Base output unit for KSP sols
    base_io_unit_inner_ = base_io_unit_outer_+num_outer_*num_hook_*mpi_time_%size
    ! Base output unit for residual sols
    base_io_unit_res_ = base_io_unit_inner_+num_inner_*num_hook_*mpi_time_%size
    ! Base output unit for error sols
    base_io_unit_err_ = base_io_unit_res_+num_res_*num_hook_*mpi_time_%size

    ! Parse input file
    call hook_read_arg_core_hook(io_unit,name_file_input_pfasst)

    ! Write data before the predictor
    if(ante_pred) then
      call pf_add_hook(pfasst_pf,-1,PF_PRE_PREDICTOR,hook_output)
    end if ! `ante_pred'

    ! Write data after the predictor
    if(post_pred) then
      call pf_add_hook(pfasst_pf,-1,PF_POST_PREDICTOR,hook_output)
    end if ! `post_pred'

    ! Write data before every iter.
    if(ante_iter) then
      call pf_add_hook(pfasst_pf,-1,PF_PRE_ITERATION,hook_output)
    end if ! `ante_iter'

    ! Write data after every iter.
    if(post_iter) then
      call pf_add_hook(pfasst_pf,-1,PF_POST_ITERATION,hook_output)
    end if ! `post_iter'

    ! Write data before every sweep
    if(ante_sweep) then
      call pf_add_hook(pfasst_pf,-1,PF_PRE_SWEEP,hook_output)
    end if ! `ante_sweep'

    ! Write data after every sweep
    if(post_sweep) then
      call pf_add_hook(pfasst_pf,-1,PF_POST_SWEEP,hook_output)
    end if ! `post_sweep'

    ! Write data before every step
    if(ante_step) then
      call pf_add_hook(pfasst_pf,-1,PF_PRE_STEP,hook_output)
    end if ! `ante_step'

    ! Write data after every step
    if(post_step) then
      call pf_add_hook(pfasst_pf,-1,PF_POST_STEP,hook_output)
    end if ! `post_step'

    ! Write data before every interpolation call
    if(ante_prolo) then
      call pf_add_hook(pfasst_pf,-1,PF_PRE_INTERP_ALL,hook_output)
    end if ! `ante_prolo'

    ! Write data after every interpolation call
    if(post_prolo) then
      call pf_add_hook(pfasst_pf,-1,PF_POST_INTERP_ALL,hook_output)
    end if ! `post_prolo'

    ! Write data before every interpolation `Q0'' call
    if(ante_ipq0) then
      call pf_add_hook(pfasst_pf,-1,PF_PRE_INTERP_Q0,hook_output)
    end if ! `ante_ipq0'

    ! Write data after every interpolation `Q0'' call
    if(post_ipq0) then
      call pf_add_hook(pfasst_pf,-1,PF_POST_INTERP_Q0,hook_output)
    end if ! `post_ipq0'

    ! Write data before every restriction call
    if(ante_rest) then
      call pf_add_hook(pfasst_pf,-1,PF_PRE_RESTRICT_ALL,hook_output)
    end if ! `ante_rest'

    ! Write data after every restriction call
    if(post_rest) then
      call pf_add_hook(pfasst_pf,-1,PF_POST_RESTRICT_ALL,hook_output)
    end if ! `post_rest'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_hook

  subroutine destroy_hook()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_hook:create_hook(...)': Create output environment
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : mpi_exit_gracefully
    use mod_pf,only : num_level
    use mod_ctx,only : num_step_space

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

    integer :: level,dealloc_stat

    ! ================================================================
    ! Work
    ! ================================================================

    ! Error values
    deallocate(data_err_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `data_err_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Error times
    deallocate(time_err_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `time_err_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Error sol names
    deallocate(name_err_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `name_err_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Residual values
    deallocate(data_res_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `data_res_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Residual times
    deallocate(time_res_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `time_res_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Residual sol names
    deallocate(name_res_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `name_res_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Inner values
    deallocate(data_inner_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `data_inner_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Inner times
    deallocate(time_inner_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `time_inner_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Inner sol names
    deallocate(name_inner_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `name_inner_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Outer values
    deallocate(data_outer_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `data_outer_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Outer times
    deallocate(time_outer_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `time_outer_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Outer sol names
    deallocate(name_outer_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `name_outer_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Hook names
    deallocate(name_hook_,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `name_hook_' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_hook

  subroutine hook_read_arg_core_hook(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_hook:hook_read_arg_core_hook(...)': Read input arguments
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHAR

    ! Project
    use mod_mpi,only : max_len_char_,num_dim_space,mpi_world_,mpi_time_,mpi_exit_gracefully
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

    ! Hooks to attach for output
    ante_pred  = .false. ! Before predictor
    post_pred  = .false. ! After predictor
    ante_iter  = .false. ! Before iteration
    post_iter  = .true. ! After iteration
    ante_sweep = .false. ! Before sweep
    post_sweep = .false. ! After sweep
    ante_step  = .false. ! Before step
    post_step  = .false. ! After step
    ante_prolo = .false. ! Before prolongation
    post_prolo = .false. ! After prolongation
    ante_ipq0 = .false. ! Before interpolation `Q0'
    post_ipq0 = .false. ! After interpolation `Q0'
    ante_rest  = .false. ! Before restriction
    post_rest  = .false. ! After restriction

    ! Residuals
    print_res = .true. ! Print maximum residual
    save_res = .false. ! Save maximum residual
    print_freq_res = 1 ! Steps between output
    save_freq_res = 1 ! Steps between output

    ! Errors if there is an exact solution
    print_err = .false. ! Shows time step tick error
    save_err = .false. ! Save time step tick error
    print_freq_err = 1 ! Steps between output
    save_freq_err = 1 ! Steps between output

    ! Numerical solution if there is no exact solution
    write_tic_sol = .false. ! Write time step tick sols using HDF
    write_toc_sol = .false. ! Write time step tock sols using HDF
    write_freq_sol = 1 ! Steps between output

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_hook)
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
          namelist_arg = "&arg_core_hook "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_core_hook,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(ante_pred,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(post_pred,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(ante_iter,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(post_iter,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(ante_sweep,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(post_sweep,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(ante_step,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(post_step,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(ante_prolo,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(post_prolo,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(ante_ipq0,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(post_ipq0,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(ante_rest,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(post_rest,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_res,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_res,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_freq_res,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_freq_res,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_err,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_err,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(print_freq_err,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(save_freq_err,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(write_tic_sol,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(write_toc_sol,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(write_freq_sol,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine hook_read_arg_core_hook

  subroutine hook_print(num_obj,name_obj,time_obj,data_obj,hook,step_time,level,sweep,iter)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_hook:hook_print(...)': Write concatenated output
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : mpi_space_,mpi_time_

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: num_obj,hook,step_time,level,sweep,iter
    character(len = max_len_name_),dimension(:),intent(in) :: name_obj
    character(len = max_len_time_),dimension(:),intent(in) :: time_obj
    character(len = max_len_data_),dimension(:),intent(in) :: data_obj

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: i_obj

    ! ================================================================
    ! Work
    ! ================================================================

    ! Print and save from space master
    if(mpi_space_%rank .eq. mpi_space_%size-1) then
      ! Print from time master
      if(mpi_time_%rank .eq. mpi_time_%size-1) then
        ! Output sol
        do i_obj = 1,num_obj
          ! Write output line
          write(6,"(a,a,a,a,a,a,a,1i0,a,1i0,a,a,a,1i0,a,1i0,a,1i0)") "# ",trim(adjustL(name_obj(i_obj)))," ",trim(adjustL(data_obj(i_obj)))," # hook ",trim(adjustL(name_hook_(hook)))," # rank_time ",mpi_time_%rank," # step_time ",step_time," # val_time ",trim(adjustL(time_obj(i_obj)))," # level ",level," # sweep ",sweep," # iter ",iter
          call flush(6)
        end do ! `i_obj'
      end if ! `mpi_time_%rank .eq. mpi_time_%size-1'
    end if ! `mpi_space_%rank .eq. mpi_space_%size-1'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine hook_print

  subroutine hook_save(num_obj,name_obj,time_obj,data_obj,base_io_unit_obj,save_freq_obj,hook,step_time,level,sweep,iter)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_hook:hook_save(...)': Write concatenated output
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : max_len_char_,mpi_space_,mpi_time_
    use mod_pf,only : num_step_time

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: num_obj,base_io_unit_obj,save_freq_obj,hook,step_time,level,sweep,iter
    character(len = max_len_name_),dimension(:),intent(in) :: name_obj
    character(len = max_len_time_),dimension(:),intent(in) :: time_obj
    character(len = max_len_data_),dimension(:),intent(in) :: data_obj

    ! ================================================================
    ! Locals
    ! ================================================================

    logical :: is_open,does_exist
    integer :: i_obj,io_unit,io_stat
    character(len = max_len_char_) :: name_file

    ! ================================================================
    ! Work
    ! ================================================================

    ! Print and save from space master
    if(mpi_space_%rank .eq. mpi_space_%size-1) then
      ! Output sol
      do i_obj = 1,num_obj
        ! Set output unit
        io_unit = base_io_unit_+base_io_unit_obj+(i_obj-1)*num_hook_*mpi_time_%size+(hook-1)*mpi_time_%size+mpi_time_%rank
        ! Check file by unit
        is_open = .false.
        inquire(unit = io_unit,opened = is_open)
        ! Open file by unit
        if(.not. is_open) then
          ! Define file name
          write(name_file,"(a,a,a,a,a,1i0,a)") "./",trim(adjustL(name_obj(i_obj))),".",trim(adjustL(name_hook_(hook))),".rank_time",mpi_time_%rank,".txt"
          ! Check file by name
          does_exist = .false.
          inquire(file = trim(adjustL(name_file)),exist = does_exist)
          ! Open file
          open(io_unit,file = trim(adjustL(name_file)),action = "write",access = "sequential",position = "append",status = "unknown",iostat = io_stat)
            ! Write header
            if(.not. does_exist) then
              write(io_unit,"(a,a,a)") "# ",trim(adjustL(name_obj(i_obj)))," # hook # rank_time # step_time # val_time # level # sweep # iter"
              call flush(io_unit)
            end if
        end if
        ! Write data
        write(io_unit,"(a,1x,a,1x,1i0,1x,1i0,1x,a,1x,1i0,1x,1i0,1x,1i0)") trim(adjustL(data_obj(i_obj))),trim(adjustL(name_hook_(hook))),mpi_time_%rank,step_time,trim(adjustL(time_obj(i_obj))),level,sweep,iter
        call flush(io_unit)
        if(step_time+1+save_freq_obj .gt. num_step_time) then
          close(io_unit)
        end if ! `state_pf_%step+1+save_freq_obj .gt. num_step_time'
      end do ! `i_obj'
    end if ! `mpi_space_%rank .eq. mpi_space_%size-1'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine hook_save

  subroutine hook_log(io_unit,data,state)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_hook:hook_log(...)': Write concatenated output
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

    integer,intent(in) :: io_unit
    character(len = *),intent(in) :: data,state

    ! ================================================================
    ! Locals
    ! ================================================================

    logical :: is_open

    ! ================================================================
    ! Work
    ! ================================================================

    ! Get file status
    inquire(unit = io_unit,opened = is_open)
    if(.not. is_open) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Output unit must be writeable.")
    end if ! `.not. is_open'

    ! Write data
    write(io_unit,"(a,a,a)") trim(state)," ",trim(data)

    ! Update file content
    flush(io_unit) 

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine hook_log

  subroutine hook_output(pfasst_pf,level_pf,state_pf,c_ctx)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_hook:hook_output(...)': Print routine for hooks
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Fortran
    use iso_c_binding,only : c_ptr

    ! External
    use pfasst,only : pf_pfasst_t,pf_level_t,pf_state_t

    ! Project
    use mod_pf,only : num_level,all_num_node_time,num_step_time
    use mod_mpi,only : max_len_char_,mpi_world_,mpi_space_,mpi_time_
    use mod_ctx,only : num_var_dep,all_f_ctx_,ctx_write
    use mod_prob,only : prob_get_err
    use mod_encap,only : encap_get_norm,encap_print,encap_pack

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    type(pf_pfasst_t),intent(inout) :: pfasst_pf
    type(pf_level_t),intent(inout) :: level_pf
    type(pf_state_t),intent(in) :: state_pf
    type(c_ptr),intent(in) :: c_ctx

    ! ================================================================
    ! Locals
    ! ================================================================

    logical :: is_open,does_exist
    integer :: i_res,i_err,io_unit,io_stat,var_dep,petsc_stat
    character(len = max_len_char_) :: val_time,name_file

    ! ================================================================
    ! Work
    ! ================================================================

    if((print_res .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,print_freq_res) .eq. 0))) .or. (save_res .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,save_freq_res) .eq. 0)))) then
      ! Time values
      write(time_res_(1),"(1es24.16e3)") (2*state_pf%t0+state_pf%dt)/2
      write(time_res_(2),"(1es24.16e3)") state_pf%t0
      write(time_res_(3),"(1es24.16e3)") state_pf%t0+state_pf%dt

      ! Field values
      write(data_res_(1),"(1es24.16e3)") level_pf%residual
      write(data_res_(2),"(1es24.16e3)") level_pf%residual_tic
      write(data_res_(3),"(1es24.16e3)") level_pf%residual_toc

      ! Print data
      if(print_res) then
        call hook_print(num_res_,name_res_,time_res_,data_res_,state_pf%hook,state_pf%step,level_pf%level-1,state_pf%sweep-1,state_pf%iter-1)
      end if ! ``print_res''

      ! Save data
      if(save_res) then
        call hook_save(num_res_,name_res_,time_res_,data_res_,base_io_unit_res_,save_freq_res,state_pf%hook,state_pf%step,level_pf%level-1,state_pf%sweep-1,state_pf%iter-1)
      end if ! ``save_res''
    end if ! `(print_res .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,print_freq_res) .eq. 0))) .or. (save_res .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,save_freq_res) .eq. 0)))'

    if((print_err .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,print_freq_err) .eq. 0))) .or. (save_err .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,save_freq_err) .eq. 0)))) then
      ! Time values
      write(time_err_(1),"(1es24.16e3)") state_pf%t0
      write(time_err_(2),"(1es24.16e3)") state_pf%t0+state_pf%dt

      ! Field values
      write(data_err_(1),"(1es24.16e3)") prob_get_err(state_pf%t0,level_pf%Q(1))
      write(data_err_(2),"(1es24.16e3)") prob_get_err(state_pf%t0+state_pf%dt,level_pf%Q(all_num_node_time(level_pf%level)))

      ! Print data
      if(print_err) then
        call hook_print(num_err_,name_err_,time_err_,data_err_,state_pf%hook,state_pf%step,level_pf%level-1,state_pf%sweep-1,state_pf%iter-1)
      end if ! ``print_err''

      ! Save data
      if(save_err) then
        call hook_save(num_err_,name_err_,time_err_,data_err_,base_io_unit_err_,save_freq_err,state_pf%hook,state_pf%step,level_pf%level-1,state_pf%sweep-1,state_pf%iter-1)
      end if ! ``save_err''
    end if ! `(print_err .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,print_freq_err) .eq. 0))) .or. (save_err .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,save_freq_err) .eq. 0)))'

    ! Write time step tick sols using HDF
    if(write_tic_sol .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,write_freq_sol) .eq. 0))) then
      ! Define file name
      write(name_file,"(a,a,a,1i0,a,1i0,a,1i0,a,1i0,a,1i0,a)") "./sol_tic.",trim(adjustL(name_hook_(state_pf%hook))),".rank_time",mpi_time_%rank,".step_time",state_pf%step,".level",level_pf%level-1,".sweep",state_pf%sweep-1,".iter",state_pf%iter-1,".hdf"

      ! Pack data
      call encap_pack(all_f_ctx_(level_pf%level)%f_ctx%hdf,level_pf%Q(1))

      ! Write file
      call ctx_write(level_pf%level,name_file)
    end if ! `write_tic_sol .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,write_freq_sol) .eq. 0))'

    ! Write time step tock sols using HDF
    if(write_toc_sol .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,write_freq_sol) .eq. 0))) then
      ! Define file name
      write(name_file,"(a,a,a,1i0,a,1i0,a,1i0,a,1i0,a,1i0,a)") "./sol_toc.",trim(adjustL(name_hook_(state_pf%hook))),".rank_time",mpi_time_%rank,".step_time",state_pf%step,".level",level_pf%level-1,".sweep",state_pf%sweep-1,".iter",state_pf%iter-1,".hdf"

      ! Pack data
      call encap_pack(all_f_ctx_(level_pf%level)%f_ctx%hdf,level_pf%Q(all_num_node_time(level_pf%level)))

      ! Write file
      call ctx_write(level_pf%level,name_file)
    end if ! `write_toc_sol .and. ((state_pf%step .eq. 0) .or. (mod(state_pf%step+1,write_freq_sol) .eq. 0))'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine hook_output

end module mod_hook
