! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

module mod_pf

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_pf': PFASST input arguments
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! External
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

  ! Number of PFASST levels
  integer,save,public :: num_level
  ! Number of PFASST iterations
  integer,save,public :: num_iter
  ! Legendre (`1') or Radau (`2') Gauss-Lobatto nodes in time (use `256+qType' if `all_num_node_time' [see below] requires proper inter-level interpolation and restriction)
  integer,save,public :: qType
  ! Window mode for convergence check
  integer,save,public :: window
  ! Absolute tolerance for iterations
  real(pfdp),save,public :: abs_res_tol
  ! Relative tolerance for iterations
  real(pfdp),save,public :: rel_res_tol
  ! Pipelined prediction
  logical,save,public :: pipeline_g
  ! Predictor activation
  logical,save,public :: pfasst_pred
  ! Calculate residual when entering a hook
  logical,save,public :: calc_residual
  ! Performance output 
  logical,save,public :: echo_timing
  ! Tau correction
  integer,save,public :: taui0
  ! Output directory
  character(len = 1),save,public :: out_dir

  ! List of arguments
  namelist /arg_core_pf/ num_level,num_iter,qType,window,abs_res_tol,rel_res_tol,pipeline_g,pfasst_pred,calc_residual,echo_timing,taui0,out_dir

  ! Number of time steps
  integer,save,public :: num_step_time
  ! Time step size (a nonpositive value is overwritten)
  real(pfdp),save,public :: size_step_time
  ! Time domain size (a positive `size_step_time' leads to this argument being inferred)
  real(pfdp),save,public :: size_dom_time
  ! Time domain origin
  real(pfdp),save,public :: orig_dom_time

  ! Number of `qType' nodes
  integer,dimension(:),allocatable,save,public :: all_num_node_time
  ! Number of predictor sweeps
  integer,dimension(:),allocatable,save,public :: all_num_sweep_pred
  ! Number of main sweeps
  integer,dimension(:),allocatable,save,public :: all_num_sweep_corr

  ! List of arguments
  namelist /arg_mantle_pf/ num_step_time,size_step_time,size_dom_time,orig_dom_time,all_num_node_time,all_num_sweep_pred,all_num_sweep_corr

  ! ================================================================
  ! Data
  ! ================================================================

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_pf,destroy_pf
  private :: pf_read_arg_core_pf,pf_read_arg_mantle_pf

contains

  subroutine create_pf(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_pf:create_pf(...)`: Allocate memory and define objects for PFASST
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHAR

    ! Project
    use mod_mpi,only : max_len_char_,num_rank_time,mpi_world_,mpi_time_,mpi_exit_gracefully

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

    integer :: num_arg_read,io_stat,mpi_stat,alloc_stat,i_level
    character(len = max_len_char_) :: cmd_line_arg,namelist_arg,io_msg

    ! ================================================================
    ! Work
    ! ================================================================

    ! Get input arguments for core of LIBPFASST
    call pf_read_arg_core_pf(io_unit,name_file_input_pfasst)

    ! Sanity check for PFASST
    if((mpi_time_%size .gt. 1) .and. (num_level .eq. 1)) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"There is nothing to parallelize in time yet we have `mpi_time_%size .gt. 1'.")
    end if ! `(mpi_time_%size .gt. 1) .and. (num_level .eq. 1)'

    ! Get input arguments for time domain
    call pf_read_arg_mantle_pf(io_unit,name_file_input_pfasst)

    ! The number of time PEs must divide the number of time steps
    if(mod(num_step_time,num_rank_time) .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"We need `mod(num_step_time,num_rank_time) .eq. 0'.")
    end if ! `mod(num_step_time,num_rank_time) .ne. 0'

    ! Time step size (a nonpositive value is overwritten)
    if(size_step_time .le. 0) then
      size_step_time = size_dom_time/num_step_time
    end if ! `size_step_time .le. 0'

    ! Time domain size (time steps size has priority)
    if(size_step_time .gt. 0) then
      size_dom_time = num_step_time*size_step_time
    end if ! `size_step_time .gt. 0'

    ! Sanity of temporal nodes
    do i_level = 1,num_level
      if((qType .eq. 1) .or. (qType .eq. 2)) then
        if((all_num_node_time(i_level) .ne. 2) .and. (all_num_node_time(i_level) .ne. 3) .and. (all_num_node_time(i_level) .ne. 5) .and. (all_num_node_time(i_level) .ne. 9)) then
          call mpi_exit_gracefully(__FILE__,__LINE__,"Time node number must be `2', `3', `5', or `9'.")
        end if ! `(all_num_node_time(i_level) .ne. 2) .and. (all_num_node_time(i_level) .ne. 3) .and. (all_num_node_time(i_level) .ne. 5) .and. (all_num_node_time(i_level) .ne. 9)'
      end if ! `(qType .eq. 1) .or. (qType .eq. 2)'
    end do ! `i_level'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_pf

  subroutine destroy_pf()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_pf:destroy_pf(...)': Set free occupied memory
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

    ! Number of sweeps for PFASST
    deallocate(all_num_sweep_corr,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_num_sweep_corr' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Number of sweeps for PFASST
    deallocate(all_num_sweep_pred,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_num_sweep_pred' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Specifics for time substeps
    deallocate(all_num_node_time,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `all_num_node_time' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_pf

  subroutine pf_read_arg_core_pf(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_pf:pf_read_arg_core_pf(...)': Read PFASST input arguments
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHAR

    ! Project
    use mod_mpi,only : max_len_char_,mpi_world_,mpi_time_,mpi_exit_gracefully

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

    integer :: num_arg_read,io_stat,mpi_stat
    character(len = max_len_char_) :: cmd_line_arg,namelist_arg,io_msg

    ! ================================================================
    ! Work
    ! ================================================================

    ! Number of PFASST levels
    num_level = 1
    ! Number of PFASST iterations
    num_iter = 8
    ! Legendre (`1') or Radau (`2') Gauss-Lobatto nodes in time (use `256+qType' if `all_num_node_time' [see below] requires proper inter-level interpolation and restriction)
    qType = 1
    ! Window mode for convergence check
    window = 1
    ! Absolute tolerance for iterations
    abs_res_tol = 0d0
    ! Relative tolerance for iterations
    rel_res_tol = 0d0
    ! Pipelined prediction
    pipeline_g = .true.
    ! Predictor activation
    pfasst_pred = .true.
    ! Calculate residual when entering a hook
    calc_residual = .true.
    ! Performance output 
    echo_timing = .false.
    ! Tau correction
    taui0 = -999999
    ! Output directory
    out_dir = ""

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_pf)
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
          namelist_arg = "&arg_core_pf "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_core_pf,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(num_level,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(num_iter,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(qType,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(window,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(abs_res_tol,1,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(rel_res_tol,1,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(pipeline_g,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(pfasst_pred,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(calc_residual,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(echo_timing,1,MPI_LOGICAL,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(taui0,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(out_dir,1,MPI_CHAR,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine pf_read_arg_core_pf

  subroutine pf_read_arg_mantle_pf(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_pf:pf_read_arg_mantle_pf(...)': Read PFASST input arguments
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_INTEGER,MPI_DOUBLE_PRECISION

    ! Project
    use mod_mpi,only : max_len_char_,mpi_world_,mpi_exit_gracefully

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

    ! Number of time steps
    num_step_time = 8
    ! Time step size (a nonpositive value is overwritten)
    size_step_time = -1.25d-1
    ! Time domain size (a positive `size_step_time' leads to this argument being inferred)
    size_dom_time = 1d0
    ! Time domain origin
    orig_dom_time = 0d0

    ! Specifics for time substeps
    allocate(all_num_node_time(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_num_node_time(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Number of `qType' nodes
    all_num_node_time = 3

    ! Number of sweeps for PFASST
    allocate(all_num_sweep_pred(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_num_sweep_pred(...)' failed.")
    end if ! `alloc_stat .ne. 0'
    
    ! Number of predictor sweeps
    all_num_sweep_pred = 1

    ! Number of sweeps for PFASST
    allocate(all_num_sweep_corr(num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `all_num_sweep_corr(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Number of main sweeps
    all_num_sweep_corr = 1

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_mantle_pf)
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
          namelist_arg = "&arg_mantle_pf "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_mantle_pf,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(num_step_time,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(size_step_time,1,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(size_dom_time,1,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(orig_dom_time,1,MPI_DOUBLE_PRECISION,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(all_num_node_time,size(all_num_node_time),MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(all_num_sweep_pred,size(all_num_sweep_pred),MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(all_num_sweep_corr,size(all_num_sweep_corr),MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine pf_read_arg_mantle_pf

end module mod_pf
