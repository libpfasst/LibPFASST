! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

module mod_mpi

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_mpi': Define parallel programming environment
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

  ! Buffer size for strings
  integer,parameter,public :: max_len_char_ = 16384
  ! Number of spatial dimensions
  integer,dimension(2),parameter,private :: avail_num_dim_space_ = (/1,2/)

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! Number of spatial dimensions
  integer,save,public :: num_dim_space
  
  ! List of arguments
  namelist /arg_core_mpi/ num_dim_space

  ! Spatial distribution of MPI ranks
  integer,dimension(:),allocatable,save,public :: num_rank_space
  ! Number of ranks for time
  integer,save,public :: num_rank_time
  
  ! List of arguments
  namelist /arg_mantle_mpi/ num_rank_space,num_rank_time

  ! ================================================================
  ! Data
  ! ================================================================

  ! MPI type
  type,public :: type_mpi
    integer :: comm,size,rank,info
  end type type_mpi

  ! World, space, and time
  type(type_mpi),save,public :: mpi_world_,mpi_space_,mpi_time_

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_mpi,destroy_mpi,mpi_sanity_check,mpi_exit_gracefully
  private :: mpi_read_arg_core_mpi,mpi_read_arg_mantle_mpi

contains

  subroutine create_mpi(mpi_world_comm,io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_mpi:create_mpi(...)': Define parallel programming world
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_Comm_size,MPI_Comm_rank,MPI_Comm_split,MPI_INFO_NULL

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    integer,intent(in) :: mpi_world_comm,io_unit
    character(len = *),intent(in) :: name_file_input_pfasst

    ! ================================================================
    ! Locals
    ! ================================================================

    integer :: mpi_stat,alloc_stat

    ! ================================================================
    ! Work
    ! ================================================================

    ! World communicator definition
    call MPI_Comm_dup(mpi_world_comm,mpi_world_%comm,mpi_stat)

    ! World size definition
    call MPI_Comm_size(mpi_world_%comm,mpi_world_%size,mpi_stat)

    ! World rank definition
    call MPI_Comm_rank(mpi_world_%comm,mpi_world_%rank,mpi_stat)

    ! Get core arguments from input
    call mpi_read_arg_core_mpi(io_unit,name_file_input_pfasst)

    ! Check for number of spatial dimensions
    call mpi_sanity_check(avail_num_dim_space_,num_dim_space)

    ! Get MPI rank related arguments from input
    call mpi_read_arg_mantle_mpi(io_unit,name_file_input_pfasst)

    ! Size of space communicator
    mpi_space_%size = product(num_rank_space)

    ! Size of time communicator
    mpi_time_%size = num_rank_time

    ! Require remainder-free integer division
    if(mpi_space_%size*mpi_time_%size .ne. mpi_world_%size) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Split of `mpi_world_%comm' failed. Change either the size of `MPI_COMM_WORLD', `num_rank_space', or `num_rank_time'. Keep in mind the dependence of `num_rank_space' on `num_dim_space'.")
    end if ! `mpi_space_%size*mpi_time_%size .ne. mpi_world_%size'

    ! Time rank
    mpi_time_%rank = mpi_world_%rank/mpi_space_%size

    ! Space communicator
    call MPI_Comm_split(mpi_world_%comm,mpi_time_%rank,mpi_world_%rank,mpi_space_%comm,mpi_stat)

    ! Space rank
    call MPI_Comm_rank(mpi_space_%comm,mpi_space_%rank,mpi_stat)

    ! Time communicator
    call MPI_Comm_split(mpi_world_%comm,mpi_space_%rank,mpi_world_%rank,mpi_time_%comm,mpi_stat)

    ! World info
    mpi_world_%info = MPI_INFO_NULL

    ! Space info
    mpi_space_%info = MPI_INFO_NULL

    ! Time info
    mpi_time_%info = MPI_INFO_NULL

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_mpi

  subroutine destroy_mpi()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_mpi:destroy_mpi(...)': Free occupied memory
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================ 
    ! Modules
    ! ================================================================ 

    ! External
    use mpi,only : mpi_comm_free

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

    integer :: dealloc_stat,mpi_stat

    ! ================================================================
    ! Work
    ! ================================================================

    ! Number of ranks for space
    deallocate(num_rank_space,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `num_rank_space' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! Destroy MPI communicators
    call mpi_comm_free(mpi_time_%comm,mpi_stat)
    call mpi_comm_free(mpi_space_%comm,mpi_stat)
    call mpi_comm_free(mpi_world_%comm,mpi_stat)

    ! Print job completion statement to screen
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      write(6,"(a,a,a,1i0,a)") "`",trim(__FILE__),":",__LINE__,"': Simulation ends"
      call flush(6)
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_mpi

  subroutine mpi_read_arg_core_mpi(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_mpi:mpi_read_arg_core_mpi(...)': Read input arguments
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_INTEGER

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

    ! Number of spatial dimensions
    num_dim_space = 1

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_mpi)
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
          namelist_arg = "&arg_core_mpi "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_core_mpi,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(num_dim_space,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine mpi_read_arg_core_mpi

  subroutine mpi_read_arg_mantle_mpi(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_mpi:mpi_read_arg_mantle_mpi(...)': Read input arguments
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_INTEGER

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

    ! Number of ranks for space
    allocate(num_rank_space(num_dim_space),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `num_rank_space(...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Number of ranks for space
    num_rank_space = 1

    ! Number of ranks for time
    num_rank_time = 1

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_mantle_mpi)
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
          namelist_arg = "&arg_mantle_mpi "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_mantle_mpi,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(num_rank_space,size(num_rank_space),MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)
    call MPI_bCast(num_rank_time,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine mpi_read_arg_mantle_mpi

  subroutine mpi_sanity_check(avail_ints,curr_int)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_mpi:mpi_sanity_check(...)': Make sure `curr_int' is contained in `avail_ints'
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================ 
    ! Modules
    ! ================================================================ 

    ! ================================================================ 
    ! Safety
    ! ================================================================ 

    implicit none

    ! ================================================================ 
    ! Dummies
    ! ================================================================ 

    integer,dimension(:),intent(in) :: avail_ints
    integer,intent(in) :: curr_int

    ! ================================================================ 
    ! Locals
    ! ================================================================ 

    logical,dimension(size(avail_ints)) :: ints_mask

    ! ================================================================ 
    ! Work
    ! ================================================================ 

    where(curr_int .eq. avail_ints)
      ints_mask = .true.
    elsewhere
      ints_mask = .false.
    end where ! `curr_int .eq. avail_ints'

    if(.not. any(ints_mask)) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Value for integer not allowed.")
    end if ! `any(ints_mask) .eqv. .false.'

    ! ================================================================ 
    ! Exit
    ! ================================================================ 

    return

  end subroutine mpi_sanity_check

  subroutine mpi_exit_gracefully(file_name,line_number,msg)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_mpi:mpi_exit_gracefully(...)': Abort simulation at error
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    character(len = *),intent(in) :: file_name,msg
    integer,intent(in) :: line_number

    ! ================================================================
    ! Locals
    ! ================================================================

    ! ================================================================
    ! Work
    ! ================================================================

    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      write(6,"(a,a,a,1i0,a,a)") "`",trim(file_name),":",line_number,"': ",trim(msg)
      call flush(6)
      stop
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine mpi_exit_gracefully

end module mod_mpi
