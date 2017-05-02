! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

module mod_interp

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_interp': Set up objects for PFASST level transfer routines
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

  ! Possible interpolation accuracies
  integer,dimension(4),parameter,private :: avail_interp_accs_ = (/2,4,6,8/)

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! Current interpolation accuracy
  integer,save,public :: interp_acc
  
  ! List of arguments
  namelist /arg_core_interp/ interp_acc

  ! ================================================================
  ! Data
  ! ================================================================

  ! Stencil width required for inerpolation
  integer,save,public :: interp_wid
  ! Interpolation weights for all spatial dimensions (*independent* of PFASST level)
  real(pfdp),dimension(:,:),allocatable,save,public :: weights

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_interp,destroy_interp,interp_set_weights
  private :: interp_read_arg_core_interp

contains

  subroutine create_interp(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_interp:create_interp(...)': Define interpolation weights
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : num_dim_space,mpi_sanity_check,mpi_exit_gracefully

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

    integer :: alloc_stat

    ! ================================================================
    ! Work
    ! ================================================================

    ! Get input arguments
    call interp_read_arg_core_interp(io_unit,name_file_input_pfasst)

    ! Check for order of interpative accuracy
    call mpi_sanity_check(avail_interp_accs_,interp_acc)

    ! Derive value for minimal width of inerpolation stencil
    interp_wid = interp_acc/2

    ! Set aside memory for interpolation weights in all spatial directions
    allocate(weights(interp_acc+1,num_dim_space),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `weights(...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_interp

  subroutine destroy_interp()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_interp:destroy_interp(...)': Set free occupied memory
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

    ! Free memory for interpolation weights in all spatial directions
    deallocate(weights,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `weights' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_interp

  subroutine interp_read_arg_core_interp(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_interp:interp_read_arg_core_interp(...)': Parse input file for interpolation arguments
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use mpi,only : MPI_INTEGER

    ! Project
    use mod_mpi,only : max_len_char_,mpi_world_

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

    ! Default order of accuracy for interpolation stencil
    interp_acc = 2

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_interp)
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
          namelist_arg = "&arg_core_interp "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_core_interp,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(interp_acc,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine interp_read_arg_core_interp

  subroutine interp_set_weights()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_interp:interp_set_weights(...)': Assign values to interpolation weights
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_mpi,only : num_dim_space

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

    ! ================================================================
    ! Work
    ! ================================================================

    ! Define interpolation weights (1st half)
    select case(interp_acc)
    case(2)
      weights(1,:) = 1d0

      weights = weights/2d0
    case(4)                                                        
      weights(1,:) = -1d0                                         
      weights(2,:) = 9d0                                          
                                          
      weights = weights/1.6d+1
    case(6)                                                        
      weights(1,:) = 3d0                                          
      weights(2,:) = -2.5d+1                                      
      weights(3,:) = 1.5d+2                                       
                                          
      weights = weights/2.56d+2
    case(8)                                                        
      weights(1,:) = -5d0                                         
      weights(2,:) = 4.9d+1                                       
      weights(3,:) = -2.45d+2                                     
      weights(4,:) = 1.225d+3                                     
                                          
      weights = weights/2.048d+3
    end select ! `interp_acc'

    ! Define interpolation weights (2nd half)
    weights(interp_wid+1:interp_acc,:) = weights(interp_wid:1:-1,:)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine interp_set_weights

end module mod_interp
