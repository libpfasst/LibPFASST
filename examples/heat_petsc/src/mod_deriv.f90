! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

module mod_deriv

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `mod_deriv': Define containers with derivative specifics
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

  ! Possible derivative accuaracies
  integer,dimension(4),parameter,private :: avail_deriv_accs_ = (/2,4,6,8/)

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! Current derivative accuracy
  integer,save,public :: deriv_acc
  
  ! List of arguments
  namelist /arg_core_deriv/ deriv_acc

  ! ================================================================
  ! Data
  ! ================================================================

  ! Width of derivative stencil
  integer,save,public :: deriv_wid
  ! 1st-order derivative coefficients per PFASST level and spatial direction
  real(pfdp),dimension(:,:,:),allocatable,save,public :: ord1_coeffs
  ! 2nd-order derivative coefficients per PFASST level and spatial direction
  real(pfdp),dimension(:,:,:),allocatable,save,public :: ord2_coeffs

  ! ================================================================
  ! Interfaces
  ! ================================================================

  public :: create_deriv,destroy_deriv,deriv_set_coeffs
  private :: deriv_read_arg_core_deriv

contains

  subroutine create_deriv(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_deriv:create_deriv(...)': Set up objects for derivatives
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! Project
    use mod_pf,only : num_level
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

    ! Get arguments from input
    call deriv_read_arg_core_deriv(io_unit,name_file_input_pfasst)

    ! Check for order of derivative accuracy
    call mpi_sanity_check(avail_deriv_accs_,deriv_acc)

    ! Derive width of derivative stencil
    deriv_wid = deriv_acc/2

    ! Reserve memory for all 1st-order derivative coefficients
    allocate(ord1_coeffs(deriv_acc+1,num_dim_space,num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `ord1_coeffs(...,...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! Reserve memory for all 2nd-order derivative coefficients
    allocate(ord2_coeffs(deriv_acc+1,num_dim_space,num_level),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `ord2_coeffs(...,...,...)' failed.")
    end if ! `alloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine create_deriv

  subroutine destroy_deriv()

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_deriv:destroy_deriv(...)': Set free occupied memory
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

    ! 2nd-order derivative coefficients
    deallocate(ord2_coeffs,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `ord2_coeffs' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! 1st-order derivative coefficients
    deallocate(ord1_coeffs,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `ord1_coeffs' failed.")
    end if ! `dealloc_stat .ne. 0'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine destroy_deriv

  subroutine deriv_read_arg_core_deriv(io_unit,name_file_input_pfasst)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_deriv:deriv_read_arg_core_deriv(...)': Read arguments from input
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

    ! Default order of accuracy for derivative stencil
    deriv_acc = 2

    ! Read on first rank
    if(mpi_world_%rank .eq. mpi_world_%size-1) then
      ! Read namelist from input file
      open(unit = io_unit,file = name_file_input_pfasst,action = "read",status = "old")
        read(unit = io_unit,nml = arg_core_deriv)
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
          namelist_arg = "&arg_core_deriv "//trim(cmd_line_arg)//" /"
          read(namelist_arg,nml = arg_core_deriv,iostat = io_stat,iomsg = io_msg)
        end if ! `num_arg_read .gt. 0'
        num_arg_read = num_arg_read+1
      end do ! `num_arg_read'
    end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

    ! Spread the word
    call MPI_bCast(deriv_acc,1,MPI_INTEGER,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine deriv_read_arg_core_deriv

  subroutine deriv_set_coeffs(size_step_space)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! `mod_deriv:deriv_set_coeffs(...)': Initialize derivative coefficients by value
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    ! ================================================================
    ! Modules
    ! ================================================================

    ! External
    use pfasst,only : pfdp

    ! Project
    use mod_pf,only : num_level
    use mod_mpi,only : num_dim_space

    ! ================================================================
    ! Safety
    ! ================================================================

    implicit none

    ! ================================================================
    ! Dummies
    ! ================================================================

    real(pfdp),dimension(num_dim_space,num_level),intent(in) :: size_step_space

    ! ================================================================
    ! Locals
    ! ================================================================
    
    integer :: i_coeff

    ! ================================================================
    ! Work
    ! ================================================================

    ! 1st-order derivative coefficients (1st half)
    select case(deriv_acc)
    case(2)
      ord1_coeffs(1,:,:) = -1d0/2d0
      ord1_coeffs(2,:,:) = 0d0
    case(4)
      ord1_coeffs(1,:,:) = 1d0/1.2d+1
      ord1_coeffs(2,:,:) = -2d0/3d0
      ord1_coeffs(3,:,:) = 0d0
    case(6)
      ord1_coeffs(1,:,:) = -1d0/6d+1
      ord1_coeffs(2,:,:) = 3d0/2d+1
      ord1_coeffs(3,:,:) = -3d0/4d0
      ord1_coeffs(4,:,:) = 0d0
    case(8)
      ord1_coeffs(1,:,:) = 1d0/2.8d+2
      ord1_coeffs(2,:,:) = -4d0/1.05d+2
      ord1_coeffs(3,:,:) = 1d0/5d0
      ord1_coeffs(4,:,:) = -4d0/5d0
      ord1_coeffs(5,:,:) = 0d0
    end select ! `deriv_acc'

    ! 1st-order derivative coefficients (2nd half)
    ord1_coeffs(deriv_wid+2:deriv_acc+1,:,:) = -ord1_coeffs(deriv_wid:1:-1,:,:)

    ! 1st-order derivative coefficients (divide by step size)
    do i_coeff = 1,deriv_acc+1
      ord1_coeffs(i_coeff,:,:) = ord1_coeffs(i_coeff,:,:)/size_step_space
    end do ! `i_coeff'

    ! 2nd-order derivative coefficients (1st half)
    select case(deriv_acc)
    case(2)
      ord2_coeffs(1,:,:) = 1d0
      ord2_coeffs(2,:,:) = -2d0
    case(4)
      ord2_coeffs(1,:,:) = -1d0/1.2d+1
      ord2_coeffs(2,:,:) = 4d0/3d0
      ord2_coeffs(3,:,:) = -5d0/2d0
    case(6)
      ord2_coeffs(1,:,:) = 1d0/9d+1
      ord2_coeffs(2,:,:) = -3d0/2d+1
      ord2_coeffs(3,:,:) = 3d0/2d0
      ord2_coeffs(4,:,:) = -4.9d+1/1.8d+1
    case(8)
      ord2_coeffs(1,:,:) = -1d0/5.6d+2
      ord2_coeffs(2,:,:) = 8d0/3.15d+2
      ord2_coeffs(3,:,:) = -1d0/5d0
      ord2_coeffs(4,:,:) = 8d0/5d0
      ord2_coeffs(5,:,:) = -2.05d+2/7.2d+1
    end select ! `deriv_acc'

    ! 2nd-order derivative coefficients (2nd half)
    ord2_coeffs(deriv_wid+2:deriv_acc+1,:,:) = ord2_coeffs(deriv_wid:1:-1,:,:)

    ! 2nd-order derivative coefficients (divide by step size)
    do i_coeff = 1,deriv_acc+1
      ord2_coeffs(i_coeff,:,:) = ord2_coeffs(i_coeff,:,:)/size_step_space**2d0
    end do ! `i_coeff'

    ! ================================================================
    ! Exit
    ! ================================================================

    return

  end subroutine deriv_set_coeffs

end module mod_deriv
