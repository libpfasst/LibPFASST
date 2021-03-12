module pf_space_comm
  use pf_mod_mpi
  use pfasst
  implicit none
contains

   subroutine create_simple_communicators(nspace, ntime, space_comm, time_comm, space_color, time_color, space_dim)
    integer, intent(out) :: space_comm, time_comm
    integer, intent(out) :: space_color, time_color
    integer, intent(in) :: nspace, ntime
    integer, intent(in) :: space_dim

    integer :: nproc, rank, error
    real(pfdp) :: nspace_real

    ! check size
    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

    if (nproc .ne. (nspace * ntime)) then
       print '(a)', 'ERROR: create_simple_communicators: processor number mismatch.'
       print '(a,i4,a,i4)', '       Expecting ', &
            nspace * ntime, ' MPI processors but received ', &
            nproc
       stop
    end if

    ! for now, nspace must be a perfect square
    if (space_dim .eq. 2) then
       nspace_real = sqrt(real(nspace))
       if (nspace_real .ne. nint(nspace_real)) then
          print'(a)', 'ERROR: create_simple_communicators: nspace must be perfect square.'
          stop
       end if
    end if

    if (ntime == 1) then
       time_color = rank
       space_color = 0
       space_comm = MPI_COMM_WORLD
       time_comm  = MPI_COMM_SELF
    else if (nspace == 1) then
       time_color = 0
       space_color = rank
       space_comm = MPI_COMM_SELF
       time_comm  = MPI_COMM_WORLD
    else
       ! split by color
       !space_color = mod(rank, nspace)
       !call mpi_comm_split(MPI_COMM_WORLD, space_color, rank, space_comm, error)

       !time_color = rank / nspace
       !call mpi_comm_split(MPI_COMM_WORLD, time_color, rank, time_comm, error)

       space_color = rank / nspace
       call mpi_comm_split(MPI_COMM_WORLD, space_color, rank, space_comm, error)

       time_color = mod(rank, nspace)
       call mpi_comm_split(MPI_COMM_WORLD, time_color, rank, time_comm, error)
    end if

  end subroutine create_simple_communicators

end module pf_space_comm
