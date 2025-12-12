module pf_space_time_comm
  use pf_mod_mpi
  use pfasst
  implicit none
contains

  subroutine create_space_time_communicators(nspace, ntime, space_comm, time_comm, space_color, time_color, group_space)
    integer, intent(out) :: space_comm, time_comm
    integer, intent(out) :: space_color, time_color
    integer, intent(in) :: nspace, ntime
    logical, intent(in), optional :: group_space
    ! local variables
    integer :: gSize, gRank, err
    logical :: group_space_loc

    !> check size
    call mpi_comm_size(MPI_COMM_WORLD, gSize, err)
    call mpi_comm_rank(MPI_COMM_WORLD, gRank,  err)

    if (gSize .ne. (nspace * ntime)) then
       print '(a)', 'ERROR: create_simple_communicators: processor number mismatch.'
       print '(a,i4,a,i4)', '       Expecting ', &
            nspace * ntime, ' MPI processors but received ', &
            gSize
       stop
    end if


    !> set-up communicators
    group_space_loc = .TRUE.
    if (present(group_space)) then
        group_space_loc = group_space
    end if
     
    if (ntime == 1) then            ! full space parallel
       if (gRank .eq. 0) print *, "Full-space parallel"   
       time_color = gRank
       space_color = 0
       space_comm = MPI_COMM_WORLD
       time_comm  = MPI_COMM_SELF
    else if (nspace == 1) then      ! full time parallel
       if (gRank .eq. 0) print *, "Full-time parallel"   
       time_color = 0
       space_color = gRank
       space_comm = MPI_COMM_SELF
       time_comm  = MPI_COMM_WORLD
    else                            ! space-time parallel
        if (group_space_loc) then
            if (gRank .eq. 0) print *, "Space-time parallel with space grouping"          
            ! space procs grouped
            space_color = (gRank - mod(gRank,nspace)) / nspace
            call mpi_comm_split(MPI_COMM_WORLD, space_color, gRank, space_comm, err)
            call mpi_barrier(MPI_COMM_WORLD, err)
            
            time_color = mod(gRank, nspace)
            call mpi_comm_split(MPI_COMM_WORLD, time_color, gRank, time_comm, err)
            call mpi_barrier(MPI_COMM_WORLD, err)
        else                            
            if (gRank .eq. 0) print *, "Space-time parallel with time grouping"          
            ! time procs grouped
            space_color = mod(gRank, ntime) 
            call mpi_comm_split(MPI_COMM_WORLD, space_color, gRank, space_comm, err)
            call mpi_barrier(MPI_COMM_WORLD, err)

            time_color = (gRank - mod(gRank,ntime) / nspace)
            call mpi_comm_split(MPI_COMM_WORLD, time_color, gRank, time_comm, err)
            call mpi_barrier(MPI_COMM_WORLD, err)
        end if
    end if

  end subroutine create_space_time_communicators

end module pf_space_time_comm
