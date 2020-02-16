!!  MPI communicator routines
!
! This file is part of LIBPFASST.
!
!>  Module to hold use mpi statement
module pf_mod_mpi
  use  mpi
end module pf_mod_mpi

!> Module to implement communication routines in  MPI.
module pf_mod_comm_mpi
  use pf_mod_dtype

  use pf_mod_mpi, only: MPI_REAL16, MPI_REAL8
  implicit none
  !  For normal double precision
  integer, parameter :: myMPI_Datatype=MPI_REAL8  
  !  For  quadruple precision  (see top of pf_dtype.f90)
  ! integer, parameter :: myMPI_Datatype=MPI_REAL16  
contains

  !> Subroutine to create an MPI based PFASST communicator using the MPI communicator *mpi_comm*.
  subroutine pf_mpi_create(pf_comm, mpi_comm)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: mpi_comm
    
    integer :: ierror
    

    pf_comm%comm = mpi_comm       !! assign communicator
    
    
    !> assign number of processors
    call mpi_comm_size(mpi_comm, pf_comm%nproc, ierror)
    
    !>  assign procedure pointers
    pf_comm%post => pf_mpi_post
    pf_comm%recv => pf_mpi_recv
    pf_comm%send => pf_mpi_send
    pf_comm%wait => pf_mpi_wait
    pf_comm%broadcast   => pf_mpi_broadcast
    pf_comm%recv_status => pf_mpi_recv_status
    pf_comm%send_status => pf_mpi_send_status
  end subroutine pf_mpi_create


  !> Subroutine to set up the PFASST communicator.
  !! This should be called soon after adding levels to the PFASST controller 
  subroutine pf_mpi_setup(pf_comm, pf,ierror)
    use pf_mod_mpi, only: MPI_REQUEST_NULL
    use pf_mod_stop, only: pf_stop

    type(pf_comm_t),   intent(inout) :: pf_comm    !!  communicator 
    type(pf_pfasst_t), intent(inout) :: pf         !!  main pfasst structure
    integer,           intent(inout) :: ierror     !!  error flag

    !>  set the rank
    call mpi_comm_rank(pf_comm%comm, pf%rank, ierror)

    !>  allocate arrarys for and and receive requests
    allocate(pf_comm%recvreq(pf%nlevels),stat=ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierror)
    allocate(pf_comm%sendreq(pf%nlevels),stat=ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierror)    

    pf_comm%sendreq = MPI_REQUEST_NULL
    pf_comm%statreq = -66   !Tells the first send_status not to wait for previous one to arrive
  end subroutine pf_mpi_setup


  !> Subroutine to destroy the PFASST communicator.
  subroutine pf_mpi_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm

    deallocate(pf_comm%recvreq)
    deallocate(pf_comm%sendreq)
  end subroutine pf_mpi_destroy

  !>  Subroutine to post receive requests.
  subroutine pf_mpi_post(pf, level, tag, ierror, source)

    type(pf_pfasst_t), intent(in   ) :: pf
    class(pf_level_t), intent(inout) :: level   !!  level to send from
    integer,           intent(in   ) :: tag     !!  message tag
    integer,           intent(inout) :: ierror  !!  error flag
    integer,           intent(in)    :: source


    call mpi_irecv(level%recv, level%mpibuflen, myMPI_Datatype, &
                   source, tag, pf%comm%comm, pf%comm%recvreq(level%index), ierror)

  end subroutine pf_mpi_post


  !> Subroutine to send convergence status information
  subroutine pf_mpi_send_status(pf, tag,istatus,ierror, dest)
    use pf_mod_mpi, only: MPI_INTEGER, MPI_STATUS_SIZE, MPI_REQUEST_NULL

    type(pf_pfasst_t), intent(inout) :: pf        !!  main pfasst structure
    integer,           intent(in)    :: tag       !!  message tag
    integer,           intent(in) :: istatus      !!  status flag to send
    integer,           intent(inout) :: ierror    !!  error flag
    integer,            intent(in)    :: dest
    integer    ::  stat(MPI_STATUS_SIZE)

    integer :: message
    message = istatus
    

    if (pf%comm%statreq /= -66) then
       if (pf%debug) print*, 'DEBUG --',pf%rank, 'waiting in send_status with statreq',pf%comm%statreq
       call mpi_wait(pf%comm%statreq, stat, ierror)
       if (pf%debug) print*, 'DEBUG --',pf%rank, 'done waiting in send_status'
    end if

    call mpi_issend(message, 1, MPI_INTEGER, &
                    dest, tag, pf%comm%comm, pf%comm%statreq, ierror)

  end subroutine pf_mpi_send_status

  !> Subroutine to receive convergence status information
  subroutine pf_mpi_recv_status(pf, tag,istatus,ierror, source)
    use pf_mod_mpi, only: MPI_INTEGER, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf        !!  main pfasst structure
    integer,           intent(in)    :: tag       !!  message tag
    integer,           intent(inout) :: istatus   !!  status flag to receive
    integer,           intent(inout) :: ierror    !!  error flag
    integer,           intent(in)    :: source
    integer    ::  stat(MPI_STATUS_SIZE)

    integer :: message

    ! Get the message
    call mpi_recv(message, 1, MPI_INTEGER,source, tag, pf%comm%comm, stat, ierror)
    istatus=message

  end subroutine pf_mpi_recv_status


  !> Subroutine to send solutions
  subroutine pf_mpi_send(pf, level, tag, blocking,ierror, dest)
    use pf_mod_mpi, only:  MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf       !!  main pfasst structure
    class(pf_level_t), intent(inout) :: level    !!  level to send from
    integer,           intent(in   ) :: tag      !!  message tag
    logical,           intent(in   ) :: blocking !!  true if send is blocking
    integer,           intent(inout) :: ierror   !!  error flag
    integer,           intent(in)    :: dest

    integer ::  stat(MPI_STATUS_SIZE)

    
    if (blocking) then
       call mpi_send(level%send, level%mpibuflen, myMPI_Datatype, &
                     dest, tag, pf%comm%comm, stat, ierror)
    else

       call mpi_isend(level%send, level%mpibuflen, myMPI_Datatype, &
                      dest, tag, pf%comm%comm, pf%comm%sendreq(level%index), ierror)
    end if
  end subroutine pf_mpi_send

  !> Subroutine to receive solutions
  !! Note when blocking == .false. this is actually a wait because the
  !! nonblocking receive  should have already been posted
  subroutine pf_mpi_recv(pf, level, tag, blocking, ierror, source)
    use pf_mod_mpi, only:  MPI_STATUS_SIZE
    type(pf_pfasst_t), intent(inout) :: pf     !!  main pfasst structure
    class(pf_level_t), intent(inout) :: level  !!  level to recieve into
    integer,           intent(in   ) :: tag    !!  message tag
    logical,           intent(in   ) :: blocking  !!  true if receive is blocking
    integer,           intent(inout) :: ierror  !!  error flag
    integer,           intent(in)    :: source
    integer ::  stat(MPI_STATUS_SIZE)

    if(blocking) then
       call mpi_recv(level%recv, level%mpibuflen, myMPI_Datatype, &
            source, tag, pf%comm%comm, stat, ierror)
    else
       call mpi_wait(pf%comm%recvreq(level%index), stat, ierror)
    end if
  end subroutine pf_mpi_recv

  !
  ! Misc
  !
  subroutine pf_mpi_wait(pf, level, ierror)
    use pf_mod_mpi, only: MPI_STATUS_SIZE
    type(pf_pfasst_t), intent(in) :: pf           !!  main pfasst structure
    integer,           intent(in) :: level        !!  level on which to wait
    integer,           intent(inout) :: ierror    !!  error flag
    integer ::  stat(MPI_STATUS_SIZE)
    call mpi_wait(pf%comm%sendreq(level), stat, ierror)
  end subroutine pf_mpi_wait

  subroutine pf_mpi_broadcast(pf, y, nvar, root,ierror)
    type(pf_pfasst_t), intent(inout) :: pf      !!  main pfasst structure
    integer,           intent(in)    :: nvar    !!  size of data to broadcast
    real(pfdp),        intent(in)    :: y(nvar) !!  data to broadcast
    integer,           intent(in)    :: root    !!  rank of broadcaster
    integer,           intent(inout) :: ierror  !!  error flag
    call mpi_bcast(y, nvar, myMPI_Datatype, root, pf%comm%comm, ierror)
  end subroutine pf_mpi_broadcast

end module pf_mod_comm_mpi

