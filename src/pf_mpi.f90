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
   subroutine pf_mpi_setup(pf_comm, pf,ierror, pf_diag_comm)
      use pf_mod_mpi, only: MPI_REQUEST_NULL
      use pf_mod_stop, only: pf_stop

      type(pf_comm_t),   intent(inout) :: pf_comm    !!  communicator
      type(pf_pfasst_t), intent(inout) :: pf         !!  main pfasst structure
      integer,           intent(inout) :: ierror     !!  error flag
      type(pf_comm_t),   intent(inout), optional :: pf_diag_comm    !!  diagonal sdc communicator

      !>  set the rank
      call mpi_comm_rank(pf_comm%comm, pf%rank, ierror)

      !>  allocate arrarys for and and receive requests
      allocate(pf_comm%recvreq(pf%nlevels),stat=ierror)
      if (ierror /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierror)
      allocate(pf_comm%sendreq(pf%nlevels),stat=ierror)
      if (ierror /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierror)

      pf_comm%sendreq = MPI_REQUEST_NULL
      pf_comm%statreq = MPI_REQUEST_NULL   !Tells the first send_status not to wait for previous one to arrive

      ! add setup for diagonal comm if available
      if (present(pf_diag_comm)) then
         !>  set the rank
         call mpi_comm_rank(pf_diag_comm%comm, pf%rank_diag, ierror)

         !>  allocate arrarys for and and receive requests
         allocate(pf_diag_comm%recvreq(pf%nlevels),stat=ierror)
         if (ierror /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierror)
         allocate(pf_diag_comm%sendreq(pf%nlevels),stat=ierror)
         if (ierror /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierror)

         pf_diag_comm%sendreq = MPI_REQUEST_NULL
         pf_diag_comm%statreq = MPI_REQUEST_NULL   !Tells the first send_status not to wait for previous one to arrive
      end if

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
      use pf_mod_mpi, only: MPI_INTEGER4, MPI_STATUS_SIZE, MPI_REQUEST_NULL,MPI_SIZEOF

      type(pf_pfasst_t), intent(inout) :: pf        !!  main pfasst structure
      integer,           intent(in)    :: tag       !!  message tag
      integer,           intent(in) :: istatus      !!  status flag to send
      integer,           intent(inout) :: ierror    !!  error flag
      integer,            intent(in)    :: dest
      integer    ::  stat(MPI_STATUS_SIZE)

      integer :: message(1)
      message = istatus


      if (pf%comm%statreq /= MPI_REQUEST_NULL) then
         if (pf%debug) print*, 'DEBUG --',pf%rank, 'waiting in send_status with statreq',pf%comm%statreq
         call mpi_wait(pf%comm%statreq, stat, ierror)
         if (pf%debug) print*, 'DEBUG --',pf%rank, 'done waiting in send_status'
      end if

      if (pf%debug) print*, 'DEBUG --',pf%rank, 'begin issend_status', istatus, message,pf%comm%statreq
      call mpi_issend(message, 1, MPI_INTEGER4, &
         dest, tag, pf%comm%comm, pf%comm%statreq, ierror)
      if (pf%debug) print*, 'DEBUG --',pf%rank, 'end issend  status', istatus, message,pf%comm%statreq


   end subroutine pf_mpi_send_status

   !> Subroutine to receive convergence status information
   subroutine pf_mpi_recv_status(pf, tag,istatus,ierror, source)
      use pf_mod_mpi, only: MPI_INTEGER4, MPI_STATUS_SIZE

      type(pf_pfasst_t), intent(inout) :: pf        !!  main pfasst structure
      integer,           intent(in)    :: tag       !!  message tag
      integer,           intent(inout) :: istatus   !!  status flag to receive
      integer,           intent(inout) :: ierror    !!  error flag
      integer,           intent(in)    :: source
      integer    ::  stat(MPI_STATUS_SIZE)

      integer :: message(1)

      ! Get the message
      call mpi_recv(message, 1, MPI_INTEGER4,source, tag, pf%comm%comm, stat, ierror)
      istatus=message(1)
      if (pf%debug) print *,'DEBUG- rank=',pf%rank,'in recv_status, istatus,message', istatus, message
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
         call mpi_ssend(level%send, level%mpibuflen, myMPI_Datatype, &
            dest, tag, pf%comm%comm, stat, ierror)
      else

         call mpi_issend(level%send, level%mpibuflen, myMPI_Datatype, &
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
      if(pf%debug)  print *,pf%rank,'rank broadcasting from rank=',root
      call mpi_bcast(y, nvar, myMPI_Datatype, root, pf%comm%comm, ierror)
   end subroutine pf_mpi_broadcast

   !> core subroutine to split - only works on plain MPI comms
   subroutine pf_mpi_orthosplit(size_comm1, size_comm2, comm1, comm2, color1, color2, parent_comm, &
                                size_comm3, comm3, color3)
      use pf_mod_mpi, only: mpi_comm_size, mpi_comm_rank, mpi_comm_split, mpi_barrier, MPI_COMM_SELF
      integer, intent(out) :: comm1, comm2
      integer, intent(out) :: color1, color2
      integer, intent(in) :: size_comm1, size_comm2
      integer, intent(in) :: parent_comm
      ! optional 3D extension
      integer, intent(in),  optional :: size_comm3
      integer, intent(out), optional :: comm3
      integer, intent(out), optional :: color3 
      ! local variables
      integer :: pSize, pRank, err, required_size
      integer :: rank1, size1, rank2, size2, rank3, size3

      !> check size
      call mpi_comm_size(parent_comm, pSize, err)
      call mpi_comm_rank(parent_comm, pRank,  err)
      
      !> sanity check on size_comm1, size_comm2
      if (present(size_comm3)) then
         required_size = size_comm1 * size_comm2 * size_comm3
      else
         required_size = size_comm1 * size_comm2
      end if
      ! 
      if (pSize .ne. required_size) then
         print '(a)', 'ERROR: pf_mpi_orthosplit_comm_core: processor number mismatch.'
         print '(a,i4,a,i4)', ' Expecting ', required_size, ' MPI processors but received ', pSize
         stop
      end if

      !> do the actual splitting
      ! Think of ranks arranged in a (size_comm2 x size_comm1) grid:
      !  pSize = 12, size_comm1=4, size_comm2=3 
      !      comm1 ranks (size 4)
      !       r0   r1   r2   r3
      !     +----+----+----+----+
      !  c  |  0 |  1 |  2 |  3 |  <- comm2 rank 0 (size 3)
      !  o  +----+----+----+----+
      !  m  |  4 |  5 |  6 |  7 |  <- comm2 rank 1
      !  m  +----+----+----+----+
      !  2  |  8 |  9 | 10 | 11 |  <- comm2 rank 2
      !     +----+----+----+----+
      !
      ! 3D example, pSize = 23, size_comm1=4, size_comm2=3, size_comm3=2 
      !         layer 0 (color3=0):          layer 1 (color3=1):
      !          comm1 (size 4)               comm1 (size 4)
      !          c0  c1  c2  c3               c0  c1  c2  c3
      !         +---+---+---+---+            +---+---+---+---+
      !    c  c0| 0 | 1 | 2 | 3 |       c  c0|12 |13 |14 |15 |
      !    o    +---+---+---+---+       o    +---+---+---+---+
      !    m  c1| 4 | 5 | 6 | 7 |       m  c1|16 |17 |18 |19 |
      !    m    +---+---+---+---+       m    +---+---+---+---+
      !    2  c2| 8 | 9 |10 |11 |       2  c2|20 |21 |22 |23 |
      !         +---+---+---+---+            +---+---+---+---+
      ! comm1: groups of size_comm1 consecutive ranks (rows)
      color1 = pRank / size_comm1
      call mpi_comm_split(parent_comm, color1, pRank, comm1, err)

      ! comm2: interleaved groups (columns)
      color2 = (pRank / (size_comm1 * size_comm2)) * size_comm1 + mod(pRank, size_comm1)
      call mpi_comm_split(parent_comm, color2, pRank, comm2, err)

      ! comm3: optional third dimension (layers)
      if (present(size_comm3) .and. present(comm3) .and. present(color3)) then
         color3 = mod(pRank, (size_comm1 * size_comm2))
         call mpi_comm_split(parent_comm, color3, pRank, comm3, err)
      end if

      !> rest is just some verbosity output
      call mpi_comm_rank(comm1, rank1, err)
      call mpi_comm_size(comm1, size1, err)
      call mpi_comm_rank(comm2, rank2, err)
      call mpi_comm_size(comm2, size2, err)

      call mpi_barrier(parent_comm, err)
      
      if (pRank == 0) then
         print '(a)', ''
         print '(a)', '============ pf_mpi_orthosplit results ============'
         if (present(size_comm3)) then
            print '(a,i4,a,i4,a,i4,a,i4)', &
               ' parent size=', pSize, &
               ' comm1 size=', size_comm1, &
               ' comm2 size=', size_comm2, &
               ' comm3 size=', size_comm3
         else
            print '(a,i4,a,i4,a,i4)', &
               ' parent size=', pSize, &
               ' comm1 size=', size_comm1, &
               ' comm2 size=', size_comm2
         end if
         print '(a)', ''
         print '(a)', ' parent_rank | comm1_rank/size | comm2_rank/size' // &
                      merge(' | comm3_rank/size', '                  ', present(size_comm3))
         print '(a)', ' ----------- | --------------- | ---------------' // &
                      merge(' | ---------------', '                  ', present(size_comm3))
      end if

      call mpi_barrier(parent_comm, err)
      flush(6)
      call sleep(pRank/10)
      if (present(size_comm3) .and. present(comm3) .and. present(color3)) then
         call mpi_comm_rank(comm3, rank3, err)
         call mpi_comm_size(comm3, size3, err)
         print '(a,i6,a,i6,a,i4,a,i6,a,i4,a,i6,a,i4)', &
            '  ', pRank, '     | ', &
            rank1, ' / ', size1, '   | ', &
            rank2, ' / ', size2, '   | ', &
            rank3, ' / ', size3
      else
         print '(a,i6,a,i6,a,i4,a,i6,a,i4)', &
            '  ', pRank, '     | ', &
            rank1, ' / ', size1, '   | ', &
            rank2, ' / ', size2
      end if

      call mpi_barrier(parent_comm, err)
      if (pRank == 0) print '(a)', '================================================='

   end subroutine pf_mpi_orthosplit

end module pf_mod_comm_mpi

