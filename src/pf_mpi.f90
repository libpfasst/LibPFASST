!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

! This module implements MPI communications.

module pf_mod_mpi
  include "mpif.h"
end module pf_mod_mpi

module pf_mod_comm_mpi
  use pf_mod_dtype
  use pf_mod_timer
  implicit none
contains

  subroutine pf_mpi_create(pf_comm, mpi_comm)
    !
    ! Create an MPI based PFASST communicator using the MPI
    ! communicator *mpi_comm*.
    !
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: mpi_comm

    integer :: ierror

    pf_comm%comm = mpi_comm
    call mpi_comm_size(mpi_comm, pf_comm%nproc, ierror)

    pf_comm%post => pf_mpi_post
    pf_comm%recv => pf_mpi_recv
    pf_comm%send => pf_mpi_send
    pf_comm%wait => pf_mpi_wait
    pf_comm%broadcast   => pf_mpi_broadcast
    pf_comm%recv_status => pf_mpi_recv_status
    pf_comm%send_status => pf_mpi_send_status
    pf_comm%recv_nmoved => pf_mpi_recv_nmoved
    pf_comm%send_nmoved => pf_mpi_send_nmoved
  end subroutine pf_mpi_create

  subroutine pf_mpi_setup(pf_comm, pf)
    !
    ! Setup the PFASST communicator.
    !
    ! This should be called soon after adding levels to the PFASST
    ! controller **pf**.
    !
    use pf_mod_mpi, only: MPI_REQUEST_NULL

    type(pf_comm_t),   intent(inout) :: pf_comm
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: ierror

    call mpi_comm_rank(pf_comm%comm, pf%rank, ierror)

    allocate(pf_comm%recvreq(pf%nlevels))
    allocate(pf_comm%sendreq(pf%nlevels))

    pf_comm%sendreq = MPI_REQUEST_NULL
    pf_comm%statreq = -66
  end subroutine pf_mpi_setup

  subroutine pf_mpi_destroy(pf_comm)
    !
    ! Destroy the PFASST communicator.
    !
    type(pf_comm_t), intent(inout) :: pf_comm

    deallocate(pf_comm%recvreq)
    deallocate(pf_comm%sendreq)
  end subroutine pf_mpi_destroy

  subroutine pf_mpi_post(pf, level, tag)
    !
    ! Post receive requests.
    !
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag

    integer :: ierror

    if (pf%rank /= pf%state%first  .and. &
         pf%state%pstatus == PF_STATUS_ITERATING) then
       call mpi_irecv(level%recv, level%nvars, MPI_REAL8, &
            modulo(pf%rank-1, pf%comm%nproc), tag, pf%comm%comm, &
            pf%comm%recvreq(level%level), ierror)
    end if
  end subroutine pf_mpi_post

  subroutine pf_mpi_recv(pf, level, tag, blocking)
    !
    ! Receive new initial condition.
    !
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call start_timer(pf, TRECEIVE + level%level - 1)

    if (pf%rank /= pf%state%first .and. &
         pf%state%pstatus == PF_STATUS_ITERATING) then

       if (blocking) then
          call mpi_recv(level%recv, level%nvars, MPI_REAL8, &
               modulo(pf%rank-1, pf%comm%nproc), tag, pf%comm%comm, stat, ierror)
       else
          call mpi_wait(pf%comm%recvreq(level%level), stat, ierror)
       end if

       if (ierror .ne. 0) then
          print *, pf%rank, 'warning: mpi error during receive', ierror
       end if

       level%q0 = level%recv
    end if

    call end_timer(pf, TRECEIVE + level%level - 1)
  end subroutine pf_mpi_recv

  ! Receive status
  subroutine pf_mpi_recv_status(pf, tag)
    use pf_mod_mpi, only: MPI_INTEGER4, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag

    integer    :: ierror, stat(MPI_STATUS_SIZE)
    integer(4) :: message(8)

    pf%state%nmoved = 0

    if (pf%rank /= pf%state%first) then

       call mpi_recv(message, 8, MPI_INTEGER4, &
            modulo(pf%rank-1, pf%comm%nproc), tag, pf%comm%comm, stat, ierror)

       pf%state%pstatus = message(7)
       pf%state%nmoved  = message(8)

       if (ierror .ne. 0) then
          print *, pf%rank, 'warning: mpi error during receive status', ierror
       end if

    end if

  end subroutine pf_mpi_recv_status

  ! Receive status
  subroutine pf_mpi_recv_nmoved(pf, tag)
    use pf_mod_mpi
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag

    integer    :: ierror, stat(MPI_STATUS_SIZE)
    integer(4) :: src, message(8)

    call mpi_probe(MPI_ANY_SOURCE, tag, pf%comm%comm, stat, ierror)
    src = stat(MPI_SOURCE)

    call mpi_recv(message, 8, MPI_INTEGER4, &
         src, tag, pf%comm%comm, stat, ierror)
    pf%state%nmoved = message(8)

  end subroutine pf_mpi_recv_nmoved

  ! Send
  subroutine pf_mpi_send(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call start_timer(pf, TSEND + level%level - 1)

    if (pf%rank /= pf%state%last &
         .and. pf%state%status == PF_STATUS_ITERATING) then

       if (blocking) then
          call level%encap%pack(level%send, level%qend)
          call mpi_send(level%send, level%nvars, MPI_REAL8, &
               modulo(pf%rank+1, pf%comm%nproc), tag, pf%comm%comm, &
               stat, ierror)
       else
          call mpi_wait(pf%comm%sendreq(level%level), stat, ierror)
          call level%encap%pack(level%send, level%qend)
          call mpi_isend(level%send, level%nvars, MPI_REAL8, &
               modulo(pf%rank+1, pf%comm%nproc), tag, pf%comm%comm, &
               pf%comm%sendreq(level%level), ierror)
       end if
    end if

    call end_timer(pf, TSEND + level%level - 1)

  end subroutine pf_mpi_send

  ! Send status
  subroutine pf_mpi_send_status(pf, tag)
    use pf_mod_mpi!, only: MPI_INTEGER, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag

    integer    :: ierror, stat(MPI_STATUS_SIZE)
    integer(4) :: message(8)

    if (pf%rank /= pf%state%last) then

       ! on orga there is some weird issue with send/recv status: the
       ! first two integer4's always get set to zero on the receiving
       ! end.  hence we use 8 integer4's and put the message in the
       ! last two slots.

       message = 666
       message(7) = pf%state%status
       message(8) = pf%state%nmoved

       if (pf%window == PF_WINDOW_RING .and. pf%state%status == PF_STATUS_CONVERGED) then
          message(8) = message(8) + 1
       end if

       if (pf%comm%statreq /= -66) then
          call mpi_wait(pf%comm%statreq, stat, ierror)
       end if
       call mpi_issend(message, 8, MPI_INTEGER4, &
            modulo(pf%rank+1, pf%comm%nproc), tag, pf%comm%comm, pf%comm%statreq, ierror)
    end if

  end subroutine pf_mpi_send_status

  subroutine pf_mpi_send_nmoved(pf, tag)
    use pf_mod_mpi

    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag

    integer    :: r, ierror, stat(MPI_STATUS_SIZE), sendreq(pf%state%nmoved)
    integer(4) :: message(8)

    message = 666
    message(8) = pf%state%nmoved

    do r = 1, pf%state%nmoved
       call mpi_isend(message, 8, MPI_INTEGER4, &
            modulo(pf%rank-r, pf%comm%nproc), tag, pf%comm%comm, sendreq(r), ierror)
    end do

    do r = 1, pf%state%nmoved
       call mpi_wait(sendreq(r), stat, ierror)
    end do

  end subroutine pf_mpi_send_nmoved

  ! Send
  subroutine pf_mpi_wait(pf, level)
    use pf_mod_mpi, only: MPI_STATUS_SIZE
    type(pf_pfasst_t), intent(in) :: pf
    integer,           intent(in) :: level
    integer :: ierror, stat(MPI_STATUS_SIZE)
    call mpi_wait(pf%comm%sendreq(level), stat, ierror)
  end subroutine pf_mpi_wait

  ! Broadcast
  subroutine pf_mpi_broadcast(pf, y, nvar, root)
    use pf_mod_mpi, only: MPI_REAL8
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root
    integer :: ierror
    call start_timer(pf, TSEND)
    call mpi_bcast(y, nvar, MPI_REAL8, root, pf%comm%comm, ierror)
    call end_timer(pf, TSEND)
  end subroutine pf_mpi_broadcast

end module pf_mod_comm_mpi
