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
  use pf_mod_mpi
  use pf_mod_dtype
  use pf_mod_timer

  implicit none

  interface create
     module procedure pf_comm_mpi_create
  end interface create

  interface setup
     module procedure pf_comm_mpi_setup
  end interface setup

  interface destroy
     module procedure pf_comm_mpi_destroy
  end interface destroy

contains

  ! Create an MPI based PFASST communicator
  subroutine pf_comm_mpi_create(pf_comm, mpi_comm)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: mpi_comm

    integer :: ierror

    pf_comm%comm = mpi_comm
    call mpi_comm_size(mpi_comm, pf_comm%nproc, ierror)
  end subroutine pf_comm_mpi_create

  ! Setup
  subroutine pf_comm_mpi_setup(pf_comm, pf)
    use pf_mod_mpi, only: MPI_REQUEST_NULL

    type(pf_comm_t),   intent(inout) :: pf_comm
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: ierror

    call mpi_comm_rank(pf_comm%comm, pf%rank, ierror)

    allocate(pf_comm%recvreq(pf%nlevels))
    allocate(pf_comm%sendreq(pf%nlevels))

    pf_comm%sendreq = MPI_REQUEST_NULL
  end subroutine pf_comm_mpi_setup

  ! Destroy
  subroutine pf_comm_mpi_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm

    deallocate(pf_comm%recvreq)
    deallocate(pf_comm%sendreq)
  end subroutine pf_comm_mpi_destroy

  ! Post
  subroutine pf_comm_mpi_post(pf, level, tag)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag

    integer :: ierror

    if (pf%comm%nproc > 1 .and. pf%rank > 0) then
       call mpi_irecv(level%recv, level%nvars, MPI_REAL8, &
            pf%rank-1, tag, pf%comm%comm, pf%comm%recvreq(level%level), ierror)
    end if
  end subroutine pf_comm_mpi_post

  ! Receive
  subroutine pf_comm_mpi_recv(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call start_timer(pf, TRECEIVE + level%level - 1)

    if (pf%rank > 0) then
       if (blocking) then
          call mpi_recv(level%recv, level%nvars, MPI_REAL8, &
               pf%rank-1, tag, pf%comm%comm, stat, ierror)
       else
          call mpi_wait(pf%comm%recvreq(level%level), stat, ierror)
       end if

       if (ierror .ne. 0) then
          print *, 'WARNING: MPI ERROR DURING RECEIVE', ierror
       end if

       level%q0 = level%recv
    end if

    call end_timer(pf, TRECEIVE + level%level - 1)
  end subroutine pf_comm_mpi_recv

  ! Send
  subroutine pf_comm_mpi_send(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call start_timer(pf, TSEND + level%level - 1)

    if (pf%rank < pf%comm%nproc-1) then

       if (blocking) then
          call pack(level%send, level%qend)
          call mpi_send(level%send, level%nvars, MPI_REAL8, &
               pf%rank+1, tag, pf%comm%comm, stat, ierror)
       else
          call mpi_wait(pf%comm%sendreq(level%level), stat, ierror)
          call pack(level%send, level%qend)
          call mpi_isend(level%send, level%nvars, MPI_REAL8, &
               pf%rank+1, tag, pf%comm%comm, pf%comm%sendreq(level%level), ierror)
       end if
    end if

    call end_timer(pf, TSEND + level%level - 1)
  end subroutine pf_comm_mpi_send

  ! Send
  subroutine pf_comm_mpi_wait(pf, level)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(in) :: pf
    integer,           intent(in) :: level

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call mpi_wait(pf%comm%sendreq(level), stat, ierror)
  end subroutine pf_comm_mpi_wait


  ! Broadcast
  subroutine pf_comm_mpi_broadcast(pf, y, nvar, root)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root

    integer :: ierror

    call start_timer(pf, TSEND)

    call mpi_bcast(y, nvar, MPI_REAL8, root, pf%comm%comm, ierror)

    call end_timer(pf, TSEND)
  end subroutine pf_comm_mpi_broadcast

end module pf_mod_comm_mpi
