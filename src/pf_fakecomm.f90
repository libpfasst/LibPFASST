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

! This module implements fake communications.

module pf_mod_comm_fake
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_pfasst
  use iso_c_binding

  implicit none

contains

  ! Create a fake PFASST communicator (call only once)
  subroutine pf_fake_create(pf_comm, nprocs)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: nprocs

    pf_comm%nproc = nprocs
    allocate(pf_comm%pfs(nprocs))

    pf_comm%post => pf_fake_post
    pf_comm%recv => pf_fake_recv
    pf_comm%send => pf_fake_send
    pf_comm%wait => pf_fake_wait
    pf_comm%broadcast => pf_fake_broadcast
  end subroutine pf_fake_create

  ! Setup
  subroutine pf_fake_setup(pf_comm, pf)
    type(pf_comm_t), intent(inout) :: pf_comm
    type(pf_pfasst_t), intent(inout), target  :: pf

    pf_comm%pfs(pf%rank+1) = c_loc(pf)
  end subroutine pf_fake_setup

  ! Retrieve the PFASST object associated with the given rank
  subroutine pf_fake_get(pf_comm, rank, pf)
    type(pf_comm_t),   intent(in)  :: pf_comm
    integer,           intent(in)  :: rank 
    type(pf_pfasst_t), intent(out), pointer :: pf

    call c_f_pointer(pf_comm%pfs(rank+1), pf)
  end subroutine pf_fake_get

  ! Destroy (call only once)
  subroutine pf_fake_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm
    integer :: t

    type(pf_pfasst_t), pointer :: pf

    do t = 1, pf_comm%nproc
       call c_f_pointer(pf_comm%pfs(t), pf)
       call pf_pfasst_destroy(pf)
       deallocate(pf)
    end do

    deallocate(pf_comm%pfs)
  end subroutine pf_fake_destroy

  ! Post
  subroutine pf_fake_post(pf, level, tag)
    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    ! this is intentionally empty
  end subroutine pf_fake_post

  ! Receive
  subroutine pf_fake_recv(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    type(pf_pfasst_t), pointer :: from

    call start_timer(pf, TRECEIVE)

    if (pf%rank > 0) then
       call pf_fake_get(pf%comm, pf%rank-1, from)
       level%q0 = from%levels(level%level)%send
    end if

    call end_timer(pf, TRECEIVE)
  end subroutine pf_fake_recv

  ! Send
  subroutine pf_fake_send(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    call start_timer(pf, TSEND)

    if (pf%rank < pf%comm%nproc-1) then
       call level%encap%pack(level%send, level%qend)
    end if

    call end_timer(pf, TSEND)
  end subroutine pf_fake_send

  ! Wait
  subroutine pf_fake_wait(pf, level)
    type(pf_pfasst_t), intent(in) :: pf
    integer,           intent(in) :: level
    stop "FAKE WAIT NOT IMPLEMENTED YET"
  end subroutine pf_fake_wait

  ! Broadcast
  subroutine pf_fake_broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(inout) :: pf
    real(kind=8),      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root
    stop "FAKE BROADCAST NOT IMPLEMENTED YET"
  end subroutine pf_fake_broadcast

end module pf_mod_comm_fake
