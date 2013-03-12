!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! This module implements communications on the fake communicator

module pf_mod_comm
  use encap
  use pf_mod_dtype
  use pf_mod_pfasst
  use pf_mod_timer
  use iso_c_binding

  implicit none

  interface create
     module procedure pf_comm_create
  end interface create

  interface setup
     module procedure pf_comm_setup
  end interface setup

  interface destroy
     module procedure pf_comm_destroy
  end interface destroy

contains

  ! Create a PTHREADS based PFASST communicator (call only once)
  subroutine pf_comm_create(pf_comm, nprocs)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: nprocs

    integer :: t, l

    pf_comm%nproc = nprocs
    allocate(pf_comm%pfs(nprocs))
  end subroutine pf_comm_create

  ! Setup
  subroutine pf_comm_setup(pf_comm, pf)
    type(pf_comm_t), intent(inout) :: pf_comm
    type(pf_pfasst_t), intent(inout), target  :: pf

    integer :: n

    n = pf%rank
    pf_comm%pfs(n+1) = c_loc(pf)
  end subroutine pf_comm_setup

  ! Retrieve the PFASST object associated with the given rank
  subroutine pf_comm_get(pf_comm, rank, pf)
    type(pf_comm_t),   intent(in)  :: pf_comm
    integer,           intent(in)  :: rank 
    type(pf_pfasst_t), intent(out), pointer :: pf

    call c_f_pointer(pf_comm%pfs(rank+1), pf)
  end subroutine pf_comm_get

  ! Destroy (call only once)
  subroutine pf_comm_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm
    integer :: t, l

    type(pf_pfasst_t), pointer :: pf

    do t = 1, pf_comm%nproc
       call c_f_pointer(pf_comm%pfs(t), pf)
       call destroy(pf)
       deallocate(pf)
    end do

    deallocate(pf_comm%pfs)
  end subroutine pf_comm_destroy

  ! Post
  subroutine post(pf, level, tag)
    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
  end subroutine post

  ! Receive
  subroutine recv(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    type(pf_pfasst_t), pointer :: from
    type(c_ptr) :: pth

    call start_timer(pf, TRECEIVE)

    if (pf%rank > 0) then
       call pf_comm_get(pf%comm, pf%rank-1, from)
       level%q0 = from%levels(level%level)%send
    end if

    call end_timer(pf, TRECEIVE)
  end subroutine recv

  ! Send
  subroutine send(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    type(c_ptr)       :: pth

    call start_timer(pf, TSEND)

    if (pf%rank < pf%comm%nproc-1) then
       call pack(level%send, level%qend)
    end if

    call end_timer(pf, TSEND)
  end subroutine send

  ! Broadcast
  subroutine broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(in) :: pf
    real(pfdp),        intent(in) :: y(nvar)
    integer,           intent(in) :: nvar, root

    stop
  end subroutine broadcast

  ! Wait
  subroutine wait(pf, level)
    type(pf_pfasst_t), intent(in) :: pf
    integer,           intent(in) :: level
  end subroutine wait

end module pf_mod_comm
