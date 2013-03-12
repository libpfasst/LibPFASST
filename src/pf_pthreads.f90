!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! This module implements PTHREADS communications.

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

  ! pthread interfaces
  interface
     integer(c_int) function pthread_self() &
          bind(c, name='pthread_self')
       use iso_c_binding
     end function pthread_self
  end interface

  interface
     integer(c_int) function pthread_join(thread, valueptr) &
          bind(c, name='pthread_join')
       use iso_c_binding
       integer(c_long), intent(in), value :: thread
       type(c_ptr),     intent(in), value :: valueptr
     end function pthread_join
  end interface

  interface
     integer(c_long) function pthread_create(thread, attr, start, arg) &
          bind(c, name='pthread_create')
       use iso_c_binding
       type(c_ptr), intent(in), value :: thread, attr, start, arg
     end function pthread_create
  end interface

  interface
     subroutine pthread_exit(retval) bind(c, name='pthread_exit')
       use iso_c_binding
       type(c_ptr), intent(in), value :: retval
     end subroutine pthread_exit
  end interface

  ! pfasst pthreads interfaces
  interface
     type(c_ptr) function pf_pth_create() bind(c, name='pf_pth_create')
       use iso_c_binding
     end function pf_pth_create
  end interface

  interface
     subroutine pf_pth_destroy(pth) bind(c, name='pf_pth_destroy')
       use iso_c_binding
       type(c_ptr), intent(in), value :: pth
     end subroutine pf_pth_destroy
  end interface

  interface
     subroutine pf_pth_wait_send(pth, tag) bind(c, name='pf_pth_wait_send')
       use iso_c_binding
       type(c_ptr), intent(in), value :: pth
       integer(c_int), intent(in), value :: tag
     end subroutine pf_pth_wait_send
  end interface

  interface
     subroutine pf_pth_set_send(pth, tag) bind(c, name='pf_pth_set_send')
       use iso_c_binding
       type(c_ptr), intent(in), value :: pth
       integer(c_int), intent(in), value :: tag
     end subroutine pf_pth_set_send
  end interface

  interface
     subroutine pf_pth_wait_recv(pth, tag) bind(c, name='pf_pth_wait_recv')
       use iso_c_binding
       type(c_ptr), intent(in), value :: pth
       integer(c_int), intent(in), value :: tag
     end subroutine pf_pth_wait_recv
  end interface

  interface
     subroutine pf_pth_set_recv(pth, tag) bind(c, name='pf_pth_set_recv')
       use iso_c_binding
       type(c_ptr), intent(in), value :: pth
       integer(c_int), intent(in), value :: tag
     end subroutine pf_pth_set_recv
  end interface

  interface
     subroutine pf_pth_lock(pth) bind(c, name='pf_pth_lock')
       use iso_c_binding
       type(c_ptr), intent(in), value :: pth
     end subroutine pf_pth_lock
  end interface

  interface
     subroutine pf_pth_unlock(pth) bind(c, name='pf_pth_unlock')
       use iso_c_binding
       type(c_ptr), intent(in), value :: pth
     end subroutine pf_pth_unlock
  end interface

contains

  ! Create a PTHREADS based PFASST communicator (call only once)
  !
  ! This is not thread safe.
  subroutine pf_comm_create(pf_comm, nthreads, nlevels)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: nthreads, nlevels

    integer :: t, l

    pf_comm%nproc = nthreads
    allocate(pf_comm%pfs(0:nthreads-1))
    allocate(pf_comm%pfpth(0:nthreads-1,nlevels))

    do t = 0, nthreads-1
       do l = 1, nlevels
          pf_comm%pfpth(t,l) = pf_pth_create()
       end do
    end do
  end subroutine pf_comm_create

  ! Setup
  subroutine pf_comm_setup(pf_comm, pf)
    type(pf_comm_t), intent(inout) :: pf_comm
    type(pf_pfasst_t), intent(inout), target  :: pf

    integer :: n

    n = pf%rank
    pf_comm%pfs(n) = c_loc(pf)
  end subroutine pf_comm_setup

  ! Retrieve the PFASST object associated with the given rank
  subroutine pf_comm_get(pf_comm, rank, pf)
    type(pf_comm_t),   intent(in)  :: pf_comm
    integer,           intent(in)  :: rank 
    type(pf_pfasst_t), intent(out), pointer :: pf

    call c_f_pointer(pf_comm%pfs(rank), pf)
  end subroutine pf_comm_get

  ! Destroy (call only once)
  subroutine pf_comm_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm
    integer :: t, l

    type(pf_pfasst_t), pointer :: pf

    do t = 0, size(pf_comm%pfs)-1
       do l = 1, size(pf_comm%pfpth(t,:))
          call pf_pth_destroy(pf_comm%pfpth(t,l))
       end do
    end do

    do t = 0, pf_comm%nproc-1
       call c_f_pointer(pf_comm%pfs(t), pf)
       call destroy(pf)
       deallocate(pf)
    end do

    deallocate(pf_comm%pfs)
    deallocate(pf_comm%pfpth)
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
       pth = pf%comm%pfpth(from%rank, level%level)

       call pf_pth_wait_send(pth, tag)

       call pf_pth_lock(pth)
       level%q0 = from%levels(level%level)%send
       call pf_pth_unlock(pth)

       call pf_pth_set_recv(pth, 0)
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
       pth = pf%comm%pfpth(pf%rank, level%level)

       call pf_pth_wait_recv(pth, 0)

       call pf_pth_lock(pth)
       call pack(level%send, level%qend)
       call pf_pth_unlock(pth)

       call pf_pth_set_recv(pth, tag)
       call pf_pth_set_send(pth, tag)
    end if

    call end_timer(pf, TSEND)
  end subroutine send

  ! Wait
  subroutine wait(pf, level)
    type(pf_pfasst_t), intent(in) :: pf
    integer,           intent(in) :: level
  end subroutine wait

  ! Broadcast
  subroutine broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(in) :: pf
    real(kind=8),      intent(in) :: y(nvar)
    integer,           intent(in) :: nvar, root

    stop
  end subroutine broadcast

end module pf_mod_comm
