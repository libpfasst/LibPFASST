!!  Communication wrappers
!
! This file is part of LIBPFASST.
!
!> Module of communication wrappers
module pf_mod_comm
  use pf_mod_pfasst

  implicit none
contains


  !>  Subroutine to post a receive request for a new initial condition to be received after doing some work
  subroutine pf_post(pf, level, tag, direction)
    type(pf_pfasst_t), intent(in)    :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: direction
    integer                          :: dir
    integer ::  ierror , source


    if (pf%comm%nproc .eq. 1) return

    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction

    ierror = 0
    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                  .and. dir == 1) then
       if (pf%debug) print  '("DEBUG-rank=", I5, " begin post, tag=",I8, " pf%state%pstatus=", I2)',pf%rank,tag,pf%state%pstatus           
       source=pf%rank-1
       call pf%comm%post(pf, level, tag, ierror, source)
       if (pf%debug) print  '("DEBUG-rank=", I5, " end post, tag=",I8, " pf%state%pstatus=", I2)',pf%rank,tag,pf%state%pstatus       
    elseif (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING &
         .and. dir == 2) then
       if (pf%debug) print  '("DEBUG-rank=", I5, " begin post, tag=",I8, " pf%state%pstatus=", I2)',pf%rank,tag,pf%state%pstatus    
       source=pf%rank+1
       call pf%comm%post(pf, level, tag, ierror, source)
       if (pf%debug) print  '("DEBUG-rank=", I5, " end post, tag=",I8, " pf%state%pstatus=", I2)',pf%rank,tag,pf%state%pstatus       
    end if

    if (ierror /= 0) then
       print *, 'Rank',pf%rank
       call pf_stop(__FILE__,__LINE__,'error during post',ierror)
   endif
  end subroutine pf_post

  !>  Subroutine to send this processor's convergence status to the next processor
  subroutine pf_send_status(pf, tag, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: direction
    integer ::  dir
    integer ::  istatus
    integer ::  ierror,dest

    if (pf%comm%nproc .eq. 1) return

    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction

    if (pf%rank == 0 .and. dir == 2) return
    if (pf%rank == pf%comm%nproc-1 .and. dir == 1) return

    ierror = 0
    istatus = pf%state%status
    if (dir == 1) then
       dest=pf%rank+1
    elseif (dir == 2) then
       dest=pf%rank-1
    else
      call pf_stop(__FILE__,__LINE__,'invalid dir during send_status',dir)
    end if

    if (pf%debug) print  '("DEBUG-rank=", I5, " begin send_status, tag=",I8, " pf%state%status=", I2)',pf%rank,tag,pf%state%status           
    call pf%comm%send_status(pf, tag, istatus, ierror, dest)

    if (ierror /= 0) then
       print *, 'Rank',pf%rank
       call pf_stop(__FILE__,__LINE__,'error during send_status',ierror)
    endif
    if (pf%debug) print  '("DEBUG-rank=", I5, " end send_status, tag=",I8, " pf%state%status=", I2)',pf%rank,tag,pf%state%status           
  end subroutine pf_send_status

  !>  Subroutine to receive the convergence status from the previous processor
  subroutine pf_recv_status(pf, tag, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: direction
    integer ::  dir
    integer ::  ierror, istatus,source

    if (pf%comm%nproc .eq. 1) return

    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction

    !  Return if this is the first processor
    if (pf%rank == 0 .and. dir == 1) return
    if (pf%rank == pf%comm%nproc-1 .and. dir == 2) return


    ierror = 0
    if (dir == 1) then
       source=pf%rank-1
    elseif (dir == 2) then
       source=pf%rank+1
    else
      print *, 'Rank',pf%rank
      call pf_stop(__FILE__,__LINE__,'invalid bad dir in recv_status',dir)
   end if
   if (pf%debug) print  '("DEBUG-rank=", I5, " begin recv_status, tag=",I8, " pf%state%pstatus=", I2)',pf%rank,tag,pf%state%pstatus   
   call pf%comm%recv_status(pf, tag, istatus, ierror, source)

   if (ierror .eq. 0) then
      pf%state%pstatus = istatus
   else
      print *, 'Rank=',pf%rank
      call pf_stop(__FILE__,__LINE__,'error during recv_status',ierror)      
    endif
    if (pf%debug) print  '("DEBUG-rank=", I5, " end recv_status, tag=",I8, " pf%state%pstatus=", I2,I2)',pf%rank,tag,pf%state%pstatus,istatus           
  end subroutine pf_recv_status
  !>  Subroutine to send the solution to the next processor
  subroutine pf_send(pf, level, tag, blocking, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    integer, optional, intent(in)    :: direction
    integer                          :: dir, ierror,dest

    if (pf%comm%nproc .eq. 1) return

    dir = 1 ! default: send forward
    if(present(direction)) dir = direction

    if(pf%rank == 0 .and. dir == 2) return
    if(pf%rank == pf%comm%nproc-1 .and. dir == 1) return

    ierror = 0
    ! need to wait here to make sure last non-blocking send is done
    if(blocking .eqv. .false.) then
       if (pf%debug) print  '("DEBUG-rank=", I5, " begin wait in send, level=",I4)',pf%rank,level%index
       if (pf%save_timings > 1) call pf_start_timer(pf, T_WAIT, level%index)
       !       call pf_mpi_wait(pf, level%index, ierror)
       call pf%comm%wait(pf, level%index, ierror)       
       if (pf%save_timings > 1) call pf_stop_timer(pf, T_WAIT, level%index)
       if (ierror /= 0) then
          print *, 'Rank=',pf%rank
          call pf_stop(__FILE__,__LINE__,'error during send (wait)',ierror)
       end if
       if (pf%debug) print  '("DEBUG-rank=", I5, " end wait in send, level=",I4)',pf%rank,level%index
    end if
    
    if (pf%save_timings > 1) call pf_start_timer(pf, T_PACK, level%index)
    if(dir == 2) then
       call level%q0%pack(level%send, 2)
       dest=pf%rank-1
    else
       dest=pf%rank+1
       if(present(direction)) then  !  This is for the imk sweeper where the presence of a flag matters
          call level%qend%pack(level%send, 1)
       else
          call level%qend%pack(level%send)
       end if
    end if
    if (pf%save_timings > 1) call pf_stop_timer(pf, T_PACK, level%index)
    ierror = 0
    if (pf%debug) print  '("DEBUG-rank=", I5, " begin send, level=",I4, " tag=", I8," blocking=", L3," state%status=",I3)',pf%rank,level%index,tag,blocking,pf%state%status
    if (pf%save_timings > 1) call pf_start_timer(pf, T_SEND, level%index)
    call pf%comm%send(pf, level, tag, blocking, ierror, dest)
    if (pf%save_timings > 1) call pf_stop_timer(pf, T_SEND,level%index)
    if (ierror /= 0) then
       print *, 'Rank=',pf%rank
       call pf_stop(__FILE__,__LINE__,'error during send',ierror)
    endif

    if (pf%debug) print  '("DEBUG-rank=", I5, " end send, level=",I4, " tag=", I8," blocking=", L3," state%status=",I3)',pf%rank,level%index,tag,blocking,pf%state%status
  end subroutine pf_send

  !>  Subroutine to recieve the solution from the previous processor
  subroutine pf_recv(pf, level, tag, blocking, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    integer, optional, intent(in)    :: direction
    integer                          :: dir, ierror,source

    if (pf%comm%nproc .eq. 1) return

    dir = 1 ! default: send forward
    if(present(direction)) dir = direction

    ierror = 0

    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                  .and. dir == 1) then
       if (pf%debug) print  '("DEBUG-rank=", I5, " begin recv, blocking=" ,L4, " tag=",I8, " pf%state%pstatus=", I2)',pf%rank,blocking,tag,pf%state%pstatus       
       source=pf%rank-1
       if (pf%save_timings > 1) then
          if (blocking)  then
             call pf_start_timer(pf, T_RECEIVE, level%index)
          else
             call pf_start_timer(pf, T_WAIT, level%index)
          end if
       end if

       call pf%comm%recv(pf, level,tag, blocking, ierror, source)
       if (ierror .ne. 0) call pf_stop(__FILE__,__LINE__,'error during receive, rank=',pf%rank)
       
       if (pf%save_timings > 1) then
          if (blocking)  then
             call pf_stop_timer(pf, T_RECEIVE, level%index)
          else
             call pf_stop_timer(pf, T_WAIT, level%index)
          end if
       end if

       
!       if (pf%debug) print*,  'DEBUG --',pf%rank, 'end recv', SIZE(level%recv)


       if (pf%save_timings > 1) call pf_start_timer(pf, T_UNPACK, level%index)
       if (present(direction)) then
          call level%q0%unpack(level%recv, 1)
       else
          call level%q0%unpack(level%recv)
       end if
       if (pf%save_timings > 1) call pf_stop_timer(pf, T_UNPACK, level%index)
       if (pf%debug) print  '("DEBUG-rank=", I5, " end recv, blocking=" ,L4, " tag=",I8, " pf%state%pstatus=", I2)',pf%rank,blocking,tag,pf%state%pstatus

    elseif (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING &
         .and. dir == 2) then
       if (pf%debug) print  '("DEBUG-rank=", I5, " begin recv, blocking=" ,L4, " tag=",I8, " pf%state%pstatus=", I2)',pf%rank,blocking,tag,pf%state%pstatus       
       source=pf%rank+1
       
       if (pf%save_timings > 1) call pf_start_timer(pf, T_RECEIVE, level%index)
       call pf%comm%recv(pf, level,tag, blocking, ierror, source)
       if (ierror .ne. 0) call pf_stop(__FILE__,__LINE__,'error during receive, rank=',pf%rank)       
       if (pf%save_timings > 1) call pf_stop_timer(pf, T_RECEIVE,level%index)
       if (pf%debug) print*,  'DEBUG --',pf%rank, SIZE(level%recv)
       !  Unpack solution
       if (pf%save_timings > 1) call pf_start_timer(pf, T_UNPACK, level%index)
       if (present(direction)) then
          call level%qend%unpack(level%recv, 2)
       else
          call level%qend%unpack(level%recv)
       end if
       if (pf%save_timings > 1) call pf_stop_timer(pf, T_UNPACK, level%index)
       if (pf%debug) print  '("DEBUG-rank=", I5, " end recv, blocking=" ,L4, " tag=",I8, " pf%state%pstatus=", I2)',pf%rank,blocking,tag,pf%state%pstatus
    end if

  end subroutine pf_recv


  !>  Subroutine to broadcast the initial condition to all processors
  subroutine pf_broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: nvar, root
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer :: ierror

    if (pf%comm%nproc .eq. 1) return

    if (pf%save_timings > 1) call pf_start_timer(pf, T_BROADCAST)
    if(pf%debug) print *,'DEBUG-rank=',pf%rank,' begin broadcast'
    call pf%comm%broadcast(pf, y, nvar, root, ierror)
    if (ierror /= 0) then
       print *, 'Rank',pf%rank
       call pf_stop(__FILE__,__LINE__,'error during broadcast',ierror)
    endif
    if (pf%save_timings > 1) call pf_stop_timer(pf, T_BROADCAST)
    if(pf%debug)print *,'DEBUG-rank=',pf%rank,' end broadcast'
  end subroutine pf_broadcast

  !> Save current solution and function value so that future corrections can be computed
  subroutine save(pf, lev, flags)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: lev  !!  Level to save on
    integer, optional, intent(in)   :: flags !!  which component to save (state/adjoint)
    integer :: m, p

    !  Save the data so we can interpolate correction later
    if(lev%index < pf%state%finest_level) then
       do m = 1, lev%nnodes
          call lev%pQ(m)%copy(lev%Q(m), flags)
          if (lev%Finterp) then
             do p = 1,SIZE(lev%F(1,:))
                call lev%pF(m,p)%copy(lev%F(m,p), flags)
             end do
          end if
       end do
    end if
    
  end subroutine save

end module pf_mod_comm
