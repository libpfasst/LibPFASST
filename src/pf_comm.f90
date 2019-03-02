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
    if (pf%debug) print*,'DEBUG --', pf%rank, 'is beginning pf_post, state%pstatus=', pf%state%pstatus, 'with tag =', tag
    ierror = 0
    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                  .and. dir == 1) then
       source=pf%rank-1
       call pf%comm%post(pf, level, tag, ierror, source)
    elseif (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING &
         .and. dir == 2) then
       source=pf%rank+1
       call pf%comm%post(pf, level, tag, ierror, source)
    end if

    
    if (ierror /= 0) then
      print *, pf%rank, 'warning: error during post', ierror
      stop "pf_parallel:pf_post"
   endif
   if (pf%debug) print*,'DEBUG --', pf%rank, 'is leaving pf_post, state%pstatus=', pf%state%pstatus, 'with tag =', tag     
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
       print *, pf%rank, 'warning: bad dir during send_status', dir
      stop "pf_parallel:pf_send_status"
    end if

    if (pf%debug) print*, 'DEBUG --',pf%rank, 'begins send_status with status', istatus, 'with tag =', tag 
    call pf%comm%send_status(pf, tag, istatus, ierror, dest)
    if (pf%debug) print*, 'DEBUG --',pf%rank, 'ends send_status'
    
    if (ierror /= 0) then
      print *, pf%rank, 'warning: error during send_status', ierror
      stop "pf_parallel:pf_send_status"
    endif
    
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
    
    if (pf%debug) print*, 'DEBUG --',pf%rank, 'begin recv_status with pstatus=',pf%state%pstatus, ' tag=',tag
    ierror = 0
    if (dir == 1) then
       source=pf%rank-1
    elseif (dir == 2) then
       source=pf%rank+1       
    else
      print *, pf%rank, 'warning: bad dir in recv_status', dir
      stop "pf_parallel_oc:pf_recv_status"
   end if
   
   if (pf%debug) print*, pf%rank,  'is receiving status backwards with tag ', tag  

   call pf%comm%recv_status(pf, tag, istatus, ierror, source)

   if (ierror .eq. 0) then
      pf%state%pstatus = istatus
   else
      print *, pf%rank, 'warning: error during recv_status', ierror
      stop "pf_parallel_oc:pf_recv_status"
    endif    
   
   if (pf%debug) print *, pf%rank, 'status recvd = ', istatus 
    if (pf%debug) print*,  'DEBUG --',pf%rank, 'end recv_statuswith pstatus=',pf%state%pstatus,'tag=',tag       
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

    ierror = 0
    call start_timer(pf, TSEND + level%index - 1)
    if (pf%debug) print*,  'DEBUG --',pf%rank, 'begin send, tag=',tag,blocking,' pf%state%status =',pf%state%status

    call pf%comm%send(pf, level, tag, blocking, ierror, dest)   
    if (ierror /= 0) then
      print *, pf%rank, 'warning: error during send', ierror
      stop "pf_parallel:pf_send"
   endif

   if (pf%debug) print*,  'DEBUG --',pf%rank, 'end send, tag=',tag,blocking                  
    call end_timer(pf, TSEND + level%index - 1)
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
    call start_timer(pf, TRECEIVE + level%index - 1)
    if (pf%debug) print*,  'DEBUG --',pf%rank, 'begin recv, tag=',tag,blocking, "pf%state%pstatus=",pf%state%pstatus
    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                  .and. dir == 1) then
       source=pf%rank-1
       call pf%comm%recv(pf, level,tag, blocking, ierror, source)
       
       if (ierror .eq. 0) then
          if (present(direction)) then
             call level%q0%unpack(level%recv, 1)
          else 
             call level%q0%unpack(level%recv)
          end if
        end if
    elseif (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                     .and. dir == 2) then
       source=pf%rank+1
       call pf%comm%recv(pf, level,tag, blocking, ierror, source)

       if (ierror .eq. 0) then 
         if (present(direction)) then
           call level%qend%unpack(level%recv, 2)
         else
           call level%qend%unpack(level%recv)
        end if
  
       end if
    end if

    if (pf%debug) print*,  'DEBUG --',pf%rank, 'end recv, tag=',tag,blocking    
    if(ierror .ne. 0) then
      print *, pf%rank, 'warning: mpi error during receive', ierror
      stop "pf_parallel:pf_recv"     
    end if
    
    call end_timer(pf, TRECEIVE + level%index - 1)
  end subroutine pf_recv


  !>  Subroutine to broadcast the initial condition to all processors
  subroutine pf_broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: nvar, root
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer :: ierror

    if (pf%comm%nproc .eq. 1) return
    
    call start_timer(pf, TBROADCAST)
    if(pf%debug) print *,'beginning broadcast'
    call pf%comm%broadcast(pf, y, nvar, root, ierror)
    if (ierror /= 0) then
       print *, pf%rank, 'warning:  error during broadcast', ierror
       stop "pf_parallel:pf_broadcast"
    endif
    call end_timer(pf, TBROADCAST)
    if(pf%debug)print *,'ending broadcast'
  end subroutine pf_broadcast

  !> Save current solution and function value so that future corrections can be computed
  subroutine save(lev, flags)
    class(pf_level_t), intent(inout) :: lev  !!  Level to save on
    integer, optional, intent(in)   :: flags !!  which component to save (state/adjoint)
    integer :: m, p
    

    if (lev%Finterp) then
       if (allocated(lev%pFflt)) then
          do m = 1, lev%nnodes
             do p = 1,size(lev%F(1,:))
                call lev%pF(m,p)%copy(lev%F(m,p), flags)
             end do
             call lev%pQ(m)%copy(lev%Q(m), flags)
          end do
       end if
    else
       if (allocated(lev%pQ)) then
          do m = 1, lev%nnodes
             call lev%pQ(m)%copy(lev%Q(m), flags)
          end do
       end if
    end if
  end subroutine save
  
end module pf_mod_comm
