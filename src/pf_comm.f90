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
  subroutine pf_post(pf, level, tag, dir,stride)
    type(pf_pfasst_t), intent(in)    :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: dir      !  The direction to communicate
    integer, optional, intent(in)    :: stride   !  The distance to communicate

   !  Local variables
    integer :: dir_=1              ! default 1: send forward; set to -1 for send backwards
    integer :: stride_=1           ! default 1: send to next processor 
    integer :: ierr=0             ! error flag
    integer :: source             ! rank of sending processor
    integer :: rank, pstatus      ! local holders to be tidier

    pstatus=pf%state%pstatus   !  Tells if the previous processor is done iterating

    !  We don't need to recieve if the previous processor is done, or there is only one
    if (pf%comm%nproc .eq. 1 .or. pstatus /= PF_STATUS_ITERATING) return
    
    rank=pf%rank               !  Rank of recieving processor

    ! See what the direction is
    if(present(dir)) dir_ = dir
    if(present(stride)) stride_ = stride
    if(abs(dir_) .gt. 1) call pf_stop(__FILE__,__LINE__,'bad value for dir',val=dir_,rank=pf%rank)
    
    if (pf%debug) print  '("DEBUG-rank=", I5, " begin post, tag=",I8, " pstatus=", I2)',rank,tag,pstatus
    
    source=rank - (dir_)*(stride_)
    if (source < 0  .or. source > pf%comm%nproc-1) return 
    call pf%comm%post(pf, level, tag, ierr, source)
!    if (rank /= 0 .and.  dir == 1) then    
!       source=rank-1
!       call pf%comm%post(pf, level, tag, ierr, source)
!    elseif (rank /= pf%comm%nproc-1 .and. dir == 2) then       
!       source=rank+1
!       call pf%comm%post(pf, level, tag, ierr, source)
!    end if

    !  Check to see if there was an error
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'error during post, ierr',val=ierr,rank=rank)  
    if (pf%debug) print  '("DEBUG-rank=", I5, " end post, tag=",I8, " pstatus=", I2)',rank,tag,pstatus       
  end subroutine pf_post

  !>  Subroutine to send this processor's convergence status to the next processor
  subroutine pf_send_status(pf, tag, dir,stride)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: dir      !  The direction to communicate
    integer, optional, intent(in)    :: stride   !  The distance to communicate
    integer ::  dir_=1
    integer ::  stride_=1
    integer ::  ierr=0            ! error flag
    integer ::  dest              ! processor to send to
    integer ::  rank,istatus      ! local holders to be tidy

    if (pf%comm%nproc .eq. 1) return


    rank=pf%rank               !  Rank of sending processor
    istatus = pf%state%status  !  The status to send forward

    if(present(dir)) dir_ = dir
    if(present(stride)) stride_ = stride

    if(abs(dir_) .gt. 1) call pf_stop(__FILE__,__LINE__,'bad value for dir',val=dir_,rank=rank)

    dest=rank + (dir_)*(stride_)
    if (dest < 0  .or. dest > pf%comm%nproc-1) return 
    
    if (pf%debug) print  '("DEBUG-rank=",I5," begin send_status, tag=",I8, " istatus=",I2)',rank,tag,istatus
    
    call pf%comm%send_status(pf, tag, istatus, ierr, dest)

    if (ierr /= 0)  call pf_stop(__FILE__,__LINE__,'error during send_status, ierr',val=ierr,rank=rank)
    if (pf%debug) print  '("DEBUG-rank=",I5," end send_status, tag=",I8," istatus=",I2)',rank,tag,istatus
  end subroutine pf_send_status
  
  !>  Subroutine to receive the convergence status from the previous processor
  subroutine pf_recv_status(pf, tag, dir,stride)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: dir      !  The direction to communicate
    integer, optional, intent(in)    :: stride   !  The distance to communicate
    integer ::  dir_=1
    integer ::  stride_=1
    integer ::  ierr=0            ! error flag
    integer ::  source            ! processor to receive from
    integer ::  rank,pstatus
    
    if (pf%comm%nproc .eq. 1) return
    
    rank=pf%rank               !  Rank of sending processor
    pstatus=pf%state%pstatus 
    
    if(present(dir)) dir_ = dir
    if(present(stride)) stride_ = stride
    if(abs(dir_) .gt. 1) call pf_stop(__FILE__,__LINE__,'bad value for dir',val=dir_,rank=rank)
    

    source=pf%rank-(dir_)*(stride_)
    if (source < 0  .or. source > pf%comm%nproc-1) return 
    if (pf%debug) print  '("DEBUG-rank=",I5," begin recv_status, tag=",I8," pstatus=",I2)',rank,tag,pstatus   
    call pf%comm%recv_status(pf, tag, pstatus, ierr, source)

   if (ierr .ne. 0) call pf_stop(__FILE__,__LINE__,'error during recv_status, ierr',val=ierr,rank=pf%rank)
   pf%state%pstatus = pstatus

   if (pf%debug) print  '("DEBUG-rank=",I5," end recv_status, tag=",I8," pstatus=",I2)',rank,tag,pstatus   

 end subroutine pf_recv_status
 
  !>  Subroutine to send the solution to the next processor
  subroutine pf_send(pf, level, tag, blocking, dir,which,stride)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    integer, optional, intent(in)    :: dir
    integer, optional, intent(in)    :: which
    integer, optional, intent(in)    :: stride
    !  Local variables 
    integer                          :: dir_=1
    integer                          :: stride_=1
    integer                          :: which_=1
    integer                          :: dest
    integer                          :: ierr=0     ! error flag
    integer ::  rank,istatus

    
    if (pf%comm%nproc .eq. 1) return
    
    rank=pf%rank               !  Rank of sending processor
    istatus=pf%state%status 

    if(present(dir)) dir_ = dir
    if(present(which)) which_ = which
    if(present(stride)) stride_ = stride
    if(abs(dir_) .gt. 1) call pf_stop(__FILE__,__LINE__,'bad value for dir',val=dir_,rank=rank)

    ! need to wait here to make sure last non-blocking send is done
    if(blocking .eqv. .false.) then

       if (pf%debug) print  '("DEBUG-rank=", I5, " begin wait, level=",I4)',rank,level%index

       call pf_start_timer(pf, T_WAIT, level%index)
       call pf%comm%wait(pf, level%index, ierr)       
       call pf_stop_timer(pf, T_WAIT, level%index)
       
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'error during wait, ierr',val=ierr,rank=rank)
       if (pf%debug) print  '("DEBUG-rank=", I5, " end wait, level=",I4)',rank,level%index
    end if

    !  Compute the destination
    dest=rank+(dir_)*(stride_)
    if (dest < 0  .or. dest > pf%comm%nproc-1) return

    !  Pack the solution to send
    call pf_start_timer(pf, T_PACK, level%index)
    if(dir_ == -1) then
       call level%q0%pack(level%send, which_)
    else
       call level%qend%pack(level%send, which_) ! The imk encap now ignores the flag in pack and upack
    end if
    call pf_stop_timer(pf, T_PACK, level%index)

    !  Do the send 
    if (pf%debug) print  '("DEBUG-rank=", I5, " begin send, level=",I4, " tag=", I8," blocking=", L3," status=",I3)',rank,level%index,tag,blocking,istatus

    call pf_start_timer(pf, T_SEND, level%index)
    call pf%comm%send(pf, level, tag, blocking, ierr, dest)
    call pf_stop_timer(pf, T_SEND,level%index)

    !  Check for an error
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'error during send, ierr',val=ierr,rank=rank)
    if (pf%debug) print  '("DEBUG-rank=", I5, " end send, level=",I4, " tag=", I8," blocking=", L3," status=",I3)',rank,level%index,tag,blocking,istatus

  end subroutine pf_send

  !>  Subroutine to recieve the solution from the previous processor
  subroutine pf_recv(pf, level, tag, blocking, dir,which,stride)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer, optional, intent(in)    :: dir
    integer, optional, intent(in)    :: which
    integer, optional, intent(in)    :: stride
    !  Local variables 
    integer                          :: dir_=1
    integer                          :: stride_=1
    integer                          :: which_=1

    integer                          :: source
    integer                          :: ierr=0      ! error flag
    integer                          :: ntimer
    integer ::  rank,pstatus

    pstatus=pf%state%pstatus
    if (pf%comm%nproc .eq. 1 .or. pstatus /= PF_STATUS_ITERATING) return
    rank=pf%rank
    
    if(present(dir)) dir_ = dir
    if(present(which)) which_ = which
    if(present(stride)) stride_ = stride
    if(abs(dir_) .gt. 1) call pf_stop(__FILE__,__LINE__,'bad value for dir',val=dir_,rank=rank)


    if (pf%debug) print '("DEBUG-rank=",I5," begin recv, blocking=",L4," tag=",I8," pstatus=", I2)',rank,blocking,tag,pstatus

    source=pf%rank-(dir_)*(stride_)
    if (source < 0  .or. source > pf%comm%nproc-1) return 

    if (blocking)  then
       ntimer=T_RECEIVE
    else
       ntimer=T_WAIT

    end if
    
    call pf_start_timer(pf, ntimer, level%index)
    call pf%comm%recv(pf, level,tag, blocking, ierr, source)
    call pf_stop_timer(pf, ntimer, level%index)
    
    if (ierr .ne. 0) call pf_stop(__FILE__,__LINE__,'error during receive, ierr',val=ierr,rank=pf%rank)
    if (pf%debug) print*,  'DEBUG --',pf%rank, 'end recv', SIZE(level%recv)
    
    !  Unpack solution
    call pf_start_timer(pf, T_UNPACK, level%index)
    if (dir_ == 1) then
       call level%q0%unpack(level%recv, which_)
    else
       call level%qend%unpack(level%recv, which_)
    end if
    call pf_stop_timer(pf, T_UNPACK, level%index)
  
    if (pf%debug) print  '("DEBUG-rank=", I5, " end recv, blocking=" ,L4, " tag=",I8, " pstatus=", I2)',rank,blocking,tag,pstatus
  end subroutine pf_recv


  !>  Subroutine to broadcast the initial condition to all processors
  subroutine pf_broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: nvar, root
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer :: ierr=0      ! error flag

    if (pf%comm%nproc .eq. 1) return

    if(pf%debug) print *,'DEBUG-rank=',pf%rank,' begin broadcast'

    call pf_start_timer(pf, T_BROADCAST)
    call pf%comm%broadcast(pf, y, nvar, root, ierr)
    call pf_stop_timer(pf, T_BROADCAST)

    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'error during broadcast, ierr',val=ierr,rank=pf%rank)

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
