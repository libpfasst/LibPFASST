!!  Timing routines
!
! This file is part of LIBPFASST.
!
!> Module for setting timers
module pf_mod_timer
  use pf_mod_dtype
  use pf_mod_stop
  use pf_mod_mpi, only: MPI_Wtime
  implicit none
  !!$  ! if you add more timers here, make sure to update the PF_NUM_TIMERS in pf_dtype.f90
  character(len=10), parameter :: timer_names(PF_NUM_TIMERS) = (/ &
       'total     ',  &        ! 1
       'predictor ',  &        ! 2
       'iteration ',  &        ! 3
       'step      ',  &        ! 4
       'broadcast ',  &        ! 5
       'hooks     ',  &        ! 6
       'residual  ',  &        ! 7
       'interp    ',  &        ! 8
       'restrict  ',  &        ! 9
       'receive   ',  &        ! 10
       'send      ',  &        ! 11
       'wait      ',  &        ! 12
       'sweep     ',  &        ! 13
       'feval     ',  &        ! 14
       'fcomp     ',  &        ! 15
       'aux       '/)          ! 16

  ! Assign numbers to timer names
  integer, parameter :: &
       T_TOTAL       = 1,  &
       T_PREDICTOR   = 2,  &
       T_ITERATION   = 3,  &
       T_STEP        = 4,  &
       T_BROADCAST   = 5,  &
       T_HOOKS       = 6,  &
       T_RESIDUAL    = 7,  &
       T_INTERPOLATE = 8,  &
       T_RESTRICT    = 9,  &
       T_RECEIVE     = 10,  &
       T_SEND        = 11,  &
       T_WAIT        = 12,  &
       T_SWEEP       = 13,  &  
       T_FEVAL       = 14,  &  
       T_FCOMP       = 15,  &  
       T_AUX         = 16    
  
contains
  !>  Subroutine to start a timer
  subroutine pf_start_timer(pf, timer_index, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer_index
    integer, optional, intent(in)    :: level_index

    double precision ::  t_wall
    integer :: l   !  Either the first index or level_index
    l=1
    if (present(level_index)) l=level_index

    pf%pf_timers%timers(timer_index,l) = MPI_Wtime()
    if (pf%save_timings .eq. 3) then
       write(*, '("start timer:",a10,", rank:",i3,", step:",i4,", level:",i1,", iter: ",i3, ' &
            // '" Current t: ",f20.8, " Elapsed_time: ",f20.8)') &
            timer_names(timer_index), pf%rank, pf%state%step, l, pf%state%iter,  &
            pf%pf_timers%timers(timer_index,l),pf%pf_timers%timers(timer_index,l)-pf%pf_timers%timers(T_TOTAL,1)
    end if

  end subroutine pf_start_timer
  
  !>  Subroutine to stop a timer
  subroutine pf_stop_timer(pf, timer_index,level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer_index
    integer, optional, intent(in)    :: level_index

    double precision ::  t_wall  !  Current wall clock
    double precision ::  t_prev  !  Time that timer was started (relative to T_TOTAL)
    double precision ::  t_now   !  Current time relative to T_TOTAL)
    double precision ::  delta_t !  Elapsed time for this time
    integer :: l   !  Either the first index or level_index
    l=1
    if (present(level_index)) l=level_index

    t_wall=MPI_Wtime()
    t_prev=pf%pf_timers%timers(timer_index,l)-pf%pf_timers%timers(T_TOTAL,1)
    t_now=t_wall-pf%pf_timers%timers(T_TOTAL,1)
    delta_t = t_now - t_prev
    
    pf%pf_timers%timers(timer_index,l)=t_wall
    pf%pf_timers%runtimes(timer_index,l)= pf%pf_timers%runtimes(timer_index,l) + delta_t

    !  Echo timings 
    if (pf%save_timings .eq. 3) then
       write(*, '("stop timer:",a10,", rank:",i3,", step:",i4,", level:",i1,", iter: ",i3, ' &
            // '" Wall t: ",f20.8, " begin t: ",f20.8, " end t: ",f20.8, " Delta t: ",f20.8, " Cum: ",f20.8)') &
            timer_names(timer_index), pf%rank, pf%state%step, l, pf%state%iter,  &
            t_wall,t_prev,t_now,delta_t,pf%pf_timers%runtimes(timer_index,l)            
    end if


  end subroutine pf_stop_timer

end module pf_mod_timer
