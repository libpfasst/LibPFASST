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
       'total     ',  &   ! 1:  Total time to do all blocks
       'predictor ',  &   ! 2:  Time in predictor   
       'block     ',  &   ! 3:  Time for all blocks
       'iteration ',  &   ! 4:  Time for all iterations 
       'sweep     ',  &   ! 5:  Time for all sweeps (n-steps in parareal)
       'feval     ',  &   ! 6:  Time for explicit function evaluations
       'fcomp     ',  &   ! 7:  Time for implicit function evaluations
       'residual  ',  &   ! 8:  Time for computing residuals
       'interp    ',  &   ! 9:  Interpolation time
       'restrict  ',  &   ! 10: Restricting time
       'broadcast ',  &   ! 11: Time for broadcast (of initial conditions)
       'receive   ',  &   ! 12: Time in receive in pf_comm
       'send      ',  &   ! 13: Time in send in pf_comm
       'wait      ',  &   ! 14: Time in wait in 
       'pack      ',  &   ! 15: Time to pack solution
       'unpack    ',  &   ! 16: Time to unpack solution
       'hooks     ',  &   ! 17: Time in hooks routines
       'aux       '/)     ! 18: Extra for whatever

  ! Assign numbers to timer names
  integer, parameter :: &
       T_TOTAL       = 1,  &
       T_PREDICTOR   = 2,  &
       T_BLOCK       = 3,  &
       T_ITERATION   = 4,  &
       T_SWEEP       = 5,  &  
       T_FEVAL       = 6,  &  
       T_FCOMP       = 7,  &  
       T_RESIDUAL    = 8,  &
       T_INTERPOLATE = 9,  &
       T_RESTRICT    = 10, &
       T_BROADCAST   = 11, &
       T_RECEIVE     = 12, &
       T_SEND        = 13, &
       T_WAIT        = 14, &
       T_PACK        = 15, &
       T_UNPACK      = 16, &
       T_HOOKS       = 17, &
       T_AUX         = 18    
  
contains
  !>  Subroutine to start a timer
  subroutine pf_start_timer(pf, timer_index, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer_index
    integer, optional, intent(in)    :: level_index

    double precision ::  t_wall
    integer :: l   !  Either the first index or level_index
    if (pf%save_timings .eq. 0) return
    if ((pf%save_timings .eq. 1) .and. (timer_index .ne. T_TOTAL)) return    
    
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
