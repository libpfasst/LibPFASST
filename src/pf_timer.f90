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
!!$  ! if you add more timers here, make sure to update the timer def in pf_dtype.f90
  character(len=14), parameter :: timer_names(16) = (/ &
       'total       ',  &        ! 1
       'predictor   ',  &        ! 2
       'iteration   ',  &        ! 3
       'hooks       ',  &        ! 4
       'step        ',  &        ! 5
       'residual    ',  &        ! 6
       'broadcast   ',  &        ! 7
       'interp      ',  &        ! 8
       'restrict    ',  &        ! 9
       'receive     ',  &        ! 10
       'send        ',  &        ! 11
       'wait        ',  &        ! 12
       'sweep       ',  &        ! 13
       'feval       ',  &        ! 14
       'fcomp       ',  &        ! 15
       'aux         '/)          ! 16

  integer, parameter :: &
       T_TOTAL       = 1,  &
       T_PREDICTOR   = 2,  &
       T_ITERATION   = 3,  &
       T_HOOKS       = 4,  &
       T_STEP        = 5,  &
       T_RESIDUAL    = 6,  &
       T_BROADCAST   = 7,  &
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
  subroutine pf_start_timer(pf, timer, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer
    integer, optional, intent(in)    :: level_index

    double precision ::  t_wall
    integer :: l   !  Either the first index or level_index
    l=1
    if (present(level_index)) l=level_index
    t_wall=MPI_Wtime()
    select case(timer)
    case(T_TOTAL)
       pf%pf_timers%t_total=t_wall
    case(T_PREDICTOR)
       pf%pf_timers%t_predictor= t_wall
    case(T_ITERATION)
       pf%pf_timers%t_iteration= t_wall
    case(T_STEP)
       pf%pf_timers%t_step= t_wall
    case(T_HOOKS)
       pf%pf_timers%t_hooks(l) = t_wall
    case(T_SWEEP)
       pf%pf_timers%t_sweeps(l) = t_wall
    case(T_FEVAL)
       pf%pf_timers%t_feval(l) = t_wall
    case(T_FCOMP)
       pf%pf_timers%t_fcomp(l) = t_wall
    case(T_AUX)
       pf%pf_timers%t_aux(l) = t_wall
    case(T_BROADCAST)
       pf%pf_timers%t_broadcast = t_wall
    case(T_RESIDUAL)
       pf%pf_timers%t_residual(l) = t_wall
    case(T_INTERPOLATE)
       pf%pf_timers%t_interpolate(l) = t_wall
    case(T_RESTRICT)
       pf%pf_timers%t_restrict(l) = t_wall
    case(T_WAIT)
       pf%pf_timers%t_wait(l) = t_wall
    case(T_SEND)
       pf%pf_timers%t_send(l) = t_wall
    case(T_RECEIVE)
       pf%pf_timers%t_receive(l) = t_wall
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',timer)
    end select

  end subroutine pf_start_timer
  
  !>  Subroutine to stop a timer
  subroutine pf_stop_timer(pf, timer,level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer
    integer, optional, intent(in)    :: level_index

    double precision ::  t_wall,delta_t
    integer :: l   !  Either the first index or level_index
    l=1
    if (present(level_index)) l=level_index

    t_wall=MPI_Wtime()
    select case(timer)
    case(T_TOTAL)
       t_prev=pf%pf_timers%t_total
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_total=pf%pf_runtimes%t_total + delta_t
    case(T_PREDICTOR)
       t_prev=pf%pf_runtimes%t_predictor
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_predictor=pf%pf_runtimes%t_predictor + delta_t
    case(T_ITERATION)
       t_prev=pf%pf_timers%t_iteration       
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_iteration=pf%pf_runtimes%t_iteration + delta_t 
    case(T_STEP)
       t_prev=pf%pf_timers%t_step       
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_step=pf%pf_runtimes%t_step + delta_t 
    case(T_HOOKS)
       t_prev=pf%pf_timers%t_hooks(l)       
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_hooks(l)=pf%pf_runtimes%t_hooks(l) + delta_t   
    case(T_SWEEP)
       t_prev=pf%pf_timers%t_sweeps(l)      
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_sweeps(l)=pf%pf_runtimes%t_sweeps(l) + delta_t     
    case(T_FEVAL)
       t_prev=pf%pf_timers%t_feval(l)      
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_feval(l)=pf%pf_runtimes%t_feval(l) + delta_t     
    case(T_FCOMP)
       t_prev=pf%pf_timers%t_fcomp(l)      
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_fcomp(l)=pf%pf_runtimes%t_fcomp(l) + delta_t     
    case(T_AUX)
       t_prev=pf%pf_timers%t_aux(l)      
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_aux(l)=pf%pf_runtimes%t_aux(l) + delta_t     
    case(T_BROADCAST)
       t_prev=pf%pf_timers%t_broadcast
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_broadcast=pf%pf_runtimes%t_broadcast + delta_t   
    case(T_RESIDUAL)
       t_prev=pf%pf_timers%t_residual(l)       
       delta_t = t_wall - t_prev       
       pf%pf_runtimes%t_residual(l)=pf%pf_runtimes%t_residual(l) + delta_t   
    case(T_INTERPOLATE)
       t_prev=pf%pf_timers%t_interpolate(l)       
       delta_t = t_wall - t_prev       
       pf%pf_runtimes%t_interpolate(l)=pf%pf_runtimes%t_interpolate(l) + delta_t   
    case(T_RESTRICT)
       t_prev=pf%pf_timers%t_restrict(l)       
       delta_t = t_wall - t_prev       
       pf%pf_runtimes%t_restrict(l)=pf%pf_runtimes%t_restrict(l) + delta_t    
    case(T_WAIT)
       t_prev=pf%pf_timers%t_wait(l)       
       delta_t = t_wall - t_prev       
       pf%pf_runtimes%t_wait(l)=pf%pf_runtimes%t_wait(l) + delta_t     
    case(T_SEND)
       t_prev=pf%pf_timers%t_send(l)       
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_send(l)=pf%pf_runtimes%t_send(l) + delta_t     
    case(T_RECEIVE)
       t_prev=pf%pf_timers%t_receive(l)       
       delta_t = t_wall - t_prev
       pf%pf_runtimes%t_receive(l)=pf%pf_runtimes%t_receive(l) + delta_t      
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',timer)
    end select

    if (pf%save_timings .eq. 3) then
       write(*, '("timer:",a16,", rank: ",i3,", step: ",i4, ", level: ", i3,' &
            // '", iter: ",i3, " Current t: ",f20.8, " Elapse_begin: ",f20.8, " Elapse_end ",f20.8, " Delta t: ",f20.8)') &
            timer_names(timer), pf%rank, &
            pf%state%step, l, pf%state%iter,  &
            t_wall,t_prev-pf%pf_timers%t_total,t_wall-pf%pf_timers%t_total,delta_t
    end if

  end subroutine pf_stop_timer

end module pf_mod_timer
