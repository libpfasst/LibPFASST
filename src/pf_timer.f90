!!  Timing routines
!
! This file is part of LIBPFASST.
!
!> Module for setting timers
module pf_mod_timer
  use pf_mod_dtype
  use pf_mod_mpi, only: MPI_Wtime
  implicit none
  !  List of timers 
  integer, parameter :: &
       TTOTAL       = 1,  &
       TPREDICTOR   = 2,  &
       TITERATION   = 3,  &
       THOOKS       = 4,  &
       TSTEP        = 5,  &
       TRESIDUAL    = 6,  &
       TBROADCAST   = 7,  &
       TINTERPOLATE = 10,  &
       TRESTRICT    = 20,  &
       TRECEIVE     = 30,  &
       TSEND        = 40,  &
       TLEVEL       = 50,  &
       TAUX         = 60

  ! if you add more timers here, make sure to update the timer arrays in pf_dtype.f90

  character(len=14), parameter :: timer_names(62) = (/ &
       'total       ',  &        ! 1
       'predictor   ',  &
       'iteration   ',  &
       'hooks       ',  &
       'step        ',  &        ! 5
       'residual    ',  &
       'broadcast   ',  &
       '8           ',  &
       '9           ',  &
       'interpL1    ',  &    ! 10
       'interpL2    ',  &
       'interpL3    ',  &
       'interpL4    ',  &
       'interpL5    ',  &
       'interpL6    ',  &
       'interpL7    ',  &
       'interpL8    ',  &
       'interpL9    ',  &
       'interpL10   ',  &        
       'restrictL1  ',  &  ! 20
       'restrictL2  ',  &
       'restrictL3  ',  &
       'restrictL4  ',  &
       'restrictL5  ',  &
       'restrictL6  ',  &
       'restrictL7  ',  &
       'restrictL8  ',  &
       'restrictL9  ',  &
       'restrictL10 ',  &        
       'recvL1      ',  &      ! 30
       'recvL2      ',  &
       'recvL3      ',  &
       'recvL4      ',  &
       'recvL5      ',  &
       'recvL6      ',  &
       'recvL7      ',  &
       'recvL8      ',  &
       'recvL9      ',  &
       'recvL10     ',  &  
       'sendL1      ',  &      ! 40
       'sendL2      ',  &
       'sendL3      ',  &
       'sendL4      ',  &
       'sendL5      ',  &
       'sendL6      ',  &
       'sendL7      ',  &
       'sendL8      ',  &
       'sendL9      ',  &
       'sendL10     ',  &
       'sweepL1     ',  &       ! 50
       'sweepL2     ',  &
       'sweepL3     ',  &
       'sweepL4     ',  &
       'sweepL5     ',  &
       'sweepL6     ',  &
       'sweepL7     ',  &
       'sweepL8     ',  &
       'sweepL9     ',  &
       'sweepL10    ',  &     
       'exp         ', &        ! 60
       'omega       ', &
       'feval       '/)

contains
  !>  Subroutine to start a timer
  subroutine start_timer(pf, timer)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer

    pf%timers(timer)=MPI_Wtime()

  end subroutine start_timer

  !>  Subroutine to stop a timer
  subroutine end_timer(pf, timer)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer

    double precision :: t

    t=MPI_Wtime()

    pf%runtimes(timer) = pf%runtimes(timer) + t - pf%timers(timer)

    if (pf%echo_timings) then
       write(*, '("timer:",a16,", rank: ",i3,", step: ",i4, ", level: ", i3,' &
            // '", iter: ",i3, f23.8,f23.8,f23.8)') &
            timer_names(timer), pf%rank, &
            pf%state%step, pf%state%level_index, pf%state%iter,  &
            t-pf%timers(timer), pf%runtimes(timer), t-pf%timers(TTOTAL)
    end if

  end subroutine end_timer

end module pf_mod_timer
