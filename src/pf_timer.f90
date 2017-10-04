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

module pf_mod_timer
  use pf_mod_dtype
  implicit none

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
       'interp0     ',  &        ! 10
       'interp1     ',  &
       'interp2     ',  &
       'interp3     ',  &
       'interp4     ',  &
       'interp5     ',  &
       'interp6     ',  &
       'interp7     ',  &
       'interp8     ',  &
       'interp9     ',  &
       'restrict0   ',  &        ! 20
       'restrict1   ',  &
       'restrict2   ',  &
       'restrict3   ',  &
       'restrict4   ',  &
       'restrict5   ',  &
       'restrict6   ',  &
       'restrict7   ',  &
       'restrict8   ',  &
       'restrict9   ',  &
       'recv0       ',  &        ! 30
       'recv1       ',  &
       'recv2       ',  &
       'recv3       ',  &
       'recv4       ',  &
       'recv5       ',  &
       'recv6       ',  &
       'recv7       ',  &
       'recv8       ',  &
       'recv9       ',  &
       'send0       ',  &        ! 40
       'send1       ',  &
       'send2       ',  &
       'send3       ',  &
       'send4       ',  &
       'send5       ',  &
       'send6       ',  &
       'send7       ',  &
       'send8       ',  &
       'send9       ',  &
       'sweep0      ',  &        ! 50
       'sweep1      ',  &
       'sweep2      ',  &
       'sweep3      ',  &
       'sweep4      ',  &
       'sweep5      ',  &
       'sweep6      ',  &
       'sweep7      ',  &
       'sweep8      ',  &
       'sweep9      ',  &
       'exponentiate', &        ! 60
       'omega       ',  &
       'feval       '/)

contains

  subroutine start_timer(pf, timer)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer

    call system_clock(pf%timers(timer))
  end subroutine start_timer

  subroutine end_timer(pf, timer)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: timer

    integer(8) :: t, rate

    call system_clock(t, rate)
    pf%runtimes(timer) = pf%runtimes(timer) + t - pf%timers(timer)

    if (pf%echo_timings) then
       write(601+pf%rank, '("timer:",a16,", rank: ",i3,", step: ",i4,' &
       !write(*, '("timer:",a16,", rank: ",i3,", step: ",i3,' &
            // '", iter: ",i3,", cycle: ",i3,", time (rate ",i12,"Hz): ",i18,i18,i18,i18)') &
            timer_names(timer), pf%rank, &
            pf%state%step, pf%state%iter, pf%state%cycle, rate, &
            t-pf%timers(timer), t-pf%timers(TTOTAL), pf%timers(timer), t
       call flush(601+pf%rank)
    end if

  end subroutine end_timer

end module pf_mod_timer
