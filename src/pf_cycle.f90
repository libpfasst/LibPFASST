!
! Copyright (C) 2013 Matthew Emmett and Michael Minion.
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

! Cycle building routines.

module pf_mod_cycle
  use pf_mod_dtype
  implicit none
contains

  !
  ! Build cycle (list of PFASST operations).
  !
  subroutine pf_cycle_build(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: l, c

    type(pf_stage_t), pointer :: stages(:)

    select case(pf%ctype)
    case (SDC_CYCLE_FULL)

       stop 'SDC_CYCLE_FULL NOT IMPLEMENTED YET'

    case (SDC_CYCLE_V)

       allocate(pf%cycles%start(2*(pf%nlevels-1)-1))
       allocate(pf%cycles%pfasst(2*pf%nlevels-1))
       allocate(pf%cycles%end(1))

       ! start: transfer from coarsest to finest, sweeping on the way up
       stages => pf%cycles%start
       c = 1
       do l = 2, pf%nlevels-1
          stages(c)%type = SDC_CYCLE_INTERP
          stages(c)%F    = l
          stages(c)%G    = l-1
          c = c + 1

          stages(c)%type = SDC_CYCLE_SWEEP
          stages(c)%F    = l
          stages(c)%G    = -1
          c = c + 1
       end do

       stages(c)%type = SDC_CYCLE_INTERP
       stages(c)%F    = pf%nlevels
       stages(c)%G    = pf%nlevels-1

       ! pfasst: v-cycle from finest, but end in middle
       stages => pf%cycles%pfasst

       c = 1
       do l = pf%nlevels, 2, -1
          stages(c)%type = SDC_CYCLE_DOWN
          stages(c)%F    = l
          stages(c)%G    = l-1
          c = c + 1
       end do

       stages(c)%type = SDC_CYCLE_BOTTOM
       stages(c)%F    = 1
       stages(c)%G    = -1
       c = c + 1

       do l = 2, pf%nlevels
          stages(c)%type = SDC_CYCLE_UP
          stages(c)%F    = l
          stages(c)%G    = l-1
          c = c + 1
       end do

       ! end: sweep on finest
       stages => pf%cycles%end

       stages(1)%type = SDC_CYCLE_SWEEP
       stages(1)%F    = pf%nlevels
       stages(1)%G    = -1
    end select

  end subroutine pf_cycle_build

end module pf_mod_cycle
