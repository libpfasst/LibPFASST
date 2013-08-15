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
       
       ! start: transfer from coarsest to finest, sweeping on the way up
       if (pf%nlevels > 1) then
          allocate(pf%cycles%start(2*(pf%nlevels-1)-1))

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
       else

          allocate(pf%cycles%start(0))

       end if

       ! pfasst: v-cycle from middle
       allocate(pf%cycles%pfasst(2*pf%nlevels-1))
       stages => pf%cycles%pfasst

       c = 1
       do l = pf%nlevels-1, 2, -1
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
       ! XXX: get rid of this crap...
       allocate(pf%cycles%end(1))
       stages => pf%cycles%end

       stages(1)%type = SDC_CYCLE_SWEEP
       stages(1)%F    = pf%nlevels
       stages(1)%G    = -1

    end select

  end subroutine pf_cycle_build

  subroutine pf_cycle_print_stage(stage)
    type(pf_stage_t), intent(in) :: stage

    select case(stage%type)
    case (SDC_CYCLE_UP)
       print *, "==> cycle up:     F/G: ", stage%F, stage%G
    case (SDC_CYCLE_DOWN)
       print *, "==> cycle down:   F/G: ", stage%F, stage%G
    case (SDC_CYCLE_BOTTOM)
       print *, "==> cycle bottom: F:   ", stage%F
    case (SDC_CYCLE_SWEEP)
       print *, "==> cycle sweep:  F:   ", stage%F
    case (SDC_CYCLE_INTERP)
       print *, "==> cycle interp: F/G: ", stage%F, stage%G
    end select
  end subroutine pf_cycle_print_stage

  subroutine pf_cycle_print(pf)
    type(pf_pfasst_t), intent(in) :: pf
    
    integer :: c

    print *, 'STARTING (POST PREDICTOR) STAGES'
    if (pf%nlevels > 1) then
       do c = 1, size(pf%cycles%start)
          call pf_cycle_print_stage(pf%cycles%start(c))
       end do
    end if

    print *, 'PFASST (ITERATION) STAGES'
    do c = 1, size(pf%cycles%pfasst)
       call pf_cycle_print_stage(pf%cycles%pfasst(c))
    end do

    print *, 'END STAGES'
    do c = 1, size(pf%cycles%end)
       call pf_cycle_print_stage(pf%cycles%end(c))
    end do


  end subroutine pf_cycle_print


end module pf_mod_cycle
