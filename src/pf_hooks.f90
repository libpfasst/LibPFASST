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

module pf_mod_hooks
  use pf_mod_dtype
  implicit none

  integer, parameter :: &
       PF_POST_SWEEP     = 1, & ! after each sdc sweep
       PF_POST_ITERATION = 2, & ! after each pfasst iteration
       PF_POST_STEP      = 3, & ! after each time step
       PF_PRE_ITERATION  = 4, & ! before each pfasst iteration
       PF_POST_PREDICTOR = 5, &
       PF_PRE_BLOCK      = 6

contains

  ! Add a procedure to the hook on the given level
  subroutine pf_add_hook(pf, level, hook, proc)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level
    integer,           intent(in)    :: hook
    procedure(pf_hook_p)             :: proc

    pf%nhooks(level) = pf%nhooks(level) + 1
    pf%hooks(level,pf%nhooks(level))%proc => proc
    pf%hooks(level,pf%nhooks(level))%hook = hook
  end subroutine pf_add_hook

  ! Call hooks associated with the hook and level
  subroutine call_hooks(pf, level, hook)
    use pf_mod_timer
    type(pf_pfasst_t), intent(inout), target :: pf
    integer,           intent(in)            :: level, hook

    integer :: i, l

    call start_timer(pf, THOOKS)

    if (level == -1) then
       do l = 1, pf%nlevels
          do i = 1, pf%nhooks(l)
             if (pf%hooks(l,i)%hook == hook) then
                call pf%hooks(l,i)%proc(pf, pf%levels(l), pf%state, pf%levels(l)%ctx)
             end if
          end do
       end do
    else
       l = level
       do i = 1, pf%nhooks(l)
          if (pf%hooks(l,i)%hook == hook) then
             call pf%hooks(l,i)%proc(pf, pf%levels(l), pf%state, pf%levels(l)%ctx)
          end if
       end do
    end if

    call end_timer(pf, THOOKS)
  end subroutine call_hooks

end module pf_mod_hooks
