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
       PF_PRE_BLOCK      = 1, &
       PF_POST_BLOCK     = 2, &
       PF_PRE_PREDICTOR  = 3, &
       PF_POST_PREDICTOR = 4, &
       PF_PRE_ITERATION  = 5, &
       PF_POST_ITERATION = 6, &
       PF_PRE_SWEEP      = 7, &
       PF_POST_SWEEP     = 8, &
       PF_PRE_STEP       = 9, &
       PF_POST_STEP      = 10, &
       PF_MAX_HOOK       = 10


  integer, parameter :: &
       PF_HOOK_LOG_ONE  = 1, &
       PF_HOOK_LOG_ALL  = 7, &
       PF_HOOK_LOG_LAST = 10

  character(len=20), parameter :: hook_names(PF_HOOK_LOG_LAST) = (/ &
       'pre-block          ',  &
       'post-block         ',  &
       'pre-predictor      ',  &
       'post-predictor     ',  &
       'pre-iteration      ',  &
       'post-iteration     ',  &
       'pre-sweep          ',  &
       'post-sweep         ',  &
       'pre-step           ',  &
       'post-step          ' /)

contains

  ! Add a procedure to the hook on the given level
  subroutine pf_add_hook(pf, level, hook, proc)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level
    integer,           intent(in)    :: hook
    procedure(pf_hook_p)             :: proc

    pf%nhooks(level,hook) = pf%nhooks(level,hook) + 1
    pf%hooks(level,hook,pf%nhooks(level,hook))%proc => proc
  end subroutine pf_add_hook

  ! Call hooks associated with the hook and level
  subroutine call_hooks(pf, level, hook)
    use pf_mod_timer
    type(pf_pfasst_t), intent(inout), target :: pf
    integer,           intent(in)            :: level, hook

    integer :: i, l

    call start_timer(pf, THOOKS)

    pf%state%hook = hook

    if (level == -1) then
       do l = 1, pf%nlevels
          do i = 1, pf%nhooks(l,hook)
             call pf%hooks(l,hook,i)%proc(pf, pf%levels(l), pf%state, pf%levels(l)%ctx)
          end do
       end do
    else
       l = level
       do i = 1, pf%nhooks(l,hook)
          call pf%hooks(l,hook,i)%proc(pf, pf%levels(l), pf%state, pf%levels(l)%ctx)
       end do
    end if

    call end_timer(pf, THOOKS)
  end subroutine call_hooks

end module pf_mod_hooks
