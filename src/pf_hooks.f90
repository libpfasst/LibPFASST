!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module pf_mod_hooks
  use pf_mod_dtype
  implicit none

  integer, parameter :: &
       PF_POST_SWEEP = 1, &     ! call hooks after each sdc sweep
       PF_POST_ITERATION = 2, & ! call hooks after each pfasst iteration
       PF_POST_STEP = 3, &      ! call hooks after each time step
       PF_PRE_ITERATION = 4, &  ! call hooks before each pfasst iteration
       PF_POST_PREDICTOR = 5, &    
       PF_PRE_BLOCK = 6

  interface add_hook
     module procedure pf_add_hook
  end interface add_hook

contains

  ! Add a procedure to the hook on the given level
  subroutine pf_add_hook(pf, level, hook, proc)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level
    integer,           intent(in)    :: hook
    procedure(hook_proc)             :: proc

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
