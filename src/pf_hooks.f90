!! Module for calling of user defined routines from various places in the pfasst algorithm
!
! This file is part of LIBPFASST.
!
!>  Module for the calling of user defined routines from various places in the pfasst algorithm
module pf_mod_hooks
  use pf_mod_dtype
  implicit none

  !>  Define hook indices
  integer, parameter :: &
       PF_PRE_PREDICTOR     = 1, &
       PF_POST_PREDICTOR    = 2, &
       PF_PRE_ITERATION     = 3, &
       PF_POST_ITERATION    = 4, &
       PF_PRE_SWEEP         = 5, &
       PF_POST_SWEEP        = 6, &
       PF_PRE_STEP          = 7, &
       PF_POST_STEP         = 8, &
       PF_PRE_INTERP_ALL    = 9, &
       PF_POST_INTERP_ALL   = 10, &
       PF_PRE_INTERP_Q0     = 11, &
       PF_POST_INTERP_Q0    = 12, &
       PF_PRE_RESTRICT_ALL  = 13, &
       PF_POST_RESTRICT_ALL = 14, &
       PF_PRE_CONVERGENCE   = 15, &
       PF_POST_CONVERGENCE  = 16, &
       PF_MAX_HOOK          = 16

  integer, parameter :: &
       PF_HOOK_LOG_ONE  = 1, &
       PF_HOOK_LOG_ALL  = 7, &
       PF_HOOK_LOG_LAST = PF_MAX_HOOK

  !>  Define hook names
  character(len=20), parameter :: hook_names(PF_HOOK_LOG_LAST) = (/ &
       'pre-predictor      ',  &
       'post-predictor     ',  &
       'pre-iteration      ',  &
       'post-iteration     ',  &
       'pre-sweep          ',  &
       'post-sweep         ',  &
       'pre-step           ',  &
       'post-step          ',  &
       'pre-interp-all     ',  &
       'post-interp-all    ',  &
       'pre-interp-q0      ',  &
       'post-interp-q0     ',  &
       'pre-restrict-all   ',  &
       'post-restrict-all  ',  &
       'pre-convergence    ',  &
       'post-convergence   ' /)

contains

  !> Subroutine to add a procedure to the hook on the given level
  subroutine pf_add_hook(pf, level_index, hook, proc)
    type(pf_pfasst_t), intent(inout) :: pf            !! main pfasst structure
    integer,           intent(in)    :: level_index     !! which pfasst level to add hook
    integer,           intent(in)    :: hook          !! which hook to add
    procedure(pf_hook_p)             :: proc          !! precudre to call from hook

    integer :: l   !

    if (level_index == -1) then  ! Do to all levels
       do l = 1, pf%nlevels
          pf%nhooks(l,hook) = pf%nhooks(l,hook) + 1
          pf%hooks(l,hook,pf%nhooks(l,hook))%proc => proc
       end do
    else  ! Do to just level level_index
       pf%nhooks(level_index,hook) = pf%nhooks(level_index,hook) + 1
       pf%hooks(level_index,hook,pf%nhooks(level_index,hook))%proc => proc
    end if

  end subroutine pf_add_hook

  !> Subroutine to call hooks associated with the hook and level
  subroutine call_hooks(pf, level_index, hook)
    use pf_mod_timer
    type(pf_pfasst_t), intent(inout), target :: pf         !! main pfasst structure
    integer,           intent(in)            :: level_index  !! which pfasst level to call hook
    integer,           intent(in)            :: hook       !! which hook to call

    integer :: i  !!  hook loop index
    integer :: l  !!  level loop index

    call start_timer(pf, THOOKS)

    pf%state%hook = hook
    if (level_index == -1) then  ! Do to all levels
       do l = 1, pf%nlevels
          do i = 1, pf%nhooks(l,hook)
             call pf%hooks(l,hook,i)%proc(pf,l)
          end do
       end do
    else  ! Do to just level level_index
       do i = 1, pf%nhooks(level_index,hook)
          call pf%hooks(level_index,hook,i)%proc(pf,level_index)
       end do
    end if

    call end_timer(pf, THOOKS)
  end subroutine call_hooks

  !>  Subroutine defining log hook
  subroutine pf_logger_hook(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index
    
    print '("PF:: trank: ",i4,", step: ",i6,", iter: ",i3,", level: ",i2," location: ",a)', &
         pf%rank, pf%state%step, pf%state%iter, level_index, hook_names(pf%state%hook)
  end subroutine pf_logger_hook

  !>  Subroutine to add log hook
  subroutine pf_logger_attach(pf)
    type(pf_pfasst_t), intent(inout) :: pf
    
    integer :: l, h
    
    do h = PF_HOOK_LOG_ONE, PF_HOOK_LOG_ALL-1
       call pf_add_hook(pf, 1, h, pf_logger_hook)
    end do
    
    do l = 1, pf%nlevels
       do h = PF_HOOK_LOG_ALL, PF_HOOK_LOG_LAST
          call pf_add_hook(pf, l, h, pf_logger_hook)
       end do
    end do
  end subroutine pf_logger_attach

end module pf_mod_hooks
