!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_imex
  use encap
  use output
  implicit none
contains

  subroutine echo_error_hook(pf, level, state, ctx)
    use solutions, only: exact
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    class(pf_context_t), intent(inout) :: ctx

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)

    qend => get_array(level%qend%q)

    call exact(state%t0+state%dt, level%nvars, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(qend-yexact))

  end subroutine echo_error_hook


  subroutine echo_error_hook(pf, level, state, ctx)
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    class(pf_context_t), intent(inout) :: ctx

    real(pfdp), pointer :: qend(:)
    qend => get_array(level%qend%q)

    call dump_numpy()

  end subroutine echo_error_hook

end module hooks
