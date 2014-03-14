!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
contains

  subroutine echo_error_hook(pf, level, state, levelctx)
    use solutions, only: exact
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: levelctx

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)

    qend => array1(level%qend)

    call exact(state%t0+state%dt, level%nvars, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(qend-yexact))

  end subroutine echo_error_hook

  subroutine echo_residual_hook(pf, level, state, levelctx)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: levelctx

    real(pfdp), pointer :: r(:)

    r => array1(level%R(level%nnodes-1))

    print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
         state%step+1, state%iter, level%level, maxval(abs(r))
  end subroutine echo_residual_hook

end module hooks
