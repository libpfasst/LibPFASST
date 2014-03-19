!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
contains

  subroutine echo_error_hook(pf, level, state, levelctx)
    use pf_mod_utils
    use solutions, only: exact
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: levelctx

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)
    real(pfdp) :: relres,res,relerr,err,t

    qend => array1(level%qend)
    t = state%t0+state%dt   
    call exact(t, level%nvars, yexact)

    err = maxval(abs(qend-yexact))
    relerr = maxval(abs(qend-yexact))/maxval(abs(yexact))

    call pf_residual(pf, level, state%dt)
    res= level%residual
    relres = res/maxval(abs(yexact))
    print '(" lev:",i5," step:",i5," t=",es10.3," iter:",i3," Err:",es13.6," RErr:",es13.6," Res:",es13.6," RRes:",es13.6)', &
               level%level,state%step+1, t,state%iter, err,relerr,res,relres 

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
