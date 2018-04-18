!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
contains

  subroutine echo_error(pf, level, state)
    use iso_c_binding
    use feval, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(c_double) :: yexact(level%mpibuflen)
    real(pfdp), pointer :: qend(:),r(:)

    qend => array1(level%qend)
    r => array1(level%R(level%nnodes-1))

    call exact(state%t0+state%dt, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
         state%step+1, state%iter,level%index, maxval(abs(qend-yexact)),maxval(abs(r))
    call flush
  end subroutine echo_error


  subroutine echo_residual(pf, level, state)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(pfdp), pointer :: r(:)

    r => array1(level%R(level%nnodes-1))

    print '("resid: time: ", f8.4," rank: ",i3.3," step: ",i5.5," iter: ",i3.3," level: ",i1.1," resid: ",es18.10e4)', &
         state%t0+state%dt, pf%rank, state%step+1, state%iter, level%index, maxval(abs(r))
    call flush
  end subroutine echo_residual

end module hooks
