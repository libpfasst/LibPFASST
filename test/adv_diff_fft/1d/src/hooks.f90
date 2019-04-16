!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level, state)
    use pf_my_sweeper, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(pfdp) :: yexact(level%mpibuflen)
    real(pfdp), pointer :: y_end(:)
    real(pfdp) :: maxerr

    !>  compute the error at last end point
    y_end => get_array1d(level%qend)
    call exact(state%t0+state%dt, yexact)
    maxerr = maxval(abs(y_end-yexact))

    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
         state%step+1, state%iter,level%index, maxerr,level%residual
    call flush(6)
  end subroutine echo_error
end module hooks
