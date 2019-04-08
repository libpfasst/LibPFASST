!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_ndarray
  use pf_my_sweeper
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

    !> Get the solution at the end of this step    
    y_end => get_array1d(level%qend)

    !>  compute the exact solution
    call exact(state%t0+state%dt, yexact)

    !>  compute error
    maxerr = maxval(abs(y_end-yexact))
    
    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
         state%step+1, state%iter,level%index, maxerr,level%residual
    call flush(6)
  end subroutine echo_error

end module hooks
