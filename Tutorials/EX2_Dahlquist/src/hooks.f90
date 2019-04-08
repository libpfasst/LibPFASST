!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use encap
  use pf_my_sweeper
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level, state)
    use pf_my_sweeper, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(pfdp) :: yexact
    real(pfdp) :: maxerr
    class(scalar_encap), pointer :: y_end

    !> Get the solution at the end of this step
    y_end => cast_as_scalar(level%qend)

    !>  compute the exact solution
    call exact(state%t0+state%dt, yexact)
    !>  compute error
    maxerr = abs(y_end%y-yexact)
    
    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
         state%step+1, state%iter,level%index, maxerr,level%residual
    call flush(6)
  end subroutine echo_error


end module hooks
