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
  subroutine echo_error(pf, level_index)
    use pf_my_sweeper, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(scalar_encap), pointer :: y_end
    real(pfdp) :: yexact
    real(pfdp) :: maxerr
    real(pfdp) :: time,resid
    integer ::   step,rank,iter
    time=pf%state%t0+pf%state%dt
    step=pf%state%step+1
    rank=pf%rank
    iter=pf%state%iter

    !> Get the solution at the end of this step
    y_end => cast_as_scalar(pf%levels(level_index)%qend)

    !>  compute the exact solution
    call exact(time, yexact)
    !>  compute error
    maxerr = abs(y_end%y-yexact)
    resid=pf%levels(level_index)%residual
    
    print '("time: ", f10.4," step: ", i7.7," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
         time, step, rank, iter,level_index, maxerr,resid
    call flush(6)
  end subroutine echo_error


end module hooks
