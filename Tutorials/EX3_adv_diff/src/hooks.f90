!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pfasst
  use pf_mod_ndarray
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    use my_sweeper, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: yexact(pf%levels(level_index)%lev_shape(1))
    real(pfdp), pointer :: y_end(:)
    real(pfdp) :: maxerr
    real(pfdp) ::   time,resid
    integer ::   step,rank,iter
    time=pf%state%t0+pf%state%dt
    step=pf%state%step+1
    rank=pf%rank
    iter=pf%state%iter
    resid=pf%levels(level_index)%residual

    !>  compute the error at last end point
    y_end => get_array1d(pf%levels(level_index)%qend)
    call exact(pf%state%t0+pf%state%dt, yexact)
    maxerr = maxval(abs(y_end-yexact))
    
    print '("time: ", f10.4," step: ", i7.7," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
         time, step, rank, iter,level_index, maxerr,resid

    call pf_set_error(pf,level_index,maxerr)
    
  end subroutine echo_error
end module hooks
