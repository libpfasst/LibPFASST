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
    real(pfdp) :: residual
    real(pfdp) :: t
    
    !>  compute the error at last end point
    y_end => get_array1d(pf%levels(level_index)%qend)
    call exact(pf%state%t0+pf%state%dt, yexact)
    maxerr = maxval(abs(y_end-yexact))
    residual=pf%levels(level_index)%residual
    t = pf%state%t0+pf%state%dt

    !print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
    !     pf%state%step+1, pf%state%iter,level_index, maxerr,residual
    print '("error: time: ", f8.4," step: ",i8.1," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," resid: ",es14.7," deltaq0: ",es14.7)', &
         t,pf%state%step+1, pf%rank, pf%state%iter,level_index,maxerr,pf%levels(level_index)%residual,pf%levels(level_index)%max_delta_q0    

    call pf_set_error(pf,level_index,maxerr)
    
  end subroutine echo_error
end module hooks
