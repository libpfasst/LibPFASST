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
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: y_exact
    real(pfdp) :: maxerr
    real(pfdp) :: residual
    real(pfdp) :: t
    real(pfdp), pointer :: y_end(:)

    !> Get the solution at the end of this step
    y_end => get_array1d(cast_as_ndarray(pf%levels(level_index)%qend))

    !>  compute the exact solution
    t = pf%state%t0+pf%state%dt
    call exact(t, y_exact)
    !>  compute error
    maxerr = abs(y_end(1) - y_exact)
    residual = pf%levels(level_index)%residual
    
    !print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
    !     pf%state%step+1, pf%state%iter,level_index, maxerr,residual
    print '("error: time: ", f8.4," step: ",i8.1," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," resid: ",es14.7," deltaq0: ",es14.7)', &
         t,pf%state%step+1, pf%rank, pf%state%iter,level_index,maxerr,pf%levels(level_index)%residual,pf%levels(level_index)%max_delta_q0
    call flush(6)
    
    call pf_set_error(pf,level_index,maxerr)

  end subroutine echo_error

  subroutine exact(t, yex)
    use probin, only: lam1,lam2
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex

    yex = exp((lam1 + lam2) * t)
  end subroutine exact

end module hooks
