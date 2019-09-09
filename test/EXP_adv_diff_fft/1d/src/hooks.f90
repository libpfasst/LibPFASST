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
  subroutine echo_error(pf, level_index)
    use pf_my_sweeper, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: yexact(pf%levels(level_index)%mpibuflen)
    real(pfdp), pointer :: y_end(:)
    real(pfdp) :: maxerr,err0

    y_end => get_array1d(pf%levels(level_index)%qend)

    !>  compute the exact solution
    call exact(pf%state%t0+pf%state%dt, yexact)
    !>  compute error
    maxerr = maxval(abs(y_end-yexact))

    y_end => get_array1d(pf%levels(level_index)%q0)

    !>  compute the exact solution
    call exact(pf%state%t0, yexact)
    !>  compute error
    err0 = maxval(abs(y_end-yexact))

    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7,es14.7," res: ",es14.7)', &
         pf%state%step+1, pf%state%iter,level_index, maxerr,err0,pf%levels(level_index)%residual
    call flush(6)
  end subroutine echo_error

end module hooks
