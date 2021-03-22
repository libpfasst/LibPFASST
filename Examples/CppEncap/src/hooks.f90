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

    real(pfdp) :: yexact
    real(pfdp) :: maxerr
    real(pfdp) :: residual
    class(scalar_encap), pointer :: y_end

    !> Get the solution at the end of this step
    y_end => cast_as_scalar(pf%levels(level_index)%qend)

    !>  compute the exact solution
    call exact(pf%state%t0+pf%state%dt, yexact)
    !>  compute error
    maxerr = abs(y_end%getval() - yexact)
    residual=pf%levels(level_index)%residual
    
    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
         pf%state%step+1, pf%state%iter,level_index, maxerr,residual
    call flush(6)
  end subroutine echo_error


end module hooks
