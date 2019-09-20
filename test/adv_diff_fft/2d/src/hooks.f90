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

    real(pfdp), pointer :: y_end(:,:),y_ex(:,:)
    real(pfdp) :: maxerr
    type(ndarray), target :: y_exact      !<  the initial condition
    y_end => get_array2d(pf%levels(level_index)%qend)
    call ndarray_build(y_exact, [ pf%levels(level_index)%lev_shape])
    y_ex => get_array2d(y_exact)    
    
    !>  compute the exact solution
    call exact(pf%state%t0+pf%state%dt, y_ex)
    !>  compute error
    maxerr = maxval(abs(y_end-y_ex))
    
    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
         pf%state%step+1, pf%state%iter,level_index, maxerr,pf%levels(level_index)%residual
    call flush(6)
    call ndarray_destroy(y_exact)    
  end subroutine echo_error

end module hooks
