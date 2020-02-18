!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_rutils
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    use probin, only: grid_size
    use pf_my_stepper, only: my_stepper_t, as_my_stepper
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(my_stepper_t), pointer :: stepper
    real(pfdp) :: maxerr,t
    type(pf_ndarray_t) :: y_ex      !<  the initial condition
    
    call ndarray_build(y_ex, pf%levels(level_index)%lev_shape)
    stepper => as_my_stepper(pf%levels(level_index)%ulevel%stepper)    

    
    !>  compute the exact solution
    t=pf%state%t0+pf%state%dt
    call exact(stepper%fft_tool,t, y_ex)

    !>  compute error
    call y_ex%axpy(-1.0d0,pf%levels(level_index)%qend)

    !>  compute error
    maxerr = y_ex%norm()
    print '("error: time: ", f8.4," step: ",i3.3," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," resid: ",es14.7," deltaq0: ",es14.7)', &
         t,pf%state%step+1, pf%rank, pf%state%iter,level_index,maxerr,pf%levels(level_index)%residual,pf%levels(level_index)%max_delta_q0


    call flush(6)

    call ndarray_destroy(y_ex)

    
    
  end subroutine echo_error

end module hooks
