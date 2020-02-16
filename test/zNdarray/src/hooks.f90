!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pfasst
  use pf_mod_zndarray
  use pf_mod_zutils
  implicit none
contains
 
  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    use probin, only: grid_size
    use pf_my_sweeper, only: my_sweeper_t, as_my_sweeper
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(my_sweeper_t), pointer :: sweeper
    real(pfdp) :: maxerr,t
    type(pf_zndarray_t) :: y_ex      !<  the initial condition
    
    call zndarray_build(y_ex, pf%levels(level_index)%lev_shape)
    sweeper => as_my_sweeper(pf%levels(level_index)%ulevel%sweeper)    

    
    !>  compute the exact solution
    t=pf%state%t0+pf%state%dt
    call exact(sweeper%fft_tool,t, y_ex)

    !>  compute error
    call y_ex%axpy(-1.0d0,pf%levels(level_index)%qend)

    !>  compute error
    maxerr = y_ex%norm()
    print '("error: time: ", f8.4," step: ",i3.3," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," sweep: ",i2.2," error: ",es14.7," resid: ",es14.7," deltaq0: ",es14.7)', &
         t,pf%state%step+1, pf%rank, pf%state%iter,level_index,pf%state%sweep,maxerr,pf%levels(level_index)%residual,pf%levels(level_index)%max_delta_q0

    
    call flush(6)

    call zndarray_destroy(y_ex)

  end subroutine echo_error
  
   subroutine set_error(pf, level_index)
    use probin, only: grid_size
    use pf_my_sweeper, only: my_sweeper_t, as_my_sweeper
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(my_sweeper_t), pointer :: sweeper
    real(pfdp) :: maxerr,t
    type(pf_zndarray_t) :: y_ex      !<  the initial condition
    
    call zndarray_build(y_ex, pf%levels(level_index)%lev_shape)
    sweeper => as_my_sweeper(pf%levels(level_index)%ulevel%sweeper)    

    !>  compute the exact solution
    t=pf%state%t0+pf%state%dt
    call exact(sweeper%fft_tool,t, y_ex)

    !>  compute error
    call y_ex%axpy(-1.0d0,pf%levels(level_index)%qend)

    !>  compute error
    maxerr = y_ex%norm()

    call zndarray_destroy(y_ex)

    
    call pf_set_error(pf,level_index,maxerr)
  end subroutine set_error

end module hooks
