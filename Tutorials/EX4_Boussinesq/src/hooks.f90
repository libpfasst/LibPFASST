!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
  use pf_mod_zndsysarray
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    use pf_my_sweeper, only: exact,my_sweeper_t, as_my_sweeper
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(my_sweeper_t), pointer :: sweeper
    complex(pfdp), pointer :: omega_end(:,:),rho_end(:,:)
    complex(pfdp), pointer :: omega_ex(:,:), rho_ex(:,:)
    real(pfdp) :: maxerr_om,maxerr_rho,t,t0


    sweeper => as_my_sweeper(pf%levels(level_index)%ulevel%sweeper)    

    !  Get the solution at the end of time step
    omega_end => get_array2d(pf%levels(level_index)%qend,1)
    rho_end => get_array2d(pf%levels(level_index)%qend,2)

    !  Use sweeper temp space for exact solutions 
    omega_ex => sweeper%tmp
    rho_ex => sweeper%tmphat
    rho_ex => sweeper%u

    !>  compute the exact solution
    t=pf%state%t0+pf%state%dt
    t0=pf%state%t0
    call sweeper%exact(t, omega_ex,rho_ex)
    
    !>  compute error
    maxerr_om = maxval(abs(omega_end-omega_ex))
    maxerr_rho = maxval(abs(rho_end-rho_ex))
    
    print '("error: time: ", f8.4," step: ",i3.3," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7,es14.7," resid: ",es14.7)', &
         t,pf%state%step+1, pf%rank, pf%state%iter,level_index,maxerr_om,maxerr_rho,pf%levels(level_index)%residual    


    call flush(6)

    
  end subroutine echo_error

end module hooks
