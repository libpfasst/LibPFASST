!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
  use pf_mod_AMReX_mfab
  use pf_mod_rutils
  implicit none
contains
 
  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    use pf_my_sweeper, only: my_sweeper_t, as_my_sweeper
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    !>  Local
    class(my_sweeper_t), pointer :: sweeper
    type(pf_amrex_mfab_t) :: y_ex      

    real(pfdp) :: maxerr, time, resid, relresid
    integer    :: rank, iter, step 

    ! compute vars & set pointers
    step     = pf%state%step+1
    time     = pf%state%t0+pf%state%dt
    iter     = pf%state%iter
    rank     = pf%rank
    resid    = pf%levels(level_index)%residual
    relresid = pf%levels(level_index)%residual_rel
    
    !> build mfab to store exact solution
    select type (am_fac => pf%levels(level_index)%ulevel%factory)
    type is (pf_AMReX_mfab_factory_t)
      !>  build mfab to store exact solution
      call amrex_mfab_build(y_ex, pf%levels(level_index)%lev_shape, am_fac%ba, am_fac%dm, am_fac%geom)
    end select
    sweeper => as_my_sweeper(pf%levels(level_index)%ulevel%sweeper)    

    !>  compute the exact solution
    call exact(time,y_ex)
    
    !>  compute error
    call y_ex%axpy(-1.0d0,pf%levels(level_index)%qend) ! computes exact - final(qend)
    maxerr = y_ex%norm()   
    
    print '("error: time: ", f8.4," step: ",i8.1," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," resid: ",es14.7," relresid: ",es14.7)', &
         time,step, rank, iter,level_index,maxerr,resid,relresid


    call flush(6)

    call amrex_mfab_destroy(y_ex)

    call pf_set_error(pf,level_index,maxerr)
  
  end subroutine echo_error

end module hooks
