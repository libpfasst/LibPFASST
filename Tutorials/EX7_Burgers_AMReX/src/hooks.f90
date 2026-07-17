!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
  use pf_mod_AMReX_mfab
  use pf_mod_rutils
  use pf_mod_mpi
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
    real(pfdp) :: errL2(1), errL1(1), errLinf(1)
    
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
    errL2  = y_ex%normL2()
    errL1  = y_ex%normL1()
    errLinf= y_ex%normLinf()   
    
    print '("error: time: ", f4.2," step: ",i2.1," rank: ",i2.2," iter: ",i2.2," level: ",i2.2," err0: ",es12.5," err1: ",es12.5," err2: ",es12.5," resid: ",es12.5," relresid: ",es12.5)', &
         time,step, rank, iter,level_index,errLinf,errL1,errL2,resid,relresid


    call flush(6)

    call amrex_mfab_destroy(y_ex)

    call pf_set_error(pf,level_index,errL2(1))
  
  end subroutine echo_error

  !>  Output the error and residual in the solution
  subroutine write_plotfile(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    !>  Local
    real(pfdp)          :: time
    integer             :: rank, iter, step 
    character(len=4)    :: step_str, iter_str
    type(pf_amrex_mfab_t) :: y_ex      

    ! compute vars & set pointers
    step     = pf%state%step+1
    time     = pf%state%t0+pf%state%dt
    iter     = pf%state%iter

    ! Convert integers to strings BEFORE using them
    write(step_str, '(I0)') step
    write(iter_str, '(I0)') iter
    
    ! write
    select type (qend_ptr => pf%levels(level_index)%qend)
    type is (pf_amrex_mfab_t)
      call qend_ptr%txtprint("dat/"//trim(pf%outdir)//"/plt_file_intermediate_step"//trim(step_str)//"_iter"//trim(iter_str), &
                              time, step, iter)      
    end select
    
    select type (am_fac => pf%levels(level_index)%ulevel%factory)
    type is (pf_AMReX_mfab_factory_t)
      !>  build mfab to store exact solution
      call amrex_mfab_build(y_ex, pf%levels(level_index)%lev_shape, am_fac%ba, am_fac%dm, am_fac%geom)
    end select
    !>  compute the exact solution
    call exact(time,y_ex)
    call y_ex%txtprint("dat/"//trim(pf%outdir)//"/plt_file_real_step"//trim(step_str)//"_iter"//trim(iter_str), &
                              time, step, iter) 
    call amrex_mfab_destroy(y_ex)


  end subroutine write_plotfile

end module hooks
