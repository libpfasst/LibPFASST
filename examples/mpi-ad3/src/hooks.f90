!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  use probin, only: nprob,poutmod, fbase, foutbase, N_Vcycles,kfreq
  implicit none
contains

  subroutine echo_error_hook(pf, level, state, levelctx)
    use pf_mod_utils
    use solutions, only: exact,exact_ode
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: levelctx

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)
    real(pfdp) :: res,max_y,err,t,ODE_err
    integer :: un
    character(len=64) :: fout 
    character(len=7) :: stepstring
    character(len=3) :: Vstring
    character(len=3) :: Kstring
!    qend => array1(level%qend)
    qend => array1(level%Q(level%Nnodes))
    t = state%t0+state%dt   
    call exact(t, level%nvars, yexact)
!    max_y=maxval(abs(yexact))
    max_y=maxval(abs(qend))
    err = maxval(abs(qend-yexact))
    call exact_ode(t, level%nvars, yexact)
    ODE_err = maxval(abs(qend-yexact))


    call pf_residual(pf, level, state%dt)
    res= level%residual
    print '(" lev:",i5," step:",i5," t=",es10.3," iter:",i3," Max_y:",es13.6," Err:",es13.6," ODEERR:",es13.6," Res:",es13.6)', &
               level%level,state%step+1, t,state%iter, max_y,err,ODE_err,res

    un = 1000+state%step+1
    write(stepstring,"(I0.3)") state%step+1
    write(Vstring,"(I0.2)") N_Vcycles
    write(Kstring,"(I0.3)") kfreq
    fout = trim(foutbase)//'N_V'//trim(Vstring)//'N_step'//trim(stepstring)//'K'//trim(Kstring)//'.m'
    open(unit=un, file = fout, status = 'unknown', action = 'write', position='append')
    write(un,*)  level%level,state%step+1,t,state%iter,max_y,err,ODE_err,res
    close(un)

  end subroutine echo_error_hook

  subroutine echo_residual_hook(pf, level, state, levelctx)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: levelctx

    real(pfdp), pointer :: r(:)

    r => array1(level%R(level%nnodes-1))

    print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
         state%step+1, state%iter, level%level, maxval(abs(r))
  end subroutine echo_residual_hook

end module hooks
