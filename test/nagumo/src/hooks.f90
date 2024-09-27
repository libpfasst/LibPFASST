!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray_oc
  use probin, only: nprob,poutmod, fbase, foutbase, N_Vcycles,kfreq
  implicit none
contains

subroutine echo_error_hook(pf, level_index)
    use pf_mod_utils
    use solutions
    type(pf_pfasst_t),   intent(inout) :: pf
    integer, intent(in) :: level_index

!     class(pf_level_t),    intent(inout) :: level
!     type(pf_state_t),    intent(in)    :: state

!     real(c_double) :: yexact(product(level%shape)), pexact(product(level%shape))
    real(pfdp), pointer :: qend(:), q0(:)
    real(pfdp) :: res,max_y,max_p, err_y, err_p, t
    integer :: un, Nx, Nxy, m
    character(len=64) :: fout 
    character(len=7) :: stepstring
    character(len=3) :: Vstring
    character(len=3) :: Kstring
    real(pfdp) :: pnorms(pf%levels(level_index)%nnodes-1), ynorms(pf%levels(level_index)%nnodes-1)


    qend => get_array1d_oc(pf%levels(level_index)%Q(pf%levels(level_index)%nnodes),1) 
    q0   => get_array1d_oc(pf%levels(level_index)%Q(1),2) 
    t = pf%state%t0+pf%state%dt   
!     call exact_y(t, level%nvars, yexact)
    max_y=maxval(abs(qend))
!     err_y = maxval(abs(qend-yexact))
    
!     call exact_p(state%t0, level%nvars, pexact)
    max_p = maxval(abs(q0))
!     err_p = maxval(abs(q0-pexact))
    
    
    do m = 1, pf%levels(level_index)%nnodes-1
      pnorms(m) = pf%levels(level_index)%R(m)%norm(2)
      ynorms(m) = pf%levels(level_index)%R(m)%norm(1)
    end do

    res = min(pf%levels(level_index)%residual, pf%levels(level_index)%residual_rel)
    !res = pf%levels(level_index)%residual
    print '(" rank:",i5," lev:",i5," step:",i5," iter:",i3," res:",es13.6)', &
          pf%rank, level_index, pf%state%step+1, pf%state%iter, res


!     un = 1000+state%step+1
!     write(stepstring,"(I0.3)") state%step+1
!     write(Vstring,"(I0.2)") N_Vcycles
!     write(Kstring,"(I0.3)") kfreq
!     fout = trim(foutbase)//'N_V'//trim(Vstring)//'N_step'//trim(stepstring)//'K'//trim(Kstring)//'.m'
!     open(unit=un, file = fout, status = 'unknown', action = 'write', position='append')
!     write(un,*)  level%level,state%step+1,t,state%iter,max_y, err_y, max_p, err_p, res
!     close(un)

  end subroutine echo_error_hook


  subroutine echo_residual_hook(pf, level_index)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp), pointer :: ry(:), rp(:)

    ry => get_array1d_oc(pf%levels(level_index)%R(pf%levels(level_index)%nnodes-1),1)
    rp => get_array1d_oc(pf%levels(level_index)%R(2),2)

    print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," res_y: ",es14.7," res_p: ",es14.7)', &
         pf%state%step+1, pf%state%iter, level_index, maxval(abs(ry)), maxval(abs(rp))
  end subroutine echo_residual_hook


end module hooks
