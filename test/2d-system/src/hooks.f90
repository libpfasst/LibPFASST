!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
  use pf_mod_ndsysarray
  implicit none
contains

  subroutine echo_error(pf, level, state)
    use iso_c_binding
    use feval, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    type(ndsysarray) :: y_exact

    real(pfdp), pointer :: u_err(:,:),v_err(:,:),r(:,:)
    real(pfdp)  :: toterr

    r => get_array2d(level%R(level%nnodes-1),1)

    call ndsysarray_build(y_exact,  level%shape)
    call exact(state%t0+state%dt, y_exact)




    call y_exact%axpy(-1.0_pfdp,level%qend)
    u_err => get_array2d(y_exact,1)
    v_err => get_array2d(y_exact,2)

    toterr =maxval(abs(u_err)) + maxval(abs(v_err))

    print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
         state%step+1, state%iter,level%index, toterr,maxval(abs(r))
    call flush
  end subroutine echo_error



end module hooks
