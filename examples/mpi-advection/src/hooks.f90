!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
contains

  subroutine echo_error(pf, level, state)
    use iso_c_binding
    use feval, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)

    qend => array1(level%qend)

    call exact(state%t0+state%dt, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(qend-yexact))
  end subroutine echo_error


  subroutine echo_residual(pf, level, state)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    real(pfdp), pointer :: r(:)
    integer :: np
    r => array1(level%R(level%nnodes-1))

    print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
         state%step+1, state%iter, level%level, maxval(abs(r))


!!$    r => array1(level%Q(1))
!!$    print *,'Q1',r
!!$    r => array1(level%Q(level%nnodes))
!!$    print *,'Qend',r
!!$    r => array1(level%F(1,1))
!!$    print *,'F1',r
!!$    r => array1(level%F(level%nnodes,1))
!!$    print *,'F1end',r
!!$    r => array1(level%S(level%nnodes-1))
!!$    print *,'S',r
!!$
!!$    np = size(level%F(level%nnodes,:))
!!$    print *,'npieces',np
!!$    if (np > 1) then
!!$       r => array1(level%F(1,2))
!!$       print *,'F21',r
!!$       r => array1(level%F(level%nnodes,2))
!!$       print *,'F2end',r
!!$    end if
  end subroutine echo_residual


end module hooks
