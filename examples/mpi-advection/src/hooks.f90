!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pfasst
  use encap_array1d
  implicit none
contains

  subroutine echo_error(pf, level, state, ctx)
    use iso_c_binding
    use feval, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)

    qend => array(level%qend)

    call exact(state%t0+state%dt, level%nvars, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(qend-yexact))
  end subroutine echo_error

end module hooks
