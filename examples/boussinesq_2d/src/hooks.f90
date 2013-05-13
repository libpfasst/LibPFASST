!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pfasst
  use encap_array1d
  use spatialdiscretization, only : WriteData, Nx, Ny, nr_fields
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
    print '("error: step: ",i4.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(qend-yexact))
  end subroutine echo_error

  subroutine echo_residual(pf, level, state, ctx)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    real(pfdp), pointer :: r(:)
    type(c_ptr) :: residual

   r => array(level%R(level%nnodes-1))

    print '("resid: step: ",i3.3," iter: ",i3.3," level: ",i2.2," resid: ",es14.7)', &
         state%step+1, state%iter, level%level, maxval(abs(r))

  end subroutine echo_residual

  subroutine output(pf, level, state, ctx)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    real(pfdp), pointer :: y(:)
    double precision, allocatable, dimension(:) :: residuals

    y => array(level%qend)

    allocate(residuals(pf%niters)) ! must match maxit
    residuals = 0.0

    Call WriteData( residuals, DBLE(0.0), DBLE(0.0), DBLE(0.0), pf%niters, state%step+1, RESHAPE(y, (/ nr_fields, Ny, Nx /)) )

  end subroutine output
end module hooks
