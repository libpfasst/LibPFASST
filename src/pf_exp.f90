!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module pf_mod_sweep
  use pf_mod_dtype
  implicit none
  integer, parameter :: npieces = 1
contains

  ! Perform on SDC sweep on level F and set qend appropriately.
  subroutine sweep(pf, t0, dt, F)
    use pf_mod_timer
    use feval, only : eval_f1

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:F%nnodes-1)
    type(pf_encap_t) :: S(F%nnodes-1)

    call start_timer(pf, TLEVEL+F%level-1)

    ! compute integrals and add fas correction
    do m = 1, F%nnodes-1
       call create(S(m), F%level, .false., F%nvars, F%shape, F%ctx)
       call setval(S(m), 0.0d0)
       do n = 1, F%nnodes
          call axpy(S(m), dt*F%smat(m,n,1), F%fSDC(n,1))
       end do
       if (associated(F%tau)) then
          call axpy(S(m), 1.0d0, F%tau(m))
       end if
    end do

    ! do the time-stepping
    call unpack(F%qSDC(1), F%q0)

    call eval_f1(F%qSDC(1), t0, F%level, F%ctx, F%fSDC(1,1))

    t = t0
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do m = 1, F%nnodes-1
       t = t + dtsdc(m)

       call copy(F%qSDC(m+1), F%qSDC(m))
       call axpy(F%qSDC(m+1), dtsdc(m), F%fSDC(m,1))
       call axpy(F%qSDC(m+1), 1.0d0, S(m))

       call eval_f1(F%qSDC(m+1), t, F%level, F%ctx, F%fSDC(m+1,1))
    end do

    call copy(F%qend, F%qSDC(F%nnodes))

    ! done
    do m = 1, F%nnodes-1
       call destroy(S(m))
    end do

    call end_timer(pf, TLEVEL+F%level-1)
  end subroutine sweep

  ! Evaluate function values
  subroutine sdceval(t, m, F)
    use pf_mod_dtype
    use feval, only: eval_f1

    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F

    call eval_f1(F%qSDC(m), t, F%level, F%ctx, F%fSDC(m,1))
  end subroutine sdceval

  ! Initialize smats
  subroutine sdcinit(F)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: F
    real(pfdp) :: dsdc(F%nnodes-1)

    integer :: m

    allocate(F%smat(F%nnodes-1,F%nnodes,npieces))

    F%smat(:,:,1) = F%s0mat

    dsdc = F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1)
    do m = 1, F%nnodes-1
       F%smat(m,m,1) = F%smat(m,m,1)   - dsdc(m)
    end do
  end subroutine sdcinit


  ! Compute SDC integral
  subroutine sdc_integrate(qSDC, fSDC, dt, F, fintSDC)
    type(pf_level_t),  intent(in)    :: F
    type(pf_encap_t),  intent(in)    :: qSDC(F%nnodes,npieces) ! Solution
    type(pf_encap_t),  intent(in)    :: fSDC(F%nnodes,npieces) ! Function values
    type(pf_encap_t),  intent(inout) :: fintSDC(F%nnodes-1)    ! Integrals of f
    real(pfdp),        intent(in)    :: dt

    integer :: n, m, p

    do n = 1, F%nnodes-1
       call setval(fintSDC(n), 0.0d0)
       do m = 1, F%nnodes
          do p = 1, npieces
             call axpy(fintSDC(n), dt*F%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine sdc_integrate

end module pf_mod_sweep

