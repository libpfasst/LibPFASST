!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

module pf_mod_sweep
  use pf_mod_dtype
  implicit none
  integer, parameter :: npieces = 2
contains

  ! Perform on SDC sweep on level F and set qend appropriately.
  subroutine sweep(pf, t0, dt, F)
    use pf_mod_timer
    use feval, only : eval_f1, eval_f2, comp_f2

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:F%nnodes-1)
    type(pf_encap_t) :: S(F%nnodes-1), rhs

    call start_timer(pf, TLEVEL+F%level-1)

    ! compute integrals and add fas correction
    do m = 1, F%nnodes-1
       call create(S(m), F%level, .false., F%nvars, F%shape, F%ctx)
       call setval(S(m), 0.0d0)
       do n = 1, F%nnodes
          call axpy(S(m), dt*F%smat(m,n,1), F%fSDC(n,1))
          call axpy(S(m), dt*F%smat(m,n,2), F%fSDC(n,2))
       end do
       if (associated(F%tau)) then
          call axpy(S(m), 1.0d0, F%tau(m))
       end if
    end do

    ! do the time-stepping
    call unpack(F%qSDC(1), F%q0)

    call eval_f1(F%qSDC(1), t0, F%level, F%ctx, F%fSDC(1,1))
    call eval_f2(F%qSDC(1), t0, F%level, F%ctx, F%fSDC(1,2))

    call create(rhs, F%level, .false., F%nvars, F%shape, F%ctx)

    t = t0
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do m = 1, F%nnodes-1
       t = t + dtsdc(m)

       if (associated(F%gen_imex_rhs)) then
          call F%gen_imex_rhs(rhs, F%qSDC(m), dtsdc(m), F%fSDC(m,1), S(m), F%level, F%ctx)
       else
          call copy(rhs, F%qSDC(m))
          call axpy(rhs, dtsdc(m), F%fSDC(m,1))
          call axpy(rhs, 1.0d0, S(m))
       end if

       call comp_f2(F%qSDC(m+1), t, dtsdc(m), rhs, F%level, F%ctx, F%fSDC(m+1,2))
       call eval_f1(F%qSDC(m+1), t, F%level, F%ctx, F%fSDC(m+1,1))
    end do

    call copy(F%qend, F%qSDC(F%nnodes))

    ! done
    call destroy(rhs)
    do m = 1, F%nnodes-1
       call destroy(S(m))
    end do

    call end_timer(pf, TLEVEL+F%level-1)
  end subroutine sweep

  ! Evaluate function values
  subroutine sdceval(t, m, F)
    use feval, only: eval_f1, eval_f2

    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F

    call eval_f1(F%qSDC(m), t, F%level, F%ctx, F%fSDC(m,1))
    call eval_f2(F%qSDC(m), t, F%level, F%ctx, F%fSDC(m,2))
  end subroutine sdceval

  ! Initialize smats
  subroutine sdcinit(F)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: F
    real(pfdp) :: dsdc(F%nnodes-1)

    integer :: m

    allocate(F%smat(F%nnodes-1,F%nnodes,npieces))

    F%smat(:,:,1) = F%s0mat
    F%smat(:,:,2) = F%s0mat

    dsdc = F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1)
    do m = 1, F%nnodes-1
       F%smat(m,m,1)   = F%smat(m,m,1)   - dsdc(m)
       F%smat(m,m+1,2) = F%smat(m,m+1,2) - dsdc(m)
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

