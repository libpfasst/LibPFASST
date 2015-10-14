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

module pf_mod_explicit
  use pf_mod_dtype
  use pf_mod_utils
  implicit none
  integer, parameter, private :: npieces = 1

  type, extends(pf_sweeper_t), abstract :: pf_explicit_t
     real(pfdp), allocatable :: SdiffE(:,:)
   contains
     procedure(pf_f1eval_p), deferred :: f1eval
     procedure :: sweep => explicit_sweep
     procedure :: initialize => explicit_initialize
     procedure :: evaluate => explicit_evaluate
     procedure :: integrate => explicit_integrate
     procedure :: residual => explicit_residual
     procedure :: evaluate_all => explicit_evaluate_all
  end type pf_explicit_t

  interface
     subroutine pf_f1eval_p(this, y, t, level, levelctx, f1)
       import pf_explicit_t, c_ptr, c_int, pfdp
       class(pf_explicit_t), intent(in) :: this
       type(c_ptr),    intent(in), value :: y, f1, levelctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f1eval_p
  end interface

contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine explicit_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer

    class(pf_explicit_t), intent(in) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:lev%nnodes-1)

    call start_timer(pf, TLEVEL+lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1
       call lev%encap%setval(lev%S(m), 0.0_pfdp)
       do n = 1, lev%nnodes
          call lev%encap%axpy(lev%S(m), dt*this%SdiffE(m,n), lev%F(n,1))
       end do
       if (associated(lev%tau)) then
          call lev%encap%axpy(lev%S(m), 1.0_pfdp, lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call lev%encap%unpack(lev%Q(1), lev%q0)

    call this%f1eval(lev%Q(1), t0, lev%level, lev%levelctx, lev%F(1,1))

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       call lev%encap%copy(lev%Q(m+1), lev%Q(m))
       call lev%encap%axpy(lev%Q(m+1), dtsdc(m), lev%F(m,1))
       call lev%encap%axpy(lev%Q(m+1), 1.0_pfdp, lev%S(m))

       call this%f1eval(lev%Q(m+1), t, lev%level, lev%levelctx, lev%F(m+1,1))
    end do

    call lev%encap%copy(lev%qend, lev%Q(lev%nnodes))

    ! done
    call end_timer(pf, TLEVEL+lev%level-1)
  end subroutine explicit_sweep

  ! Evaluate function values
  subroutine explicit_evaluate(this, lev, t, m)
    use pf_mod_dtype

    class(pf_explicit_t), intent(in) :: this
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: lev

    call this%f1eval(lev%Q(m), t, lev%level, lev%levelctx, lev%F(m,1))
  end subroutine explicit_evaluate

  ! Initialize matrix
  subroutine explicit_initialize(this, lev)
    use pf_mod_dtype
    class(pf_explicit_t), intent(inout) :: this
    type(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)

    integer :: m,nnodes

    nnodes = lev%nnodes
    allocate(this%SdiffE(nnodes-1,nnodes))  ! S-FE

    this%SdiffE = lev%s0mat

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       this%SdiffE(m,m)   = this%SdiffE(m,m)   - dsdc(m)
    end do
  end subroutine explicit_initialize


  ! Compute SDC integral
  subroutine explicit_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_explicit_t), intent(in) :: this
    type(pf_level_t), intent(in)    :: lev
    type(c_ptr),      intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),       intent(in)    :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call lev%encap%setval(fintSDC(n), 0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, npieces
             call lev%encap%axpy(fintSDC(n), dt*lev%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine explicit_integrate

  subroutine explicit_residual(this, Lev, dt)
    class(pf_explicit_t), intent(in)  :: this
    type(pf_level_t),  intent(inout) :: Lev
    real(pfdp),        intent(in)    :: dt

    integer :: m, n

    call pf_generic_residual(this, Lev, dt)
  end subroutine explicit_residual

  subroutine explicit_evaluate_all(this, Lev, t)
    class(pf_explicit_t), intent(in)  :: this
    type(pf_level_t),  intent(inout) :: Lev
    real(pfdp),        intent(in)    :: t(:)

    call pf_generic_evaluate_all(this, Lev, t)
  end subroutine explicit_evaluate_all

end module pf_mod_explicit

