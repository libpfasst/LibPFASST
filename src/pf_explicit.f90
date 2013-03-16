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
  implicit none
  integer, parameter, private :: npieces = 1

  interface
     subroutine pf_f1eval_p(y, t, level, ctx, f1)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, f1, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f1eval_p
  end interface

  type :: pf_explicit_t
     procedure(pf_f1eval_p),   pointer, nopass :: f1eval
  end type pf_explicit_t

contains

  ! Perform on SDC sweep on level F and set qend appropriately.
  subroutine explicit_sweep(pf, F, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:F%nnodes-1)
    type(pf_explicit_t), pointer :: exp

    call c_f_pointer(F%sweeper%ctx, exp)

    call start_timer(pf, TLEVEL+F%level-1)

    ! compute integrals and add fas correction
    do m = 1, F%nnodes-1
       call F%encap%setval(F%S(m), 0.0d0)
       do n = 1, F%nnodes
          call F%encap%axpy(F%S(m), dt*F%smat(m,n,1), F%F(n,1))
       end do
       if (associated(F%tau)) then
          call F%encap%axpy(F%S(m), 1.0d0, F%tau(m))
       end if
    end do

    ! do the time-stepping
    call F%encap%unpack(F%Q(1), F%q0)

    call exp%f1eval(F%Q(1), t0, F%level, F%ctx, F%F(1,1))

    t = t0
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do m = 1, F%nnodes-1
       t = t + dtsdc(m)

       call F%encap%copy(F%Q(m+1), F%Q(m))
       call F%encap%axpy(F%Q(m+1), dtsdc(m), F%F(m,1))
       call F%encap%axpy(F%Q(m+1), 1.0d0, F%S(m))

       call exp%f1eval(F%Q(m+1), t, F%level, F%ctx, F%F(m+1,1))
    end do

    call F%encap%copy(F%qend, F%Q(F%nnodes))

    ! done
    call end_timer(pf, TLEVEL+F%level-1)
  end subroutine explicit_sweep

  ! Evaluate function values
  subroutine explicit_evaluate(F, t, m)
    use pf_mod_dtype

    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F

    type(pf_explicit_t), pointer :: exp

    call c_f_pointer(F%sweeper%ctx, exp)

    call exp%f1eval(F%Q(m), t, F%level, F%ctx, F%F(m,1))
  end subroutine explicit_evaluate

  ! Initialize smats
  subroutine explicit_initialize(F)
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
  end subroutine explicit_initialize


  ! Compute SDC integral
  subroutine explicit_integrate(F, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in) :: F
    type(c_ptr),      intent(in) :: qSDC(:, :), fSDC(:, :), fintSDC(:)
    real(pfdp),       intent(in) :: dt

    integer :: n, m, p

    do n = 1, F%nnodes-1
       call F%encap%setval(fintSDC(n), 0.0d0)
       do m = 1, F%nnodes
          do p = 1, npieces
             call F%encap%axpy(fintSDC(n), dt*F%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine explicit_integrate

  ! Create implicit sweeper
  subroutine explicit_create(sweeper, f1eval)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval

    type(pf_explicit_t), pointer :: exp

    allocate(exp)
    exp%f1eval => f1eval

    sweeper%sweep => explicit_sweep
    sweeper%evaluate => explicit_evaluate
    sweeper%initialize => explicit_initialize
    sweeper%integrate => explicit_integrate

    sweeper%ctx = c_loc(exp)
  end subroutine explicit_create


end module pf_mod_explicit

