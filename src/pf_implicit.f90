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

module pf_mod_implicit
  use pf_mod_dtype
  implicit none
  integer, parameter, private :: npieces = 1

  interface
     subroutine pf_f2eval_p(y, t, level, ctx, f2)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, f2, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2eval_p
  end interface

  interface
     subroutine pf_f2comp_p(y, t, dt, rhs, level, ctx, f2)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, rhs, f2, ctx
       real(pfdp),     intent(in)        :: t, dt
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2comp_p
  end interface

  interface
     subroutine pf_imp_rhs_p(rhs, q0, S, level, ctx)
       import c_ptr
       type(c_ptr), intent(in), value :: rhs, q0, S, ctx
       integer,     intent(in)        :: level
     end subroutine pf_imp_rhs_p
  end interface

  type :: pf_implicit_t
     procedure(pf_f2eval_p),  pointer, nopass :: f2eval
     procedure(pf_f2comp_p),  pointer, nopass :: f2comp
     procedure(pf_imp_rhs_p), pointer, nopass :: gen_rhs
  end type pf_implicit_t

contains

  ! Perform one SDC sweep on level F and set qend appropriately.
  subroutine implicit_sweep(pf, F, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:F%nnodes-1)
    type(c_ptr) :: rhs

    type(pf_implicit_t), pointer :: imp

    call c_f_pointer(F%sweeper%ctx, imp)

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

    call imp%f2eval(F%Q(1), t0, F%level, F%ctx, F%F(1,1))

    call F%encap%create(rhs, F%level, SDC_KIND_SOL_FEVAL, F%nvars, F%shape, F%ctx, F%encap%ctx)

    t = t0
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do m = 1, F%nnodes-1
       t = t + dtsdc(m)

       if (associated(imp%gen_rhs)) then
          call imp%gen_rhs(rhs, F%Q(m), F%S(m), F%level, F%ctx)
       else
          call F%encap%copy(rhs, F%Q(m))
          call F%encap%axpy(rhs, 1.0d0, F%S(m))
       end if

       call imp%f2comp(F%Q(m+1), t, dtsdc(m), rhs, F%level, F%ctx, F%F(m+1,1))
    end do

    call F%encap%copy(F%qend, F%Q(F%nnodes))

    ! done
    call F%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+F%level-1)
  end subroutine implicit_sweep

  ! Evaluate function values
  subroutine implicit_evaluate(F, t, m)
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F

    type(pf_implicit_t), pointer :: imp
    call c_f_pointer(F%sweeper%ctx, imp)

    call imp%f2eval(F%Q(m), t, F%level, F%ctx, F%F(m,1))
  end subroutine implicit_evaluate

  ! Initialize smats
  subroutine implicit_initialize(F)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: F

    real(pfdp) :: dsdc(F%nnodes-1)

    integer :: m

    allocate(F%smat(F%nnodes-1,F%nnodes,npieces))

    F%smat(:,:,1) = F%s0mat

    dsdc = F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1)
    do m = 1, F%nnodes-1
       F%smat(m,m+1,1) = F%smat(m,m+1,1) - dsdc(m)
    end do
  end subroutine implicit_initialize


  ! Compute SDC integral
  subroutine implicit_integrate(F, qSDC, fSDC, dt, fintSDC)
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
  end subroutine implicit_integrate

  ! Create implicit sweeper
  subroutine implicit_create(sweeper, f2eval, f2comp)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp

    type(pf_implicit_t), pointer :: imp

    allocate(imp)
    imp%f2eval => f2eval
    imp%f2comp => f2comp

    sweeper%sweep      => implicit_sweep
    sweeper%evaluate   => implicit_evaluate
    sweeper%initialize => implicit_initialize
    sweeper%integrate  => implicit_integrate

    sweeper%ctx = c_loc(imp)
  end subroutine implicit_create

end module pf_mod_implicit

