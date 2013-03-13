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

module pf_mod_imex
  use pf_mod_dtype
  implicit none
  integer, parameter :: npieces = 2

  interface
     subroutine pf_f1eval_p(y, t, level, ctx, f1)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, f1, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f1eval_p
  end interface

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
     subroutine pf_imex_rhs_p(rhs, q0, dt, f1, S, level, ctx)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value :: rhs, q0, f1, S, ctx
       real(pfdp),  intent(in)        :: dt
       integer,     intent(in)        :: level
     end subroutine pf_imex_rhs_p
  end interface

  type :: pf_imex_t
     procedure(pf_f1eval_p),   pointer, nopass :: f1eval
     procedure(pf_f2eval_p),   pointer, nopass :: f2eval
     procedure(pf_f2comp_p),   pointer, nopass :: f2comp
     procedure(pf_imex_rhs_p), pointer, nopass :: gen_rhs
  end type pf_imex_t

contains

  ! Perform on SDC sweep on level F and set qend appropriately.
  subroutine imex_sweep(pf, F, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:F%nnodes-1)
    type(c_ptr) :: S(F%nnodes-1), rhs

    type(pf_imex_t), pointer :: imex

    call c_f_pointer(F%sweeper%ctx, imex)

    call start_timer(pf, TLEVEL+F%level-1)

    ! compute integrals and add fas correction
    do m = 1, F%nnodes-1
       call F%encap%create(S(m), F%level, .false., F%nvars, F%shape, F%ctx)
       call F%encap%setval(S(m), 0.0d0)
       do n = 1, F%nnodes
          call F%encap%axpy(S(m), dt*F%smat(m,n,1), F%fSDC(n,1))
          call F%encap%axpy(S(m), dt*F%smat(m,n,2), F%fSDC(n,2))
       end do
       if (associated(F%tau)) then
          call F%encap%axpy(S(m), 1.0d0, F%tau(m))
       end if
    end do

    ! do the time-stepping
    call F%encap%unpack(F%qSDC(1), F%q0)

    call imex%f1eval(F%qSDC(1), t0, F%level, F%ctx, F%fSDC(1,1))
    call imex%f2eval(F%qSDC(1), t0, F%level, F%ctx, F%fSDC(1,2))

    call F%encap%create(rhs, F%level, .false., F%nvars, F%shape, F%ctx)

    t = t0
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do m = 1, F%nnodes-1
       t = t + dtsdc(m)

       if (associated(imex%gen_rhs)) then
          call imex%gen_rhs(rhs, F%qSDC(m), dtsdc(m), F%fSDC(m,1), S(m), F%level, F%ctx)
       else
          call F%encap%copy(rhs, F%qSDC(m))
          call F%encap%axpy(rhs, dtsdc(m), F%fSDC(m,1))
          call F%encap%axpy(rhs, 1.0d0, S(m))
       end if

       call imex%f2comp(F%qSDC(m+1), t, dtsdc(m), rhs, F%level, F%ctx, F%fSDC(m+1,2))
       call imex%f1eval(F%qSDC(m+1), t, F%level, F%ctx, F%fSDC(m+1,1))
    end do

    call F%encap%copy(F%qend, F%qSDC(F%nnodes))

    ! done
    call F%encap%destroy(rhs)
    do m = 1, F%nnodes-1
       call F%encap%destroy(S(m))
    end do

    call end_timer(pf, TLEVEL+F%level-1)
  end subroutine imex_sweep

  ! Evaluate function values
  subroutine imex_evaluate(F, t, m)
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F

    type(pf_imex_t), pointer :: imex
    call c_f_pointer(F%sweeper%ctx, imex)

    call imex%f1eval(F%qSDC(m), t, F%level, F%ctx, F%fSDC(m,1))
    call imex%f2eval(F%qSDC(m), t, F%level, F%ctx, F%fSDC(m,2))
  end subroutine imex_evaluate

  ! Initialize smats
  subroutine imex_initialize(F)
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
  end subroutine imex_initialize

  ! Compute SDC integral
  subroutine imex_integrate(F, fSDC, dt, fintSDC)
    type(pf_level_t),  intent(in) :: F
    type(c_ptr),       intent(in) :: fSDC(:, :), fintSDC(:)
    real(pfdp),        intent(in) :: dt

    integer :: n, m, p

    do n = 1, F%nnodes-1
       call F%encap%setval(fintSDC(n), 0.0d0)
       do m = 1, F%nnodes
          do p = 1, npieces
             call F%encap%axpy(fintSDC(n), dt*F%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine imex_integrate

  ! Create IMEX sweeper
  subroutine imex_create(sweeper, f1eval, f2eval, f2comp)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp

    type(pf_imex_t), pointer :: imex

    allocate(imex)
    imex%f1eval => f1eval
    imex%f2eval => f2eval
    imex%f2comp => f2comp

    sweeper%sweep => imex_sweep
    sweeper%evaluate => imex_evaluate
    sweeper%initialize => imex_initialize
    sweeper%integrate => imex_integrate

    sweeper%ctx = c_loc(imex)
  end subroutine imex_create


end module pf_mod_imex

