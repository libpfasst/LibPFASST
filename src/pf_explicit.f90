!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

module pf_mod_explicit
  use pf_mod_dtype
  use pf_mod_utils
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
     real(pfdp), ALLOCATABLE :: SdiffE(:,:)
  end type pf_explicit_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine explicit_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:Lev%nnodes-1)
    type(pf_explicit_t), pointer :: exp

    call c_f_pointer(Lev%sweeper%sweeperctx, exp)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0_pfdp)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*exp%SdiffE(m,n), Lev%F(n,1))
       end do
       if (associated(Lev%tau)) then
          call Lev%encap%axpy(Lev%S(m), 1.0_pfdp, Lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%copy(Lev%Q(1), Lev%q0)

    call exp%f1eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,1))

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%copy(Lev%Q(m+1), Lev%Q(m))
       call Lev%encap%axpy(Lev%Q(m+1), dtsdc(m), Lev%F(m,1))
       call Lev%encap%axpy(Lev%Q(m+1), 1.0_pfdp, Lev%S(m))

       call exp%f1eval(Lev%Q(m+1), t, Lev%level, Lev%ctx, Lev%F(m+1,1))
    end do

    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine explicit_sweep

  ! Evaluate function values
  subroutine explicit_evaluate(Lev, t, m)
    use pf_mod_dtype

    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_explicit_t), pointer :: exp

    call c_f_pointer(Lev%sweeper%sweeperctx, exp)

    call exp%f1eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,1))
  end subroutine explicit_evaluate

  ! Initialize matrix
  subroutine explicit_initialize(Lev)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: Lev
    real(pfdp) :: dsdc(Lev%nnodes-1)

    integer :: m,nnodes
    type(pf_explicit_t), pointer :: exp
    call c_f_pointer(Lev%sweeper%sweeperctx, exp)

    nnodes = Lev%nnodes
    allocate(exp%SdiffE(nnodes-1,nnodes))  ! S-FE

    exp%SdiffE = Lev%s0mat

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       exp%SdiffE(m,m)   = exp%SdiffE(m,m)   - dsdc(m)
    end do
  end subroutine explicit_initialize


  ! Compute SDC integral
  subroutine explicit_integrate(Lev, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in)    :: Lev
    type(c_ptr),      intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),       intent(in)    :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, Lev%nnodes-1
       call Lev%encap%setval(fintSDC(n), 0.0_pfdp)
       do m = 1, Lev%nnodes
          do p = 1, npieces
             call Lev%encap%axpy(fintSDC(n), dt*Lev%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine explicit_integrate

  ! Create explicit sweeper
  subroutine pf_explicit_create(sweeper, f1eval)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval

    type(pf_explicit_t), pointer :: exp

    allocate(exp)
    exp%f1eval => f1eval

    sweeper%npieces = npieces
    sweeper%sweep      => explicit_sweep
    sweeper%evaluate   => explicit_evaluate
    sweeper%initialize => explicit_initialize
    sweeper%integrate  => explicit_integrate
    sweeper%destroy    => pf_explicit_destroy
    sweeper%evaluate_all => pf_generic_evaluate_all
    sweeper%residual     => pf_generic_residual

    sweeper%sweeperctx = c_loc(exp)
  end subroutine pf_explicit_create

  subroutine pf_explicit_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_explicit_t), pointer :: explicit
    call c_f_pointer(sweeper%sweeperctx, explicit)
    deallocate(explicit%SdiffE)
    deallocate(explicit)
  end subroutine pf_explicit_destroy

end module pf_mod_explicit

