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
  use pf_mod_utils
  implicit none
  integer, parameter, private :: npieces = 1

  interface
     subroutine pf_f2eval_p(y, t, level, levelctx, f2)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, f2, levelctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2eval_p
  end interface

  interface
     subroutine pf_f2comp_p(y, t, dt, rhs, level, levelctx, f2)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, rhs, f2, levelctx
       real(pfdp),     intent(in)        :: t, dt
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2comp_p
  end interface

  type :: pf_implicit_t
     procedure(pf_f2eval_p),  pointer, nopass :: f2eval
     procedure(pf_f2comp_p),  pointer, nopass :: f2comp

     real(pfdp), ALLOCATABLE :: SdiffI(:,:)
  end type pf_implicit_t

contains

  ! Perform one SDC sweep on level Lev and set qend appropriately.
  subroutine implicit_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:Lev%nnodes-1)
    type(c_ptr) :: rhs

    type(pf_implicit_t), pointer :: imp

    call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0_pfdp)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*imp%SdiffI(m,n), Lev%F(n,1))
       end do
       if (associated(Lev%tau)) then
          call Lev%encap%axpy(Lev%S(m), 1.0_pfdp, Lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%unpack(Lev%Q(1), Lev%q0)

    call imp%f2eval(Lev%Q(1), t0, Lev%level, Lev%levelctx, Lev%F(1,1))

    call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%levelctx, Lev%encap%encapctx)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%copy(rhs, Lev%Q(m))
       call Lev%encap%axpy(rhs, 1.0_pfdp, Lev%S(m))

       call imp%f2comp(Lev%Q(m+1), t, dtsdc(m), rhs, Lev%level, Lev%levelctx, Lev%F(m+1,1))
    end do

    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine implicit_sweep

  ! Evaluate function values
  subroutine implicit_evaluate(Lev, t, m)
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_implicit_t), pointer :: imp
    call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    call imp%f2eval(Lev%Q(m), t, Lev%level, Lev%levelctx, Lev%F(m,1))
  end subroutine implicit_evaluate

  ! Initialize matrix
  subroutine implicit_initialize(Lev)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)

    integer :: m,nnodes
    type(pf_implicit_t), pointer :: imp
    call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    nnodes = Lev%nnodes
    allocate(imp%SdiffI(nnodes-1,nnodes))  !  S-BE

    imp%SdiffI = Lev%s0mat

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       imp%SdiffI(m,m+1) = imp%SdiffI(m,m+1) - dsdc(m)
    end do

  end subroutine implicit_initialize


  ! Compute SDC integral
  subroutine implicit_integrate(Lev, qSDC, fSDC, dt, fintSDC)
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
  end subroutine implicit_integrate

  ! Create implicit sweeper
  subroutine pf_implicit_create(sweeper, f2eval, f2comp)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp

    type(pf_implicit_t), pointer :: imp

    allocate(imp)
    imp%f2eval => f2eval
    imp%f2comp => f2comp

    sweeper%npieces = npieces
    sweeper%sweep      => implicit_sweep
    sweeper%evaluate   => implicit_evaluate
    sweeper%initialize => implicit_initialize
    sweeper%integrate  => implicit_integrate
    sweeper%evaluate_all => pf_generic_evaluate_all
    sweeper%residual     => pf_generic_residual

    sweeper%sweeperctx = c_loc(imp)
  end subroutine pf_implicit_create

  subroutine pf_implicit_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_implicit_t), pointer :: imp
    call c_f_pointer(sweeper%sweeperctx, imp)
    deallocate(imp%SdiffI)
    deallocate(imp)
  end subroutine pf_implicit_destroy

end module pf_mod_implicit

