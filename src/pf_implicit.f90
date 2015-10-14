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

  type, extends(pf_sweeper_t), abstract :: pf_implicit_t
     real(pfdp), allocatable :: SdiffI(:,:)
   contains
     procedure(pf_f2eval_p), deferred :: f2eval
     procedure(pf_f2comp_p), deferred :: f2comp
     procedure :: sweep => implicit_sweep
     procedure :: initialize => implicit_initialize
     procedure :: evaluate => implicit_evaluate
     procedure :: integrate => implicit_integrate
     procedure :: residual => implicit_residual
     procedure :: evaluate_all => implicit_evaluate_all
  end type pf_implicit_t

  interface
     subroutine pf_f2eval_p(this, y, t, level, levelctx, f2)
       import pf_implicit_t, c_ptr, c_int, pfdp
       class(pf_implicit_t), intent(in) :: this
       type(c_ptr),    intent(in), value :: y, f2, levelctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2eval_p

     subroutine pf_f2comp_p(this, y, t, dt, rhs, level, levelctx, f2)
       import pf_implicit_t, c_ptr, c_int, pfdp
       class(pf_implicit_t), intent(in) :: this
       type(c_ptr),    intent(in), value :: y, rhs, f2, levelctx
       real(pfdp),     intent(in)        :: t, dt
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2comp_p
  end interface

contains

  ! Perform one SDC sweep on level Lev and set qend appropriately.
  subroutine implicit_sweep(this, pf, Lev, t0, dt)
    use pf_mod_timer

    class(pf_implicit_t), intent(in) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:Lev%nnodes-1)
    type(c_ptr) :: rhs

    ! type(pf_implicit_t), pointer :: imp

    ! call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0_pfdp)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*this%SdiffI(m,n), Lev%F(n,1))
       end do
       if (associated(Lev%tau)) then
          call Lev%encap%axpy(Lev%S(m), 1.0_pfdp, Lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%unpack(Lev%Q(1), Lev%q0)

    call this%f2eval(Lev%Q(1), t0, Lev%level, Lev%levelctx, Lev%F(1,1))

    call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%levelctx, Lev%encap%encapctx)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%copy(rhs, Lev%Q(m))
       call Lev%encap%axpy(rhs, 1.0_pfdp, Lev%S(m))

       call this%f2comp(Lev%Q(m+1), t, dtsdc(m), rhs, Lev%level, Lev%levelctx, Lev%F(m+1,1))
    end do

    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine implicit_sweep

  ! Evaluate function values
  subroutine implicit_evaluate(this, Lev, t, m)
    class(pf_implicit_t), intent(in) :: this
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: Lev

    ! type(pf_implicit_t), pointer :: imp
    ! call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    call this%f2eval(Lev%Q(m), t, Lev%level, Lev%levelctx, Lev%F(m,1))
  end subroutine implicit_evaluate

  ! Initialize matrix
  subroutine implicit_initialize(this, Lev)
    use pf_mod_dtype
    class(pf_implicit_t), intent(inout) :: this
    
    type(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)

    integer :: m,nnodes
    ! type(pf_implicit_t), pointer :: imp
    ! call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    nnodes = Lev%nnodes
    allocate(this%SdiffI(nnodes-1,nnodes))  !  S-BE

    this%SdiffI = Lev%s0mat

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       this%SdiffI(m,m+1) = this%SdiffI(m,m+1) - dsdc(m)
    end do

  end subroutine implicit_initialize


  ! Compute SDC integral
  subroutine implicit_integrate(this, Lev, qSDC, fSDC, dt, fintSDC)
    class(pf_implicit_t), intent(in) :: this
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

  subroutine implicit_residual(this, Lev, dt)
    class(pf_implicit_t), intent(in)  :: this
    type(pf_level_t),  intent(inout) :: Lev
    real(pfdp),        intent(in)    :: dt

    integer :: m, n

    call pf_generic_residual(this, Lev, dt)
  end subroutine implicit_residual

  subroutine implicit_evaluate_all(this, Lev, t)
    class(pf_implicit_t), intent(in)  :: this
    type(pf_level_t),  intent(inout) :: Lev
    real(pfdp),        intent(in)    :: t(:)

    call pf_generic_evaluate_all(this, Lev, t)
  end subroutine implicit_evaluate_all

end module pf_mod_implicit

