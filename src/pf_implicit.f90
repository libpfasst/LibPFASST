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

  type, extends(pf_sweeper_t), abstract :: pf_implicit_t
     real(pfdp), allocatable :: SdiffI(:,:)
   contains
     procedure(pf_f2eval_p), deferred :: f2eval
     procedure(pf_f2comp_p), deferred :: f2comp
     procedure :: sweep        => implicit_sweep
     procedure :: initialize   => implicit_initialize
     procedure :: evaluate     => implicit_evaluate
     procedure :: integrate    => implicit_integrate
     procedure :: residual     => implicit_residual
     procedure :: evaluate_all => implicit_evaluate_all
     procedure :: destroy      => implicit_destroy
     procedure :: implicit_destroy
  end type pf_implicit_t

  interface
     subroutine pf_f2eval_p(this, y, t, level, f2)
       import pf_implicit_t, pf_encap_t, c_int, pfdp
       class(pf_implicit_t), intent(inout) :: this
       class(pf_encap_t),    intent(in   ) :: y
       class(pf_encap_t),    intent(inout) :: f2
       real(pfdp),           intent(in   ) :: t
       integer(c_int),       intent(in   ) :: level
     end subroutine pf_f2eval_p

     subroutine pf_f2comp_p(this, y, t, dt, rhs, level, f2)
       import pf_implicit_t, pf_encap_t, c_int, pfdp
       class(pf_implicit_t), intent(inout) :: this
       class(pf_encap_t),    intent(in   ) :: rhs
       class(pf_encap_t),    intent(inout) :: y, f2
       real(pfdp),           intent(in   ) :: t, dt
       integer(c_int),       intent(in   ) :: level
     end subroutine pf_f2comp_p
  end interface

contains

  ! Perform one SDC sweep on level lev and set qend appropriately.
  subroutine implicit_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer
    class(pf_implicit_t), intent(inout) :: this
    type(pf_pfasst_t),    intent(inout) :: pf
    real(pfdp),           intent(in   ) :: dt, t0
    class(pf_level_t),    intent(inout) :: lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:lev%nnodes-1)
    class(pf_encap_t), allocatable :: rhs

    call start_timer(pf, TLEVEL+lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1
       call lev%S(m)%setval(0.0_pfdp)
       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%SdiffI(m,n), lev%F(n,1))
       end do
       if (allocated(lev%tau)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%unpack(lev%q0)

    call this%f2eval(lev%Q(1), t0, lev%level, lev%F(1,1))

    call lev%ulevel%factory%create0(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       call rhs%copy(lev%Q(m))
       call rhs%axpy(1.0_pfdp, lev%S(m))

       call this%f2comp(lev%Q(m+1), t, dtsdc(m), rhs, lev%level, lev%F(m+1,1))
    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    call lev%ulevel%factory%destroy0(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call end_timer(pf, TLEVEL+lev%level-1)
  end subroutine implicit_sweep

  ! Evaluate function values
  subroutine implicit_evaluate(this, lev, t, m)
    class(pf_implicit_t), intent(inout) :: this
    real(pfdp),           intent(in   ) :: t
    integer,              intent(in   ) :: m
    class(pf_level_t),    intent(inout) :: lev

    call this%f2eval(lev%Q(m), t, lev%level, lev%F(m,1))
  end subroutine implicit_evaluate

  ! Initialize matrix
  subroutine implicit_initialize(this, lev)
    use pf_mod_dtype
    class(pf_implicit_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m,nnodes

    this%npieces = 1
    
    nnodes = lev%nnodes
    allocate(this%SdiffI(nnodes-1,nnodes))  !  S-BE

    this%SdiffI = lev%s0mat

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       this%SdiffI(m,m+1) = this%SdiffI(m,m+1) - dsdc(m)
    end do

  end subroutine implicit_initialize

  ! Destroy the matrices
  subroutine implicit_destroy(this, lev)
    class(pf_implicit_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%SdiffI)
  end subroutine implicit_destroy


  ! Compute SDC integral
  subroutine implicit_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_implicit_t), intent(inout) :: this
    class(pf_level_t),    intent(in   ) :: lev
    class(pf_encap_t),    intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),           intent(in   ) :: dt
    class(pf_encap_t),    intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine implicit_integrate

  subroutine implicit_residual(this, lev, dt)
    class(pf_implicit_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev
    real(pfdp),           intent(in   ) :: dt
    call pf_generic_residual(this, lev, dt)
  end subroutine implicit_residual

  subroutine implicit_evaluate_all(this, lev, t)
    class(pf_implicit_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev
    real(pfdp),           intent(in   ) :: t(:)
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine implicit_evaluate_all

end module pf_mod_implicit
