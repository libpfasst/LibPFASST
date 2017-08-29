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

  type, extends(pf_sweeper_t), abstract :: pf_explicit_t
     real(pfdp), allocatable :: SdiffE(:,:)
   contains
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure :: sweep        => explicit_sweep
     procedure :: initialize   => explicit_initialize
     procedure :: evaluate     => explicit_evaluate
     procedure :: integrate    => explicit_integrate
     procedure :: residual     => explicit_residual
     procedure :: evaluate_all => explicit_evaluate_all
     procedure :: destroy      => explicit_destroy
     procedure :: explicit_destroy
  end type pf_explicit_t

  interface
    subroutine pf_f_eval_p(this,y, t, level, f)
      import pf_explicit_t, pf_encap_t, c_int, pfdp
      class(pf_explicit_t),  intent(inout) :: this
      class(pf_encap_t), intent(in   ) :: y
      real(pfdp),        intent(in   ) :: t
      integer(c_int),    intent(in   ) :: level
      class(pf_encap_t), intent(inout) :: f
    end subroutine pf_f_eval_p
  end interface
contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine explicit_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer
    class(pf_explicit_t), intent(inout) :: this
    type(pf_pfasst_t),    intent(inout) :: pf
    real(pfdp),           intent(in   ) :: dt, t0
    class(pf_level_t),    intent(inout) :: lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:lev%nnodes-1)

    call start_timer(pf, TLEVEL+lev%index-1)

    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1
       call lev%S(m)%setval(0.0_pfdp)
       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%SdiffE(m,n), lev%F(n,1))
       end do
       if (allocated(lev%tau)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%copy(lev%q0)

    call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1))

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       call lev%Q(m+1)%copy(lev%Q(m))
       call lev%Q(m+1)%axpy(dtsdc(m), lev%F(m,1))
       call lev%Q(m+1)%axpy(1.0_pfdp, lev%S(m))

       call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1))
    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    ! done
    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine explicit_sweep

  ! Evaluate function values
  subroutine explicit_evaluate(this, lev, t, m)
    use pf_mod_dtype
    class(pf_explicit_t), intent(inout) :: this
    real(pfdp),           intent(in   ) :: t
    integer,              intent(in   ) :: m
    class(pf_level_t),    intent(inout) :: lev

    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1))
  end subroutine explicit_evaluate

  ! Initialize matrix
  subroutine explicit_initialize(this, lev)
    use pf_mod_dtype
    class(pf_explicit_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m,nnodes

    this%npieces = 1

    nnodes = lev%nnodes
    allocate(this%SdiffE(nnodes-1,nnodes))  ! S-FE

    this%SdiffE = lev%s0mat

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       this%SdiffE(m,m)   = this%SdiffE(m,m)   - dsdc(m)
    end do
  end subroutine explicit_initialize

  ! Destroy the matrices
  subroutine explicit_destroy(this, lev)
    class(pf_explicit_t), intent(inout) :: this
    class(pf_level_t), intent(inout)    :: lev
    
    deallocate(this%SdiffE)
  end subroutine explicit_destroy

  ! Compute SDC integral
  subroutine explicit_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_explicit_t), intent(inout) :: this
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
  end subroutine explicit_integrate

  subroutine explicit_residual(this, lev, dt)
    class(pf_explicit_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev
    real(pfdp),           intent(in   ) :: dt
    call pf_generic_residual(this, lev, dt)
  end subroutine explicit_residual

  subroutine explicit_evaluate_all(this, lev, t)
    class(pf_explicit_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev
    real(pfdp),           intent(in   ) :: t(:)
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine explicit_evaluate_all

end module pf_mod_explicit
