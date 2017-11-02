!
! Copyright (C) 2017 Brandon Krull and Michael Minion.
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

module pf_mod_imk
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_imk_t
    class(pf_encap_t), allocatable :: omega(:), time_ev_op(:), dexpinvs(:)
  contains
    procedure :: sweep        => imk_sweep
    procedure :: initialize   => imk_initialize
    procedure :: evaluate     => imk_evaluate
    procedure :: integrate    => imk_integrate
    procedure :: residual     => imk_residual
    procedure :: evaluate_all => imk_evaluate_all
    procedure :: destroy      => imk_destroy
    procedure(pf_f_eval_p), deferred :: f_eval
    procedure(pf_f_comp_p), deferred :: f_comp
    procedure(pf_dexpinv_p), deferred :: dexpinv
    procedure(pf_propagate_p), deferred :: propagate
  end type pf_imk_t

  interface
     subroutine pf_f_eval_p(this,y, t, level, f)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t),   intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y, f
       real(pfdp),        intent(in   ) :: t
       integer(c_int),    intent(in   ) :: level
     end subroutine pf_f_eval_p
     subroutine pf_f_comp_p(this,y, t, dt, level, f)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t),   intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y, f
       real(pfdp),        intent(in   ) :: t, dt
       integer(c_int),    intent(in   ) :: level
     end subroutine pf_f_comp_p
     subroutine pf_dexpinv_p(this, omega, a)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega, a
     end subroutine pf_dexpinv_p
     subroutine pf_propagate_p(this, q, time_ev_op, q0)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: q, time_ev_op, q0
     end subroutine pf_propagate_p
  end interface

contains

  subroutine imk_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer

    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: dt, t0

    integer :: m
    real(pfdp) :: t, dtsdc(lev%nnodes-1)

    call start_timer(pf, TLEVEL+lev%index-1)

    do m = 1, lev%nnodes-1
       call lev%R(m)%copy(lev%Q(m+1))
    end do

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes
       call lev%Q(m)%copy(lev%q0)

       ! compute new Q from previous omegas
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1))
       this%dexpinvs(m) = dexpinv(omega, lev%F(m,1))
       t=t+dtsdc(m)
    end do

    call imk_integrate(this, lev, lev%Q, lev%F, dt, lev%S)

    do m = 1, lev%nnodes
       call this%propagate()
    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    call end_timer(pf, TLEVEL+lev%index-1)

  end subroutine imk_sweep

  subroutine imk_initialize(this, lev)
    class(pf_imk_t),   intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer :: m, nnodes

    nnodes = lev%nnodes
    call lev%ulevel%factory%create_array(this%omega, nnodes-1, &
         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call lev%ulevel%factory%create_array(this%time_ev_op, nnodes-1, &
         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call lev%ulevel%factory%create_array(this%dexpinvs, nnodes-1, &
         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    do m = 1, nnodes-1
       call this%omega(m)%setval(0.0_pfdp)
       call this%time_ev_op(m)%setval(0.0_pfdp)
    end do

  end subroutine imk_initialize

  subroutine imk_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_imk_t),   intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in   ) :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)

    integer :: j, m, p

    do m = 1, lev%nnodes-1
       call fintSDC(m)%setval(0.0_pfdp)
       do j = 1, lev%nnodes
          !do p = 1, this%npieces
          call fintSDC(m)%axpy(dt*lev%qmat(m,j), fSDC(j,1))
          !end do
       end do
    end do
  end subroutine imk_integrate

  function dexpinv(omega, a) result(x)
    integer, intent(in) :: dim
    complex(pfdp), intent(in) :: omega(dim,dim), a(dim,dim)

    integer :: i, nterms
    real(pfdp) :: bernoulli(10), coef
    complex(pfdp), allocatable :: x(:,:), commutator(:,:)

    bernoulli = (/-0.5, 0.16666666666666, 0.0, -0.033333333333, 0.0, &
         0.023809523809523808, 0.0, -0.03333333333333, 0.0, 0.075757575757575/)

    allocate(x, commutator, source=a)

    nterms = size(bernoulli)
    coef = 1.0
    do i = 1, nterms
       commutator = matmul(a, commutator)
       coef = coef/i
       x = x + bernoulli(i) * coef * commutator
    end do

  end function dexpinv

  subroutine imk_evaluate(this, lev, t, m)
    use pf_mod_dtype
    class(pf_imk_t),   intent(inout) :: this
    real(pfdp),        intent(in   ) :: t
    integer,           intent(in   ) :: m
    class(pf_level_t), intent(inout) :: lev

    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1))
  end subroutine imk_evaluate

  subroutine imk_evaluate_all(this, lev, t)
    class(pf_imk_t),   intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t(:)
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine imk_evaluate_all

  subroutine imk_residual(this, lev, dt)
    class(pf_imk_t),   intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: dt
    integer :: m
    ! call pf_generic_residual(this, lev, dt)

    do m = 1, lev%nnodes-1
       call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
    end do
  end subroutine imk_residual

  subroutine imk_destroy(this, lev)
      class(pf_imk_t),   intent(inout) :: this
      class(pf_level_t), intent(inout) :: lev

      ! deallocate(this%commutator_colloc_coefs)

      call lev%ulevel%factory%destroy_array(this%omega, lev%nnodes-1, &
           lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
      call lev%ulevel%factory%destroy_array(this%time_ev_op, lev%nnodes-1, &
           lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
      call lev%ulevel%factory%destroy_array(this%dexpinvs, lev%nnodes-1, &
           lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

  end subroutine imk_destroy

end module pf_mod_imk
