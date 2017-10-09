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

module pf_mod_magnus_picard
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_magpicard_t
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
     logical                 :: use_LUq_ = .false.
     class(pf_encap_t), allocatable :: omega(:), time_ev_op(:)
     real(pfdp), allocatable :: commutator_colloc_coefs(:,:)
   contains
     procedure :: sweep      => magpicard_sweep
     procedure :: initialize => magpicard_initialize
     procedure :: evaluate   => magpicard_evaluate
     procedure :: integrate  => magpicard_integrate
     procedure :: residual   => magpicard_residual
     procedure :: evaluate_all => magpicard_evaluate_all
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_compute_omega_p), deferred :: compute_omega
     procedure(pf_compute_time_ev_ops_p), deferred :: compute_time_ev_ops
     procedure(pf_propagate_solution_p), deferred :: propagate_solution
     procedure(pf_destroy_magpicard_p), deferred :: destroy
     procedure :: magpicard_destroy
  end type pf_magpicard_t

  interface
     subroutine pf_f_eval_p(this, y, t, level, f)
       import pf_magpicard_t, pf_encap_t, c_int, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y
       real(pfdp),        intent(in   ) :: t
       integer(c_int),    intent(in   ) :: level
       class(pf_encap_t), intent(inout) :: f
     end subroutine pf_f_eval_p
     subroutine pf_compute_omega_p(this, omega, integrals, f, coef)
       import pf_magpicard_t, pf_encap_t, c_int, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega, integrals
       class(pf_encap_t), intent(inout) :: f(:)
       real(pfdp), intent(in) :: coef(:)
     end subroutine pf_compute_omega_p
     subroutine pf_compute_time_ev_ops_p(this, omega, u, level)
       import pf_magpicard_t, pf_encap_t, c_int, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega
       class(pf_encap_t), intent(inout) :: u
       integer, intent(in) :: level
     end subroutine pf_compute_time_ev_ops_p
     subroutine pf_propagate_solution_p(this, sol_t0, sol_tn, u)
       import pf_magpicard_t, pf_encap_t, c_int, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: sol_t0
       class(pf_encap_t), intent(inout) :: u
       class(pf_encap_t), intent(inout) :: sol_tn
     end subroutine pf_propagate_solution_p
     subroutine pf_destroy_magpicard_p(this, Lev)
       import pf_magpicard_t, pf_level_t
       class(pf_magpicard_t), intent(inout) :: this
       class(pf_level_t), intent(inout) :: Lev
     end subroutine pf_destroy_magpicard_p
  end interface
contains

  ! Perform one SDC sweep on level Lev and set qend appropriately.
  subroutine magpicard_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer

    class(pf_magpicard_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp), intent(in) :: dt, t0
    class(pf_level_t), intent(inout) :: lev

    integer    :: m, n, p
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:lev%nnodes-1)

    call start_timer(pf, TLEVEL+lev%index-1)

    do m = 1, lev%nnodes-1
       call lev%R(m)%copy(lev%Q(m+1))
    end do
    call lev%Q(1)%copy(lev%q0)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1))
       t=t+dtsdc(m)
    end do

    call magpicard_integrate(this, lev, lev%Q, lev%F, dt, lev%S)

    if (lev%nnodes == 3) then
        this%commutator_colloc_coefs(1,:) = (/-11/480., 1/480., -1/480./)
        this%commutator_colloc_coefs(2,:) = (/-1/15., -1/60., -1/15./)
     else
        this%commutator_colloc_coefs(:,:) = 0.0_pfdp
     endif

    do m = 1, lev%nnodes-1
       call start_timer(pf, TAUX)
       call this%compute_omega(this%omega(m), lev%S(m), lev%F(:,1), &
            dt**2 * this%commutator_colloc_coefs(m,:))
       call end_timer(pf, TAUX)

       call start_timer(pf, TAUX+1)
       call this%compute_time_ev_ops(this%omega(m), this%time_ev_op(m), lev%index)
       call end_timer(pf, TAUX+1)

       call start_timer(pf, TAUX+2)
       call this%propagate_solution(lev%Q(1), lev%Q(m+1), this%time_ev_op(m))
       call end_timer(pf, TAUX+2)
    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine magpicard_sweep

  subroutine magpicard_initialize(this, lev)
    class(pf_magpicard_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer :: m,n,nnodes

    nnodes = lev%nnodes
    allocate(this%commutator_colloc_coefs(nnodes-1, nnodes))
    allocate(this%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(this%QtilI(nnodes-1,nnodes))  !  S-BE
    this%QtilI = 0.0_pfdp

    call lev%ulevel%factory%create_array(this%omega, &
          nnodes-1, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call lev%ulevel%factory%create_array(this%time_ev_op, &
         nnodes-1, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    if(this%use_LUq_) then
        call myLUq(lev%qmat, this%QtilI, Nnodes, 1)
    else
        dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
        do m = 1, nnodes-1
            do n = 1,m
                this%QtilI(m,n+1) =  dsdc(n)
            end do
            call this%omega(m)%setval(0.0_pfdp)
            call this%time_ev_op(m)%setval(0.0_pfdp)
        end do
    end if

    this%commutator_colloc_coefs(:,:) = 0.0_pfdp
    this%QdiffI = lev%qmat-this%QtilI

  end subroutine magpicard_initialize

  !> Compute SDC integral
  !>  fintSDC = \int_{t_n}^{t_m} fSDC dt
  subroutine magpicard_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_magpicard_t), intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in   ) :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)

    integer :: j, m, p

    do m = 1, lev%nnodes-1
       call fintSDC(m)%setval(0.0_pfdp)
       do j = 1, lev%nnodes
          ! print*, 'qmat=', lev%qmat(m,j)
          do p = 1, this%npieces
             call fintSDC(m)%axpy(dt*lev%qmat(m,j), fSDC(j,1))
          end do
       end do
    end do
  end subroutine magpicard_integrate

  ! Evaluate function values
  subroutine magpicard_evaluate(this, lev, t, m)
    use pf_mod_dtype
    class(pf_magpicard_t), intent(inout) :: this
    real(pfdp),           intent(in   ) :: t
    integer,              intent(in   ) :: m
    class(pf_level_t),    intent(inout) :: lev

    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1))
  end subroutine magpicard_evaluate

  subroutine magpicard_evaluate_all(this, lev, t)
    class(pf_magpicard_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev
    real(pfdp),           intent(in   ) :: t(:)
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine magpicard_evaluate_all

  subroutine magpicard_residual(this, lev, dt)
    class(pf_magpicard_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev
    real(pfdp),           intent(in   ) :: dt
    integer :: m
    ! call pf_generic_residual(this, lev, dt)

    do m = 1, lev%nnodes-1
       call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
    end do
  end subroutine magpicard_residual

  ! Destroy the matrices
  subroutine magpicard_destroy(this, lev)
      class(pf_magpicard_t),  intent(inout) :: this
      class(pf_level_t),    intent(inout) :: lev
      integer :: i, j

      deallocate(this%QdiffI)
      deallocate(this%QtilI)
      deallocate(this%commutator_colloc_coefs)

      call lev%ulevel%factory%destroy_array(this%omega, &
           lev%nnodes-1, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
      call lev%ulevel%factory%destroy_array(this%time_ev_op, &
           lev%nnodes-1, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

  end subroutine magpicard_destroy
end module pf_mod_magnus_picard
