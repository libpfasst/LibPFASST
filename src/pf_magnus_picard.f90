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
     real(pfdp), allocatable :: dtsdc(:)
     integer :: magnus_order, qtype
     real(pfdp) :: dt
     real(pfdp), allocatable :: commutator_coefs(:,:,:)
     complex(pfdp), allocatable :: commutators(:,:,:)
     class(pf_encap_t), allocatable :: omega(:), time_ev_op(:)
   contains
     procedure :: sweep      => magpicard_sweep
     procedure :: initialize => magpicard_initialize
     procedure :: evaluate   => magpicard_evaluate
     procedure :: integrate  => magpicard_integrate
     procedure :: residual   => magpicard_residual
     procedure :: evaluate_all => magpicard_evaluate_all
     procedure(pf_f_eval_p), deferred :: f_eval
     ! procedure(pf_compute_single_commutators_p), deferred :: compute_single_commutators
     procedure(pf_compute_omega_p), deferred :: compute_omega
     procedure(pf_compute_time_ev_ops_p), deferred :: compute_time_ev_ops
     procedure(pf_propagate_solution_p), deferred :: propagate_solution
     procedure(pf_destroy_magpicard_p), deferred :: destroy
     procedure :: magpicard_destroy
  end type pf_magpicard_t

  interface
     subroutine pf_f_eval_p(this, y, t, level, f)
       import pf_magpicard_t, pf_encap_t, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y
       real(pfdp),        intent(in   ) :: t
       integer,           intent(in   ) :: level
       class(pf_encap_t), intent(inout) :: f
     end subroutine pf_f_eval_p
     ! subroutine pf_compute_single_commutators_p(this, f)
     !   import pf_magpicard_t, pf_encap_t, pfdp
     !   class(pf_magpicard_t),  intent(inout) :: this
     !   class(pf_encap_t), intent(inout) :: f(:,:)
     ! end subroutine pf_compute_single_commutators_p
     subroutine pf_compute_omega_p(this, omega, integrals, f, nodes, qmat, dt, this_node, coefs)
       import pf_magpicard_t, pf_encap_t, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega
       class(pf_encap_t), intent(inout) :: f(:,:), integrals(:)
       real(pfdp), intent(in) :: coefs(:,:), nodes(:), qmat(:,:), dt
       integer, intent(in) :: this_node
     end subroutine pf_compute_omega_p
     subroutine pf_compute_time_ev_ops_p(this, time_ev_op, omega, level)
       import pf_magpicard_t, pf_encap_t, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega
       class(pf_encap_t), intent(inout) :: time_ev_op
       integer, intent(in) :: level
     end subroutine pf_compute_time_ev_ops_p
     subroutine pf_propagate_solution_p(this, sol_t0, sol_tn, u)
       import pf_magpicard_t, pf_encap_t, pfdp
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
  subroutine magpicard_sweep(this, pf, level_index, t0, dt, nsweeps)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_magpicard_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    real(pfdp), intent(in) :: dt, t0
    integer,             intent(in)    :: level_index
    integer,             intent(in)    :: nsweeps

    class(pf_level_t), pointer :: lev
    integer    :: m, nnodes, k

    real(pfdp) :: t

    lev => pf%levels(level_index)
    nnodes = lev%nnodes

    call call_hooks(pf, level_index, PF_PRE_SWEEP)
    call lev%Q(1)%copy(lev%q0)

    call start_timer(pf, TLEVEL+lev%index-1)
    do k = 1,nsweeps
       ! Copy values into residual
       do m = 1, nnodes-1
          call lev%R(m)%copy(lev%Q(m+1))
       end do

       if (k .eq. 1) then
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(m,1))
       end if

       t = t0
       do m = 2, nnodes
          t=t+ dt*this%dtsdc(m-1)
          call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1))
       end do

       ! if (this%magnus_order > 1) then
       !    call start_timer(pf, TAUX+2)
       !    call this%compute_single_commutators(lev%F)
       !    call end_timer(pf, TAUX+2)
       ! endif

       call magpicard_integrate(this, lev, lev%Q, lev%F, dt, lev%I)

       do m = 1, nnodes-1
          call start_timer(pf, TAUX)
          call this%compute_omega(this%omega(m), lev%I, lev%F, &
               lev%nodes, lev%qmat, dt, m, this%commutator_coefs(:,:,m))
          call end_timer(pf, TAUX)

          call start_timer(pf, TAUX+1)
          call this%compute_time_ev_ops(this%time_ev_op(m), this%omega(m), lev%index)
          call end_timer(pf, TAUX+1)

          call this%propagate_solution(lev%Q(1), lev%Q(m+1), this%time_ev_op(m))
       end do

       call pf_residual(pf, lev, dt)
       call call_hooks(pf, level_index, PF_POST_SWEEP)

    end do  ! Loop over sweeps

    call lev%qend%copy(lev%Q(nnodes))
    call end_timer(pf, TLEVEL+lev%index-1)

  end subroutine magpicard_sweep

  subroutine magpicard_initialize(this, lev)
    class(pf_magpicard_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer :: m, nnodes

    this%npieces = 1
    nnodes = lev%nnodes

    allocate(this%dtsdc(nnodes-1))
    this%dtsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)   !  SDC time step size (unscaled)

    call get_commutator_coefs(this%qtype, nnodes, this%dt, this%commutator_coefs)

    call lev%ulevel%factory%create_array(this%omega, nnodes-1, &
         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call lev%ulevel%factory%create_array(this%time_ev_op, nnodes-1, &
         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    do m = 1, nnodes-1
        call this%omega(m)%setval(0.0_pfdp)
        call this%time_ev_op(m)%setval(0.0_pfdp)
    end do

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
          call fintSDC(m)%axpy(dt*lev%qmat(m,j), fSDC(j,1))
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

    do m = 1, lev%nnodes-1
       call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
    end do

  end subroutine magpicard_residual

  ! Destroy the matrices
  subroutine magpicard_destroy(this, lev)
      class(pf_magpicard_t),  intent(inout) :: this
      class(pf_level_t),    intent(inout) :: lev

      deallocate(this%dtsdc, this%commutator_coefs, this%commutators)

      call lev%ulevel%factory%destroy_array(this%omega, lev%nnodes-1, &
           lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
      call lev%ulevel%factory%destroy_array(this%time_ev_op, lev%nnodes-1, &
           lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

  end subroutine magpicard_destroy

  subroutine get_commutator_coefs(qtype, nnodes, dt, coefs)
    integer, intent(in) :: qtype, nnodes
    real(pfdp), intent(in) :: dt
    real(pfdp), intent(inout) :: coefs(:,:,:)

    if (qtype == 1) then
       ! we're talking Lobatto nodes, where nnodes=3 includes, t0, t_1/2, tn
       ! need some way to differentiate whether you want full collocation or not
        if (nnodes > 2) then
            ! coefs(1:3, 1, 1) = dt**2 * [real(8)::11/480., -1/480., 1/480.]
            ! coefs(1:3, 1, 2) = dt**2 * [real(8)::1/15., 1/60., 1/15.]
            coefs(1, 1, 1) = -1/48.d0 * dt**2
            coefs(2, 1, 2) = -1/12.d0 * dt**2
        endif
    else if (qtype == 5) then
      if (nnodes >= 3) then
          coefs(1:3, 1, 1) = 1.d-3 * [real(8)::-0.708256232441739d0, 0.201427439334681d0, -0.002608155816283d0]
          coefs(1:3, 1, 2) = [real(8)::-0.035291589565775d0, 0.004482619613666d0, -0.000569367343553d0]
          coefs(1:3, 1, 3) = [real(8)::-0.078891497044705d0, -0.018131905893999d0, -0.035152700676886d0]
          coefs(1:3, 1, 4) = [real(8)::-0.071721913818656d0, -0.035860956909328d0, -0.071721913818656d0]
          coefs(:,1,:) = dt**2 * coefs(:,1,:)

          coefs(:, 2, 1) = &
               [real(8)::1.466782892818107d-6, -2.546845448743404d-6, 7.18855795894042d-7, &
                -3.065370250683271d-7, 6.962336322868984d-7, -1.96845581200288d-7,  &
                -2.262216360714434d-8, -2.72797194008496d-9, 8.54843541920492d-10]
          coefs(:, 2, 2) = &
               [real(8) ::0.001040114336531742d0, -0.001714330280871491d0, 0.0001980882752518163d0, &
                -0.00006910549596945875d0, 0.0002905401601450182d0, -0.00003465884693947625d0, &
                 0.0000924518848932026d0, 0.0000125950571649574d0, -2.4709074423913880d-6]
          coefs(:, 2, 3) = &
               [real(8)::0.004148295975360902d0, -0.006387421893168941d0, -0.003594231910817328d0, &
                 0.000997378110327084d0, 0.0001241530237557625d0, -0.0003805975423160699d0, &
                 0.003718384934573079d0, 0.001693514295056844d0, -0.001060408584538103d0]
          coefs(:, 2, 4) = &
               [real(8)::0.003453850676072909d0, -0.005584950029394391d0, -0.007128159905937654d0, &
                 0.001653439153439147d0, 0.0d0, -0.001653439153439143d0, &
                 0.007128159905937675d0, 0.005584950029394475d0, -0.003453850676072897d0]
          coefs(:,2,:) = dt**3 * coefs(:,2,:)

          coefs(1, 3, 4) = dt**4 / 60.d0
      endif
    else
      stop 'oh no! unsupported qtype'
    endif
  end subroutine get_commutator_coefs
end module pf_mod_magnus_picard
