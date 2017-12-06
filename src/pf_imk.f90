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
    class(pf_encap_t), allocatable :: omega(:), dexpinvs(:)
  contains
    procedure :: sweep        => imk_sweep
    procedure :: initialize   => imk_initialize
    procedure :: evaluate     => imk_evaluate
    procedure :: integrate    => imk_integrate
    procedure :: residual     => imk_residual
    procedure :: evaluate_all => imk_evaluate_all
    procedure :: imk_destroy
    procedure(pf_f_eval_p), deferred :: f_eval
    ! procedure(pf_f_comp_p), deferred :: f_comp
    procedure(pf_dexpinv_p), deferred :: dexpinv
    procedure(pf_propagate_p), deferred :: propagate
 end type pf_imk_t

  interface
     subroutine pf_f_eval_p(this, y, t, level, f)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t),   intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y, f
       real(pfdp),        intent(in   ) :: t
       integer(c_int),    intent(in   ) :: level
     end subroutine pf_f_eval_p
     ! subroutine pf_f_comp_p(this,y, t, dt, level, f)
     !   import pf_imk_t, pf_encap_t, c_int, pfdp
     !   class(pf_imk_t),   intent(inout) :: this
     !   class(pf_encap_t), intent(inout) :: y, f
     !   real(pfdp),        intent(in   ) :: t, dt
     !   integer(c_int),    intent(in   ) :: level
     ! end subroutine pf_f_comp_p
     subroutine pf_dexpinv_p(this, omega, a, x)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega, a, x
     end subroutine pf_dexpinv_p
     subroutine pf_propagate_p(this, q, omega, q0)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: q, q0, omega
     end subroutine pf_propagate_p
  end interface

contains

  subroutine imk_sweep(this, pf, level_index, t0, dt,nsweeps)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    integer,             intent(in)    :: level_index
    real(pfdp),        intent(in   ) :: dt
    real(pfdp),        intent(in   ) :: t0
    integer,             intent(in)    :: nsweeps

    integer :: m,k
    real(pfdp) :: t
    real(pfdp)  :: dtsdc(1:pf%levels(level_index)%nnodes-1)
    class(pf_level_t), pointer :: lev

    lev => pf%levels(level_index)
    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))

    call start_timer(pf, TLEVEL+lev%index-1)
    do k = 1,nsweeps
       call call_hooks(pf, level_index, PF_PRE_SWEEP)    
       call lev%Q(m)%copy(lev%q0)
       do m = 1, lev%nnodes-1
          call lev%R(m)%copy(lev%Q(m+1))
          call lev%Q(m+1)%copy(lev%q0)
       end do
       
       do m = 1, lev%nnodes-1
          ! compute new Q from previous omegas
          call this%propagate(lev%Q(1), this%omega(m), lev%Q(m+1))
          call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m,1))
          call this%dexpinv(this%omega(m), lev%F(m,1), this%dexpinvs(m))
          t=t+dtsdc(m)
       end do
       call pf_residual(pf, lev, dt)   
       call lev%qend%copy(lev%Q(lev%nnodes))
       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do
    
    ! call imk_integrate(this, lev, this%omega, this%dexpinvs, dt, lev%S)
    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine imk_sweep

  subroutine imk_initialize(this, lev)
    class(pf_imk_t),   intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer :: m, nnodes

    nnodes = lev%nnodes

    call lev%ulevel%factory%create_array(this%omega, nnodes-1, &
         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call lev%ulevel%factory%create_array(this%dexpinvs, nnodes-1, &
         lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    do m = 1, nnodes-1
       call this%omega(m)%setval(0.0_pfdp)
       call this%dexpinvs(m)%setval(0.0_pfdp)
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
          call fintSDC(m)%axpy(dt*lev%qmat(m,j), fSDC(j,1))
       end do
    end do

  end subroutine imk_integrate

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

      call lev%ulevel%factory%destroy_array(this%dexpinvs, lev%nnodes-1, &
           lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

  end subroutine imk_destroy

end module pf_mod_imk
