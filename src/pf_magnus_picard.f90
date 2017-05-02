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
   contains
     procedure :: sweep      => magpicard_sweep
     procedure :: initialize => magpicard_initialize
     procedure :: evaluate   => magpicard_evaluate
     procedure :: integrate  => magpicard_integrate
     procedure :: residual   => magpicard_residual
     procedure :: evaluate_all => magpicard_evaluate_all
     procedure :: magpicard_destroy
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_compute_omega_p), deferred :: compute_omega
     procedure(pf_compute_time_ev_ops_p), deferred :: compute_time_ev_ops
     procedure(pf_propagate_solution_p), deferred :: propagate_solution
     procedure(pf_destroy_magpicard_p), deferred :: destroy
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
     subroutine pf_compute_omega_p(this, omega, f1, f2, coef)
       import pf_magpicard_t, pf_encap_t, c_int, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega
       class(pf_encap_t), intent(inout) :: f1, f2
       real(pfdp), intent(in) :: coef
     end subroutine pf_compute_omega_p
     subroutine pf_compute_time_ev_ops_p(this, omega, u, u_dagger)
       import pf_magpicard_t, pf_encap_t, c_int, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: omega
       class(pf_encap_t), intent(inout) :: u, u_dagger
     end subroutine pf_compute_time_ev_ops_p
     subroutine pf_propagate_solution_p(this, sol_t0, sol_tn, u, u_dagger)
       import pf_magpicard_t, pf_encap_t, c_int, pfdp
       class(pf_magpicard_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: sol_t0
       class(pf_encap_t), intent(inout) :: u
       class(pf_encap_t), intent(inout) :: u_dagger
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
    real(pfdp), allocatable :: coefs(:)
    complex(pfdp) :: z0=(0.0,0.0), z1=(1.0,0.0), zm1=(-1.0,0.0)
    class(pf_encap_t), allocatable :: time_ev_op(:), time_ev_op_dagger(:)

    call start_timer(pf, TLEVEL+lev%level-1)

    call lev%Q(1)%copy(lev%q0)

    call this%f_eval(lev%Q(1), t0, lev%level, lev%F(1,1))

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes
      t=t+dtsdc(m)
      call this%f_eval(lev%Q(m), t, lev%level, lev%F(m,1))
   end do

    ! computes the integral of F over the different sdc nodes, stores them on lev%S
    call magpicard_integrate(this, lev, lev%Q, lev%F, dt, lev%S)

    call lev%ulevel%factory%create_array(time_ev_op, &
         lev%nnodes-1, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%create_array(time_ev_op_dagger, &
         lev%nnodes-1, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    allocate(coefs(lev%nnodes))
    coefs(1) = dt*dt/48.0_pfdp
    coefs(2) = dt*dt/12.0_pfdp

    do m = 1, lev%nnodes-1
       do p = 1, this%npieces
          call this%compute_omega(lev%S(m), lev%F(m+1,p), lev%F(1,p), coefs(m))
          call this%compute_time_ev_ops(lev%S(m), time_ev_op(m), time_ev_op_dagger(m))
          call this%propagate_solution(lev%Q(1), lev%Q(m+1), time_ev_op(m), time_ev_op_dagger(m))
       enddo
    end do

    ! Put the last node value into qend
    call lev%qend%copy(lev%Q(lev%nnodes))

    call lev%ulevel%factory%destroy_array(time_ev_op, &
         lev%nnodes-1, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%destroy_array(time_ev_op_dagger, &
         lev%nnodes-1, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call end_timer(pf, TLEVEL+lev%level-1)
  end subroutine magpicard_sweep

  subroutine magpicard_initialize(this, lev)
    class(pf_magpicard_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer :: m,n,nnodes

    this%npieces = 1
    nnodes = lev%nnodes
    allocate(this%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(this%QtilI(nnodes-1,nnodes))  !  S-BE

    this%QtilI = 0.0_pfdp

    if(this%use_LUq_) then
        call myLUq(lev%qmat, this%QtilI, Nnodes, 1)
    else
        dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
        do m = 1, nnodes-1
            do n = 1,m
                this%QtilI(m,n+1) =  dsdc(n)
            end do
        end do
    end if

    this%QdiffI = lev%qmat-this%QtilI

!    print *,'QtilI',imp%QtilI
!    print *,'QdiffI',imp%QdiffI
!    print *,'Qmat',Lev%qmat

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
          ! this loop over pieces might actually end up being the order of expansion
          do p = 1, this%npieces
             call fintSDC(m)%axpy(dt*lev%qmat(m,j), fSDC(j,p))
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

    call this%f_eval(lev%Q(m), t, lev%level, lev%F(m,1))
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
    call pf_generic_residual(this, lev, dt)
  end subroutine magpicard_residual

  ! Destroy the matrices
  subroutine magpicard_destroy(this)
      class(pf_magpicard_t),  intent(inout) :: this

      deallocate(this%QdiffI)
      deallocate(this%QtilI)
  end subroutine magpicard_destroy
end module pf_mod_magnus_picard
