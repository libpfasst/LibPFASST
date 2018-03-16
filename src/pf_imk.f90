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
!>  This module implements fully implicit Munthe-Kaas Runge Kutta methods using explicit SDC sweeping
!!
!!  The equation to be solved is 
!!
!! y'=A(y,t)y
!!
!! where A is a matrix and y is  a vector or matrix or if Lax_pair = true
!!
!! Y'=[A(Y,t),Y] where both A and Y are matrices
!!
!!  We solve this by finding the solution to
!!
!!  Q' = dexpinv_Q(A)
!!
!!  Using PFASST
module pf_mod_imk
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_imk_t
     class(pf_encap_t), allocatable :: Y(:)
     class(pf_encap_t), allocatable :: A(:)
     class(pf_encap_t), allocatable :: mexp(:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: dtsdc(:)
     real(pfdp), allocatable :: tsdc(:)
     real(pfdp) :: bernoullis(20)
     integer ::  qtype
     logical ::  Lax_pair
     logical ::  use_SDC
     integer ::  nterms
  contains
    procedure :: sweep        => imk_sweep
    procedure :: initialize   => imk_initialize
    procedure :: evaluate     => imk_evaluate
    procedure :: integrate    => imk_integrate
    procedure :: residual     => imk_residual
    procedure :: spreadq0     => imk_spreadq0
    procedure :: evaluate_all => imk_evaluate_all
    procedure :: imk_destroy
    procedure(pf_f_eval_p), deferred :: f_eval
    procedure(pf_dexpinv_p), deferred :: dexpinv
    procedure(pf_propagate_p), deferred :: propagate
 end type pf_imk_t

 interface
    !>  Subroutine f_eval computes A(y,t)
     subroutine pf_f_eval_p(this, y, t, level, f)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t),   intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y, f
       real(pfdp),        intent(in   ) :: t
       integer(c_int),    intent(in   ) :: level
     end subroutine pf_f_eval_p
    !>  Subroutine dexpinv computes Om'=F=dexpinv_Om(A)
     subroutine pf_dexpinv_p(this, a, omega, f)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: a
       class(pf_encap_t), intent(inout) :: omega
       class(pf_encap_t), intent(inout) :: f       !<  The resultign-level
     end subroutine pf_dexpinv_p
    !>  Subroutine propagate   computes y_m=expm(Om_m)y_0(expm(Om_m))-1 or (expm(Om_m))y_0 or
     subroutine pf_propagate_p(this, q0, q)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: q, q0
     end subroutine pf_propagate_p
  end interface

contains
  !> Perform nsweep  sweeps on level  and set qend appropriately.
  ! with the two-array encap, things are a little tricky
  ! copy default behavior : copy the solution only
  ! copy flagged behavior : copy the name of the encap
  ! setval default behavior : set the value of the name of the encap
  ! setval flagged behavior : set the value of the solution
  subroutine imk_sweep(this, pf, level_index, t0, dt,nsweeps)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf      !<  PFASST structure
    integer,             intent(in)    :: level_index  !<  which level this is
    real(pfdp),        intent(in   ) :: dt             !<  time step size
    real(pfdp),        intent(in   ) :: t0             !<  Time at beginning of time step
    integer,             intent(in)    :: nsweeps      !<  number of sweeps to do

    integer :: n,m,k                      !< Loop counters
    class(pf_level_t), pointer :: lev     !<  points to current level

    lev => pf%levels(level_index) !  Assign level pointer
    this%tsdc = t0+dt*lev%nodes   !  Set times at nodes

    call start_timer(pf, TLEVEL+lev%index-1)

    ! compute residual
    call pf_residual(pf, lev, dt)

    do k = 1,nsweeps
       call call_hooks(pf, level_index, PF_PRE_SWEEP)

       ! Store all function values  (dexpinv)
       call imk_save(lev)

       !  If first sweep, compute first function value should is just A(t_0,y_0)
       if (k .eq. 1) then
          call lev%Q(1)%copy(lev%q0, 1) !copy solution, flag
          call lev%Q(1)%setval(0.0d0) !setval omega, no flag
          call imk_evaluate(this, lev, t0, 1)
       end if

       ! Assign the old function value the difference in function values
       call lev%pF(1,1)%axpy(-1.0_pfdp,lev%F(1,1))

       do m = 1, lev%nnodes-1
          !>  Accumulate rhs in Residual
          if (this%use_sdc .eqv. .true.) then
              do n = 1, m
                call lev%R(m)%axpy(-dt*this%QtilE(m,n), lev%pF(n,1))
              end do
          end if

          !  Add the starting old omega term
          call lev%R(m)%axpy(1.0_pfdp, lev%Q(m+1))
          call lev%Q(m+1)%copy(lev%R(m)) !copy omega, no flag

          !  Evaluate the new rhs  (will also do Y)
          call imk_evaluate(this, lev, this%tsdc(m+1), m+1)
          call lev%pF(m+1,1)%axpy(-1.0_pfdp,lev%F(m+1,1))
       end do

       call pf_residual(pf, lev, dt)
       call lev%qend%copy(lev%Q(lev%nnodes), 1)  !<  copy solution, flag
       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do

    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine imk_sweep

  subroutine imk_initialize(this, lev)
    class(pf_imk_t),   intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer :: n,m, nnodes

    this%npieces = 1

    allocate(this%QtilE(nnodes-1,nnodes))
    allocate(this%dtsdc(nnodes-1))
    allocate(this%tsdc(nnodes))

    nnodes = lev%nnodes
    this%dtsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    this%bernoullis = 0.0_pfdp
    this%bernoullis(1 ) =       -1.0d0 / 2.0d0
    this%bernoullis(2 ) =        1.0d0 / 6.0d0
    this%bernoullis(4 ) =       -1.0d0 / 3.0d1
    this%bernoullis(6 ) =        1.0d0 / 4.2d1
    this%bernoullis(8 ) =       -1.0d0 / 3.0d1
    this%bernoullis(10) =        5.0d0 / 6.6d1
    this%bernoullis(12) =     -691.0d0 / 2.73d3
    this%bernoullis(14) =        7.0d0 / 6.0d0
    this%bernoullis(16) =    -3617.0d0 / 5.10d2
    this%bernoullis(18) =    43867.0d0 / 7.98d2
    this%bernoullis(20) =    -174611.0d0/330.0d0
    !>  Assign explicit approximate quadrature rule
    this%QtilE =  lev%qmatFE

    !>  Make space for temporary variables
    call lev%ulevel%factory%create_array(this%Y, nnodes, &
         lev%index,   lev%shape)

    call lev%ulevel%factory%create_array(this%A, nnodes, &
         lev%index,   lev%shape)

    call lev%ulevel%factory%create_array(this%mexp, nnodes, &
         lev%index,   lev%shape)

    do m = 1, nnodes
       call this%Y(m)%setval(0.0_pfdp)
       call this%A(m)%setval(0.0_pfdp)
       call this%mexp(m)%setval(0.0_pfdp)
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

    !  Propagate to get y=exp(Om)
    !prop needs e^{Q (omega)} and apply to Y
    if (m > 1) then
       call this%propagate(lev%q0, this%Y(m),lev%Q(m))
    else
       call this%A(1)%copy(lev%q0)
    end if

    !  Compute A(y,t)
    call this%f_eval(lev%Q(m), t, lev%index, this%A(m))
    if (this%debug) print*, 'A'
    if (this%debug) call this%A(m)%eprint()

    !  Compute the series expansion for dexpinv
    if (m > 1)  then
       call this%dexpinv(this%A(m), lev%Q(m), lev%F(m,1))
    else
       call lev%F(1,1)%copy(this%A(1))
    endif
    if (this%debug) print*, 'depxinv'
    if (this%debug) call lev%F(m,1)%eprint()
    if (this%debug) print*, 'end evaluate ------------'
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

    call lev%ulevel%sweeper%integrate(lev, lev%Q, lev%F, dt, lev%I)

    ! add tau (which is 'node to node')
    if (allocated(lev%tauQ)) then
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))
       end do
    end if

    ! subtract out Q  (not initial condition is zero
    do m = 1, lev%nnodes-1
       call lev%R(m)%copy(lev%I(m))
       call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
    end do

  end subroutine imk_residual

  subroutine imk_spreadq0(this, lev, t0)
    class(pf_imk_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t0

    integer m,p
    !  Stick initial condition into first node slot
    call lev%Q(1)%copy(lev%q0, 1)

    call lev%Q(1)%setval(0.0d0)

    !  Evaluate F at first spot
    call lev%ulevel%sweeper%evaluate(lev, t0, 1)

    ! Spread F and solution to all nodes
    do m = 2, lev%nnodes
       call lev%Q(m)%copy(lev%Q(1))
       do p = 1, lev%ulevel%sweeper%npieces
         call lev%F(m,p)%copy(lev%F(1,p))
       end do
    end do

  end subroutine imk_spreadq0

  !>  Save function values so that difference can be computed
  subroutine imk_save(lev)
    class(pf_level_t), intent(inout) :: lev  !<  Level to save on

    integer :: m, p

    do m = 1, lev%nnodes
       call lev%pF(m,1)%copy(lev%F(m,1))
    end do
  end subroutine imk_save

  subroutine imk_destroy(this, lev)
      class(pf_imk_t),   intent(inout) :: this
      class(pf_level_t), intent(inout) :: lev

      deallocate(this%QtilE)
      deallocate(this%dtsdc)
      deallocate(this%tsdc)

      call lev%ulevel%factory%destroy_array(this%Y, lev%nnodes, &
           lev%index,   lev%shape)

      call lev%ulevel%factory%destroy_array(this%A, lev%nnodes, &
           lev%index,   lev%shape)

      call lev%ulevel%factory%destroy_array(this%mexp, lev%nnodes, &
           lev%index,  lev%shape)

  end subroutine imk_destroy

end module pf_mod_imk
