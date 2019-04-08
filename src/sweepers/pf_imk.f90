!!  Implicit Munthe-Kass Runge-Kutta sweeper
!
! This file is part of LIBPFASST.
!
!>  This module implements fully implicit Munthe-Kaas Runge Kutta methods using explicit SDC sweeping
!!
!!  The equation to be solved is 
!! 
!! $$ y'=A(y,t)y  $$
!!
!! where \(A\) is a matrix and \(y)\ is  a vector or matrix or if Lax_pair = true
!!
!! $$Y'=[A(Y,t),Y]$$ where both \(A\) and \(Y\) are matrices
!!
!!  We solve this by finding the solution to
!!
!!  $$Q' = dexpinv_Q(A)$$
!!
!!  Using PFASST
module pf_mod_imk
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  !>  Implicit Munthe-Kaas Runge-Kutta sweeper type, extends abstract sweeper
  type, extends(pf_sweeper_t), abstract :: pf_imk_t
     class(pf_encap_t), allocatable :: A(:)
     real(pfdp), allocatable :: QtilE(:,:)   !!  Aproximate quadrature matric
     real(pfdp), allocatable :: QdiffE(:,:)  !!  qmat-QtilE
     real(pfdp), allocatable :: dtsdc(:)     !!  SDC substep size
     real(pfdp), allocatable :: tsdc(:)
     real(pfdp) :: bernoullis(20)  !!  Bernoulli numbers
     real(pfdp) :: t0   !!  Time at beginning of time step
     real(pfdp) :: dt   !!  Time step size
     integer ::  qtype, nterms
     logical ::  Lax_pair, use_SDC, debug, mkrk, rk
  contains
    procedure :: sweep        => imk_sweep
    procedure :: initialize   => imk_initialize
    procedure :: evaluate     => imk_evaluate
    procedure :: integrate    => imk_integrate
    procedure :: residual     => imk_residual
    procedure :: spreadq0     => imk_spreadq0
    procedure :: evaluate_all => imk_evaluate_all
    procedure :: destroy   => imk_destroy
    procedure(pf_f_eval_p), deferred :: f_eval
    procedure(pf_dexpinv_p), deferred :: dexpinv
    procedure(pf_propagate_p), deferred :: propagate
    procedure(pf_commutator_p), deferred :: commutator_p
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
       class(pf_encap_t), intent(inout) :: f       !!  The result
     end subroutine pf_dexpinv_p
    !>  Subroutine propagate   computes y_m=expm(Om_m)y_0(expm(Om_m))-1 or (expm(Om_m))y_0 or
     subroutine pf_propagate_p(this, q0, q)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: q, q0
     end subroutine pf_propagate_p
     subroutine pf_commutator_p(this, a, b, out, flags)
       import pf_imk_t, pf_encap_t, c_int, pfdp
       class(pf_imk_t), intent(inout) :: this
       class(pf_encap_t), intent(inout) :: a, b, out
       integer, intent(in), optional :: flags
     end subroutine pf_commutator_p
  end interface

contains
  !> Perform nsweep  sweeps on level  and set qend appropriately.
  ! with the two-array encap, things are a little tricky
  ! copy default behavior : copy the solution only
  ! copy flagged behavior : copy the name of the encap
  ! setval default behavior : set the value of the name of the encap
  ! setval flagged behavior : set the value of the solution
  subroutine imk_sweep(this, pf, level_index, t0, dt, nsweeps, flags)
    use pf_mod_timer
    use pf_mod_hooks

    !>  Inputs
    class(pf_imk_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf      !!  PFASST structure
    integer,             intent(in)    :: level_index  !!  which level to sweep on
    real(pfdp),        intent(in   ) :: t0             !!  Time at beginning of time step
    real(pfdp),        intent(in   ) :: dt             !!  time step size
    integer,             intent(in)    :: nsweeps      !!  number of sweeps to do
    integer, optional,   intent(in)    :: flags

    !>  Local variables
    type(pf_level_t), pointer :: lev    !!  points to current level

    this%t0 = t0
    this%dt = dt
    if (this%rk) then
       call rk_step(this, pf, t0, dt)
    else if (this%mkrk) then
       call mkrk_step(this, pf, t0, dt)
    else
       call imk_actually_sweep(this, pf, level_index, t0, dt, nsweeps)
    end if
  end subroutine imk_sweep

  subroutine rk_step(this, pf, t0, dt)
    use pf_mod_timer
    use pf_mod_hooks

    !>  Inputs
    class(pf_imk_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf      !!  PFASST structure
    real(pfdp),        intent(in   ) :: t0             !!  Time at beginning of time step
    real(pfdp),        intent(in   ) :: dt             !!  time step size

    !>  Local variables
    type(pf_level_t), pointer :: lev    !!  points to current level

    integer     :: level_index   !!  Level to work on
    integer     :: m        !!  Loop variables
    real(pfdp)  :: t        !!  Time at nodes

    t = t0 + dt
    level_index=1
    lev => pf%levels(level_index)
    do m = 1, 5
       call lev%Q(m)%setval(0.0_pfdp)
    end do

    call call_hooks(pf, 1, PF_PRE_SWEEP)

    call lev%Q(1)%copy(lev%q0, flags=1)
    call this%f_eval(lev%Q(1), t, lev%index, this%A(1))
    ! commutator_p flags=1 hack copies Q(1)%y -> Q(1)%array
    ! all subsequent RK stages are done on Q(m)%array
    call this%commutator_p(this%A(1), lev%Q(1), lev%F(1,1), flags=1)

    call lev%Q(2)%axpy(1.0_pfdp, lev%Q(1))
    call lev%Q(2)%axpy(0.5_pfdp*dt, lev%F(1,1))
    call this%f_eval(lev%Q(2), t, lev%index, this%A(2))
    call this%commutator_p(this%A(2), lev%Q(2), lev%F(2,1))

    call lev%Q(3)%axpy(1.0_pfdp, lev%Q(1))
    call lev%Q(3)%axpy(0.5_pfdp*dt, lev%F(2,1))
    call this%f_eval(lev%Q(3), t, lev%index, this%A(3))
    call this%commutator_p(this%A(3), lev%Q(3), lev%F(3,1))

    call lev%Q(4)%axpy(1.0_pfdp, lev%Q(1))
    call lev%Q(4)%axpy(dt, lev%F(3,1))
    call this%f_eval(lev%Q(4), t, lev%index, this%A(4))
    call this%commutator_p(this%A(4), lev%Q(4), lev%F(4,1))

    call lev%Q(5)%axpy(1.0_pfdp, lev%Q(1))
    call lev%Q(5)%axpy(dt/6.0_pfdp, lev%F(1,1))
    call lev%Q(5)%axpy(dt/3.0_pfdp, lev%F(2,1))
    call lev%Q(5)%axpy(dt/3.0_pfdp, lev%F(3,1))
    call lev%Q(5)%axpy(dt/6.0_pfdp, lev%F(4,1))

    ! commutator_p flags=2 hack copies Q(m)%array -> Q(m)%y
    ! only the Q argument in this case matters
    ! in order to get solution back to qend
    call this%commutator_p(this%A(1), lev%Q(1), lev%F(1,1), flags=2)
    call this%commutator_p(this%A(2), lev%Q(2), lev%F(1,1), flags=2)
    call this%commutator_p(this%A(3), lev%Q(3), lev%F(1,1), flags=2)
    call this%commutator_p(this%A(4), lev%Q(4), lev%F(1,1), flags=2)
    call this%commutator_p(this%A(5), lev%Q(5), lev%F(1,1), flags=2)

    call pf_residual(pf, level_index, dt)
    call lev%qend%copy(lev%Q(lev%nnodes), flags=1)

    call call_hooks(pf, 1, PF_POST_SWEEP)
  end subroutine rk_step

  subroutine mkrk_step(this, pf, t0, dt)
    use pf_mod_timer
    use pf_mod_hooks

    !>  Inputs
    class(pf_imk_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf      !!  PFASST structure
    real(pfdp),        intent(in   ) :: t0             !!  Time at beginning of time step
    real(pfdp),        intent(in   ) :: dt             !!  time step size

    !>  Local variables
    type(pf_level_t), pointer :: lev    !!  points to current level

    integer     :: level_index !!  Level we are working on
    integer     :: m        !!  Loop variables
    real(pfdp)  :: t        !!  Time at nodes

    level_index=1
    lev => pf%levels(level_index)

    t = t0 + dt
    do m = 1, 5
       call lev%Q(m)%setval(0.0_pfdp)
    end do

    call call_hooks(pf, 1, PF_PRE_SWEEP)

    call lev%Q(1)%copy(lev%q0, flags=1)
    call this%f_eval(lev%Q(1), t, lev%index, this%A(1))
    call this%dexpinv(this%A(1), lev%Q(1), lev%F(1,1))

    call lev%Q(2)%axpy(0.5_pfdp*dt, lev%F(1,1))
    call this%propagate(lev%q0, lev%Q(2))
    call this%f_eval(lev%Q(2), t, lev%index, this%A(2))
    call this%dexpinv(this%A(2), lev%Q(2), lev%F(2,1))

    call lev%Q(3)%axpy(0.5_pfdp*dt, lev%F(2,1))
    call this%propagate(lev%q0, lev%Q(3))
    call this%f_eval(lev%Q(3), t, lev%index, this%A(3))
    call this%dexpinv(this%A(3), lev%Q(3), lev%F(3,1))

    call lev%Q(4)%axpy(dt, lev%F(3,1))
    call this%propagate(lev%q0, lev%Q(4))
    call this%f_eval(lev%Q(4), t, lev%index, this%A(4))
    call this%dexpinv(this%A(4), lev%Q(4), lev%F(4,1))

    call lev%Q(5)%axpy(dt/6.0_pfdp, lev%F(1,1))
    call lev%Q(5)%axpy(dt/6.0_pfdp, lev%F(4,1))
    call lev%Q(5)%axpy(dt/3.0_pfdp, lev%F(2,1))
    call lev%Q(5)%axpy(dt/3.0_pfdp, lev%F(3,1))

    call this%propagate(lev%q0, lev%Q(5))

    call pf_residual(pf, level_index, dt)
    call lev%qend%copy(lev%Q(lev%nnodes), flags=1)

    call call_hooks(pf, 1, PF_POST_SWEEP)
  end subroutine mkrk_step

  subroutine imk_actually_sweep(this, pf, level_index, t0, dt, nsweeps)
    use pf_mod_timer
    use pf_mod_hooks

    !>  Inputs
    class(pf_imk_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf      !!  PFASST structure
    integer,             intent(in)    :: level_index  !!  which level to sweep on
    real(pfdp),        intent(in   ) :: t0             !!  Time at beginning of time step
    real(pfdp),        intent(in   ) :: dt             !!  time step size
    integer,             intent(in)    :: nsweeps      !!  number of sweeps to do

    !>  Local variables
    type(pf_level_t), pointer :: lev    !!  points to current level

    integer     :: m, n,k   !!  Loop variables
    real(pfdp)  :: t        !!  Time at nodes
    lev => pf%levels(level_index)   !!  Assign level pointer

    call start_timer(pf, TLEVEL+lev%index-1)

    do k = 1,nsweeps   !!  Loop over sweeps
       pf%state%sweep=k
       call call_hooks(pf, level_index, PF_PRE_SWEEP)
       ! compute integrals and add fas correction
       do m = 1, lev%nnodes-1
          call lev%I(m)%setval(0.0_pfdp)
          do n = 1, lev%nnodes
             call lev%I(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
          end do
          if (level_index < pf%nlevels) then
             call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))
          end if
       end do

       !  Recompute the first function value if this is first sweep
       if (k .eq. 1) then
          call lev%Q(1)%setval(0.0_pfdp) ! likely an unnecessary setting of Omega=0
          call this%evaluate(pf,level_index, t0, 1)
       end if

       t = t0
       ! do the sub-stepping in sweep
       do m = 1, lev%nnodes-1  !!  Loop over substeps
          t = t + dt*this%dtsdc(m)

          !>  Accumulate rhs
          call lev%Q(m+1)%setval(0.0_pfdp)
          do n = 1, m
             call lev%Q(m+1)%axpy(dt*this%QtilE(m,n), lev%F(n,1))
          end do
          !>  Add the tau term
          call lev%Q(m+1)%axpy(1.0_pfdp, lev%I(m))

          !>  Compute explicit function on new value
          call this%evaluate(pf,level_index, t, m+1)
       end do  !!  End substep loop

       call pf_residual(pf, level_index, dt)
       call lev%qend%copy(lev%Q(lev%nnodes), flags=1)

    end do  !  End loop on sweeps

    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine imk_actually_sweep

  subroutine imk_initialize(this, pf,level_index)
    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize

    integer :: m, nnodes
    type(pf_level_t), pointer  :: lev        !!  Current level
    lev => pf%levels(level_index)   !!  Assign level pointer

    this%npieces = 1
    nnodes = lev%nnodes

    allocate(this%QdiffE(nnodes-1,nnodes), this%QtilE(nnodes-1,nnodes))
    allocate(this%dtsdc(nnodes-1))
    allocate(this%tsdc(nnodes))


    this%dtsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    this%bernoullis = 0.0_pfdp
    this%bernoullis(1 ) =       -1.0_pfdp / 2.0_pfdp
    this%bernoullis(2 ) =        1.0_pfdp / 6.0_pfdp
    this%bernoullis(4 ) =       -1.0_pfdp / 3.0e1_pfdp
    this%bernoullis(6 ) =        1.0_pfdp / 4.2e1_pfdp
    this%bernoullis(8 ) =       -1.0_pfdp / 3.0e1_pfdp
    this%bernoullis(10) =        5.0_pfdp / 6.6e1_pfdp
    this%bernoullis(12) =     -691.0_pfdp / 2.73e3_pfdp
    this%bernoullis(14) =        7.0_pfdp / 6.0_pfdp
    this%bernoullis(16) =    -3617.0_pfdp / 5.10e2_pfdp
    this%bernoullis(18) =    43867.0_pfdp / 7.98e2_pfdp
    this%bernoullis(20) =    -174611.0_pfdp/330.0_pfdp
    !>  Assign explicit approximate quadrature rule
    this%QtilE =  lev%sdcmats%qmatFE
    this%QdiffE = lev%sdcmats%qmat-this%QtilE

    !>  Make space for temporary variables
    call lev%ulevel%factory%create_array(this%A, nnodes, &
         lev%index,   lev%shape)

    do m = 1, nnodes
       call this%A(m)%setval(0.0_pfdp)
    end do

  end subroutine imk_initialize

  subroutine imk_integrate(this, pf,level_index, qSDC, fSDC, dt, fintSDC, flags)
    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in   ) :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)
    integer, optional,   intent(in)    :: flags

    integer :: j, m
    type(pf_level_t), pointer  :: lev        !!  Current level
    lev => pf%levels(level_index)   !!  Assign level pointer

    do m = 1, lev%nnodes-1
       call fintSDC(m)%setval(0.0_pfdp)
       do j = 1, lev%nnodes
          call fintSDC(m)%axpy(dt*lev%sdcmats%qmat(m,j), fSDC(j,1))
       end do
    end do

  end subroutine imk_integrate

  subroutine imk_evaluate(this, pf,level_index, t, m, flags, step)
    use pf_mod_dtype
    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    real(pfdp),        intent(in   ) :: t
    integer,           intent(in   ) :: m
    integer, optional,   intent(in)    :: flags, step

    integer :: i
    real(pfdp) :: dt
    type(pf_level_t), pointer  :: lev        !!  Current level
    lev => pf%levels(level_index)   !!  Assign level pointer

    !  Propagate to get y=exp(Om)
    !prop needs e^{Q (omega)} and apply to Y
    if (this%debug) &
         print*, 'level', lev%index, 'in evaluate ', m, '------------------'

    if (this%rk) then
       ! 't' in f_evals are meaningless since i have a time-independent matrix, A
       dt = this%dt
       do i = 1, 5
         call lev%Q(i)%setval(0.0_pfdp)
       end do

       call lev%Q(1)%copy(lev%q0, flags=1)
       call this%f_eval(lev%Q(1), t, lev%index, this%A(1))
       ! commutator_p flags=1 hack copies Q(1)%y -> Q(1)%array
       ! all subsequent RK stages are done on Q(m)%array
       call this%commutator_p(this%A(1), lev%Q(1), lev%F(1,1), flags=1)

       call lev%Q(2)%axpy(1.0_pfdp, lev%Q(1))
       call lev%Q(2)%axpy(0.5_pfdp*dt, lev%F(1,1))
       call this%f_eval(lev%Q(2), t, lev%index, this%A(2))
       call this%commutator_p(this%A(2), lev%Q(2), lev%F(2,1))

       call lev%Q(3)%axpy(1.0_pfdp, lev%Q(1))
       call lev%Q(3)%axpy(0.5_pfdp*dt, lev%F(2,1))
       call this%f_eval(lev%Q(3), t, lev%index, this%A(3))
       call this%commutator_p(this%A(3), lev%Q(3), lev%F(3,1))

       call lev%Q(4)%axpy(1.0_pfdp, lev%Q(1))
       call lev%Q(4)%axpy(dt, lev%F(3,1))
       call this%f_eval(lev%Q(4), t, lev%index, this%A(4))
       call this%commutator_p(this%A(4), lev%Q(4), lev%F(4,1))

       call lev%Q(5)%axpy(1.0_pfdp, lev%Q(1))
       call lev%Q(5)%axpy(dt/6.0_pfdp, lev%F(1,1))
       call lev%Q(5)%axpy(dt/3.0_pfdp, lev%F(2,1))
       call lev%Q(5)%axpy(dt/3.0_pfdp, lev%F(3,1))
       call lev%Q(5)%axpy(dt/6.0_pfdp, lev%F(4,1))

       ! commutator_p flags=2 hack copies Q(m)%array -> Q(m)%y
       ! only the Q argument in this case matters
       ! in order to get solution back to qend
       call this%commutator_p(this%A(1), lev%Q(1), lev%F(1,1), flags=2)
       call this%commutator_p(this%A(2), lev%Q(2), lev%F(1,1), flags=2)
       call this%commutator_p(this%A(3), lev%Q(3), lev%F(1,1), flags=2)
       call this%commutator_p(this%A(4), lev%Q(4), lev%F(1,1), flags=2)
       call this%commutator_p(this%A(5), lev%Q(5), lev%F(1,1), flags=2)

    else
       if (this%debug) then
       endif
       if (m > 1) then
         if (this%debug) print*, 'propagating'
         call this%propagate(lev%q0, lev%Q(m))
       else
         if (this%debug) print*, 'copying'
         call lev%Q(1)%copy(lev%q0, flags=1)
       end if
       if (this%debug) print*, 'Q'
       if (this%debug) call lev%Q(m)%eprint()

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
    endif
    if (this%debug) print*, 'end evaluate ------------'
  end subroutine imk_evaluate

  subroutine imk_evaluate_all(this, pf,level_index, t, flags, step)
    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    real(pfdp),        intent(in   ) :: t(:)
    integer, optional,   intent(in)    :: flags, step

    integer :: m
    type(pf_level_t), pointer  :: lev        !!  Current level
    lev => pf%levels(level_index)   !!  Assign level pointer

    if (this%rk) then
       call lev%ulevel%sweeper%evaluate(pf,level_index, t(1), 1)
    else
       do m = 1, lev%nnodes
          call lev%ulevel%sweeper%evaluate(pf,level_index, t(m), m)
       end do
    end if
  end subroutine imk_evaluate_all

  subroutine imk_residual(this, pf,level_index, dt, flags)
    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),        intent(in   ) :: dt
    integer, optional,   intent(in)    :: flags
    integer :: m

    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !!  Assign level pointer    
    call lev%ulevel%sweeper%integrate(pf,level_index, lev%Q, lev%F, dt, lev%I)

    ! add tau (which is 'node to node')
    if (level_index < pf%nlevels) then
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

  subroutine imk_spreadq0(this, pf,level_index, t0, flags, step)
    class(pf_imk_t),  intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    real(pfdp),        intent(in   ) :: t0
    integer, optional,   intent(in)    :: flags, step

    integer m,p
    type(pf_level_t), pointer  :: lev        !!  Current level
    lev => pf%levels(level_index)   !!  Assign level pointer

    !  Stick initial condition into first node slot
    call lev%Q(1)%copy(lev%q0, flags=1)

    ! set initial omega to 0
    call lev%Q(1)%setval(0.0_pfdp)

    !  Evaluate F at first spot
    call lev%ulevel%sweeper%evaluate(pf,level_index, t0, 1)

    ! Spread F and solution to all nodes
    do m = 2, lev%nnodes
       call lev%Q(m)%copy(lev%Q(1), flags=1)
       call lev%Q(m)%copy(lev%Q(1))
       do p = 1, lev%ulevel%sweeper%npieces
         call lev%F(m,p)%copy(lev%F(1,p))
       end do
    end do

  end subroutine imk_spreadq0

  !>  Save function values so that difference can be computed
  subroutine imk_save(lev)
    type(pf_level_t), intent(inout) :: lev  !!  Level to save on

    integer :: m

    do m = 1, lev%nnodes
       call lev%pF(m,1)%copy(lev%F(m,1))
    end do
  end subroutine imk_save

  subroutine imk_destroy(this, pf,level_index)
    class(pf_imk_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize


    type(pf_level_t), pointer  :: lev        !!  Current level
    lev => pf%levels(level_index)   !!  Assign level pointer

      deallocate(this%QtilE, this%QdiffE)
      deallocate(this%dtsdc)
      deallocate(this%tsdc)

      call lev%ulevel%factory%destroy_array(this%A)
  end subroutine imk_destroy

end module pf_mod_imk
