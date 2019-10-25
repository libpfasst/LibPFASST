!!  Exponential integrator sweeper  module

! MODULE: pf_mod_exp
! !> @author
! Tommaso Buvoli
!
! Last Modified: Dec 28, 2018
!
!> Exponential integrator module
!!
!!  This module extends pf_sweeper_t and is used for creating an exponential sweeper 
!!  that solves equations of the form
!!         $$   y' = L y + N(t,y)  $$
!!  When extending this class, you must supply the functions phib, swpPhib, and resPhib
!!  that each compute matrix-vector products of the form
!!         $$ \sum_{i=0}^n t^i \varphi_i(tL)b_i $$
!!  in addition to the function f_eval for compluting the nonlinear term N(t,y).
!!  The complete description of these three functions is contained below.
module pf_mod_fexp

  use pf_mod_dtype
  use pf_mod_utils

  implicit none

  !> Exponential SDC sweeper type, extends abstract pf_sweeper_t
  type, extends(pf_sweeper_t), abstract :: pf_fexp_t

    real(pfdp),        allocatable :: nodes(:)   ! nodes
    real(pfdp),        allocatable :: eta(:)     ! normalized substeps (on interval [0, 1])
    class(pf_encap_t), allocatable :: F_old(:)   ! scratch space for storing nonlinear terms
    class(pf_encap_t), allocatable :: newF       ! scratch space for storing new function evaluations

  contains

    ! Specialized procedures for exponential integrator
    procedure(pf_f_eval_p),           deferred :: f_eval    ! computes nonlinear term in equation
    procedure(pf_expSweepSubstep),    deferred :: expSweepSubstep
    procedure(pf_expResidualSubstep), deferred :: expResidualSubstep

    ! Generic Functions
    procedure :: initialize   => exp_initialize
    procedure :: sweep        => exp_sweep
    procedure :: evaluate     => exp_evaluate
    procedure :: integrate    => exp_integrate
    procedure :: residual     => exp_residual
    procedure :: spreadq0     => exp_spreadq0
    procedure :: evaluate_all => exp_evaluate_all
    procedure :: destroy      => exp_destroy

    ! functions that can be accessed directly by types that inherit pf_fexp_t
    procedure :: exp_destroy
    procedure :: exp_initialize

  end type pf_fexp_t

  interface ! DESCRIPTION OF REQUIRED FUNCTIONS

    ! =================================================================================
    ! EXPSWEEPSUBSTEP: Computes the jth substep of the exponential sweeper
    !
    !     		exp(x L) y_{n,j}^{k+1} + h \varphi_1(x L) [N_{n,j}^{[k+1]} - N_{n,j}^{[k]}] + I_{n,j}
    !
    !	where x = dt * (t_{n,j+1} - t_{n,j}), I_{j} are the exponential integrals
    !
    !		I_{n,j} = x \int_{0}^{1} exp(x L(1-s)) Q_j(s) ds
    !
    !	and Q_j(s) is a Lagrange interpolating polynomial that satisfies
    !
    !		Q_j(\tau_{j,i}) = N(y(t_{n,i})) for i = 1 ... m
    !
    !	for \tau_{j,i} = (t_{n,i} - t_{n,j}) / (t{n,j+1} - t_{n,j})
    !
    !	NOTE: The operator L is not passed in as a parameter, and must be
    !       implemented appropriately within the class.
    !
    ! Arguments
    !
    !   j   	(input) Integer
    !       	substep index for determining t: t = t_{n,j+1} - t_{n,j}
    !
    !   y_jp1 (inout) pf_encap_t
    !		      stores the solution y_{n,j+1}^{k+1}
    !
    !   dt  	(input) real
    !       	stepsize
    !
    !   y_j  (input) pf_encap_t
    !		     stores the solution y_{n,j}^{k+1}
    !
    !   F    (input) pf_encap_t(:,:)
    !		     should be the lev%F matrix of dimension (nnodes-1 x 2). First component
    !        stores the nonlinear term, and the second component is used by this
    !        function to store the exponential integral terms I_{n,j}
    !
    !   Nk  	(input) pf_encap_t(:)
    !		stores N(y^{k+1}_j) for j = 1 ... m
    ! =================================================================================

    subroutine pf_expSweepSubstep(this, y_jp1, j, dt, y_j, F, Nk)
      import pf_fexp_t, pf_encap_t, pfdp
      class(pf_fexp_t),   intent(inout) :: this
      class(pf_encap_t),  intent(inout) :: y_jp1
      integer,            intent(in)    :: j
      real(pfdp),         intent(in)    :: dt
      class(pf_encap_t),  intent(in)    :: y_j
      class(pf_encap_t),  intent(in)    :: F(:,:)
      class(pf_encap_t),  intent(in)    :: Nk(:)
    end subroutine pf_expSweepSubstep

    ! =================================================================================
    ! EXPRESSUBSTEP: Computes the jth residual
    !
    ! Arguments
    !
    !   j   	(input) Integer
    !       	substep index for determining t: t = t_{n,j+1} - t_{n,j}
    !
    !   dt  	(input) real
    !       	stepsize
    !
    !   y_n   	(input) pf_encap_t
    !		stores the solution y_n^{k+1}
    !
    !   F    (input) pf_encap_t(:,:)
    !		     should be the lev%F matrix of dimension (nnodes-1 x 2). First component
    !        stores the nonlinear term, and the second component is used by this
    !        function to store the exponential integral terms I_{n,j}
    !
    ! =================================================================================

    subroutine pf_expResidualSubstep(this, y_np1, j, dt, y_n, F)
      import pf_fexp_t, pf_encap_t, pfdp
      class(pf_fexp_t),   intent(inout) :: this
      class(pf_encap_t),  intent(inout) :: y_np1
      integer,            intent(in)    :: j
      real(pfdp),         intent(in)    :: dt
      class(pf_encap_t),  intent(in)    :: y_n
      class(pf_encap_t),  intent(in)    :: F(:,:)
    end subroutine pf_expResidualSubstep

    ! =================================================================================
    ! f_eval: computes the equations nonlinear term N(t,y)
    !
    ! Arguments
    !
    !   y     (input) pf_encap_t
    !         solution y(t)
    !
    !   t     (input) DOUBLE
    !         time t
    !
    !   level_index (input) INTEGER
    !         current level index
    !
    !   f     (output) pf_encap_t
    !         N(t,y)
    ! =================================================================================

    subroutine pf_f_eval_p(this, y, t, level_index, n)
      import pf_fexp_t, pf_encap_t, pfdp
      class(pf_fexp_t),   intent(inout) :: this
      class(pf_encap_t), intent(in)    :: y
      real(pfdp),        intent(in)    :: t
      integer,           intent(in)    :: level_index
      class(pf_encap_t), intent(inout) :: n
    end subroutine pf_f_eval_p

  end interface

contains

  ! =================================================================================
  ! INITIALIZE: initializes the following internal parameters
  !      w        DOUBLE(:,:,:)     contains FD weights for computing local derivatives at t_{n,j}
  !      nodes    DOUBLE(:)         sdc nodes
  !      eta      DOUBLE(:)         normalized substeps (t_{n,j+1} - t_{n,j})/h
  !      npieces  INTEGER           number of RHS peices (always will be one)
  !      newF     pf_encap_t        stores new function evaluations
  !      b        pf_encap_t(:)     stores vectors b for computing phi products
  ! =================================================================================

  subroutine exp_initialize(this, pf,level_index)

    ! arguments
    class(pf_fexp_t),   intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index


    ! local variables
    integer :: i, nnodes
    real(pfdp), allocatable :: q(:)

    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)
    nnodes = lev%nnodes
    allocate(this%eta(nnodes - 1))
    allocate(this%nodes(nnodes))
    allocate(q(nnodes))
    ! set nodes and substeps
    this%nodes = lev%sdcmats%qnodes
    this%eta = this%nodes(2 : nnodes) - this%nodes(1 : nnodes - 1) ! substeps
    ! set number of rhs components
    this%npieces = 2
    ! initialize temporary storage objects
    call lev%ulevel%factory%create_single(this%newF, lev%index, lev%lev_shape)
    call lev%ulevel%factory%create_array(this%F_old, nnodes, lev%index, lev%lev_shape)

  end subroutine exp_initialize

  ! SWEEP: exponential sweep subroutine ===============================================
  subroutine exp_sweep(this, pf, level_index, t0, dt,nsweeps, flags)

    use pf_mod_timer
    use pf_mod_hooks

    ! arguments
    class(pf_fexp_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout), target :: pf   !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  which level to sweep on
    real(pfdp),        intent(in)    :: t0           !!  Time at beginning of time step
    real(pfdp),        intent(in)    :: dt           !!  time step size
    integer,           intent(in)    :: nsweeps      !!  number of sweeps to do
    integer, optional, intent(in)    :: flags

    ! local variables
    type(pf_level_t), pointer :: lev
    integer                    :: m, nnodes, j, k
    real(pfdp)                 :: t

    lev => pf%levels(level_index)
    nnodes = lev%nnodes
    ! error sweeps

    do k = 1, nsweeps
       call call_hooks(pf, level_index, PF_PRE_SWEEP)
       if (pf%save_timings > 1) call pf_start_timer(pf, T_SWEEP,level_index)
      pf%state%sweep=k

      ! NOTE: ensure that lev%F has been properly initialized here      
      do j = 1, nnodes
        call this%F_old(j)%copy(lev%F(j,1))  ! Save old f
      end do
      if (k .eq. 1) then ! Is this necessary? it seems that lev%F(j,1) = F(lev%(Q,q)) or F_old definition above line would be incorrect or is it zero initially?
         call lev%Q(1)%copy(lev%q0)
         if (pf%save_timings > 1) call pf_start_timer(pf, T_FEVAL,level_index)         
         call this%f_eval(lev%Q(1), t0, level_index, lev%F(1,1))      ! compute F_j^{[k+1]}
         if (pf%save_timings > 1) call pf_stop_timer(pf, T_FEVAL,level_index)         
      end if
      t = t0
      do j = 1, nnodes - 1
        t = t0 + dt * this%nodes(j+1)

        call this%expSweepSubstep(lev%Q(j+1), j, dt, lev%Q(j), lev%F, this%F_old)	! compute exp(h_j L) y_{j} + \varphi_1(h_j L) * (F^{k+1}_j - F^{k}_j) + I_j
        !  Now we have to add in the tauQ
        if (level_index < pf%state%finest_level) then
          call lev%Q(j+1)%axpy(1.0_pfdp, lev%tauQ(j))
          if (j > 1) then     ! The tau is not node to node, so subtract out
            call lev%Q(j+1)%axpy(-1.0_pfdp, lev%tauQ(j-1))
          end if
       end if
       if (pf%save_timings > 1) call pf_start_timer(pf, T_FEVAL,level_index)       
       call this%f_eval(lev%Q(j+1), t, level_index, lev%F(j+1,1))      			! compute F_j^{[k+1]}
       if (pf%save_timings > 1) call pf_stop_timer(pf, T_FEVAL,level_index)
      end do  !  Substepping over nodes


      call pf_residual(pf, level_index, dt)
      call lev%qend%copy(lev%Q(lev%nnodes))
      if (pf%save_timings > 1) call pf_stop_timer(pf, T_SWEEP,level_index)      
      call call_hooks(pf, level_index, PF_POST_SWEEP)

   end do  !  Sweeps
   
  end subroutine exp_sweep

  ! =================================================================================
  ! INTEGRATE: computes the integrals for the exponential Picard residual
  !
  !         r_j = A_j - B_j
  !
  !     where
  !
  !         A_j = \left[ \exp(t h L) y(t_n) - \int_{t_n}^{t_{n,j}} \exp(L(t - t_n)) P(t) dt \right]
  !         B_j = y(t_{n,j})
  !
  !     NOTE: This procedure computes the expression \hat{r}_j = A_j - y(t_n). The
  !     term y(t_n) is subtracted from result since the generic calling function
  !     compute the residual as
  !
  !        r_j = y(t_n) + exp_integrate() - y(t_{n,j})
  !
  !     thus incorrectly adding the term y_n
  ! =================================================================================

  subroutine exp_integrate(this, pf,level_index, qSDC, fSDC, dt, fintsdc, flags)
    class(pf_fexp_t),   intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    class(pf_encap_t), intent(in   ) :: qSDC(:)      !!  Solution values
    class(pf_encap_t), intent(in   ) :: fSDC(:, :)   !!  RHS Function values
    real(pfdp),        intent(in   ) :: dt           !!  Time step
    class(pf_encap_t), intent(inout) :: fintsdc(:)   !!  Integral from t_n to t_m
    integer, optional, intent(in   ) :: flags

    ! local variables
    integer :: i, nnodes
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    nnodes = lev%nnodes
    do i = 1, nnodes - 1 ! loop over integrals : compute \int_{t_{n,i}}^{t_{n, i + 1}}
      call this%expResidualSubstep(fintsdc(i), i, dt, qSDC(i), fSDC)
      call fintsdc(i)%axpy(-1.0_pfdp,qSDC(i))
      if (i > 1) then
        call fintsdc(i)%axpy(1.0_pfdp,fintsdc(i-1))
      end if
    end do

  end subroutine exp_integrate

  ! RESIDUAL: compute  residual (generic) ====================================
  subroutine exp_residual(this, pf, level_index, dt, flags)

    class(pf_fexp_t),  intent(inout)  :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),             intent(in)    :: dt
    integer, intent(in), optional    :: flags

    integer :: m
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    !>  Compute the integral of F from t_n to t_m at each node
    call lev%ulevel%sweeper%integrate(pf,level_index, lev%Q, lev%F, dt, lev%I, flags)

    !> add tau if it exists
    if (level_index < pf%state%finest_level) then
      do m = 1, lev%nnodes-1
        call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), flags)
      end do
    end if

    !> subtract out the solution value
    do m = 1, lev%nnodes-1
      call lev%R(m)%copy(lev%I(m))
      call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
      call lev%R(m)%axpy(1.0_pfdp, lev%Q(1))
    end do


  end subroutine exp_residual

  ! SPREADQ: spread solution (generic) ======================================
  subroutine exp_spreadq0(this, pf,level_index, t0, flags, step)
    class(pf_fexp_t),  intent(inout)  :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),        intent(in   ) :: t0
    integer, optional,   intent(in)  :: flags, step

    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    call pf_generic_spreadq0(this, pf,level_index, t0)

  end subroutine exp_spreadq0

  ! EVALUATE: evaluate the nonlinear term at node m ========================
  subroutine exp_evaluate(this, pf,level_index, t, m, flags, step)
    ! arguments
    class(pf_fexp_t),   intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),        intent(in   ) :: t    !!  Time at which to evaluate
    integer,           intent(in   ) :: m    !!  Node at which to evaluate
    integer, intent(in), optional    :: flags, step

    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer
    if (pf%save_timings > 1) call pf_start_timer(pf, T_FEVAL,level_index)
    call this%f_eval(lev%Q(m), t, level_index, lev%F(m,1))
    if (pf%save_timings > 1) call pf_stop_timer(pf, T_FEVAL,level_index)
  end subroutine exp_evaluate

  ! EVALUATE_ALL: evaluate the nonlinear term at all nodes =================
  subroutine exp_evaluate_all(this, pf,level_index, t, flags, step)
    ! arguments
    class(pf_fexp_t),  intent(inout)  :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),        intent(in   ) :: t(:)     !!  Array of times at each node
    integer, intent(in), optional    :: flags, step

    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    call pf_generic_evaluate_all(this,pf, level_index, t)

  end subroutine exp_evaluate_all

  ! DEALLOCATE: deallocate sweeper variables
  subroutine exp_destroy(this, pf,level_index)
    ! arguments
    class(pf_fexp_t),   intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer


    deallocate(this%eta)
    deallocate(this%newF)
    call lev%ulevel%factory%destroy_array(this%F_old)
  end subroutine exp_destroy


end module pf_mod_fexp
