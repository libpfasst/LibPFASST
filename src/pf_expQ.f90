! =====================================================================================
! MODULE: pf_mod_exp
! !> @author
!> Tommaso Buvoli
!
! Last Modified: Dec 28, 2018
!
! Description
!
!!  this module extends pf_sweeper_t and is used for creating an exponential sweeper 
!!  that solves equations of the form
!!         $$   y' = L y + N(t,y)  $$
!!  When extending this class, you must supply the functions phib, swpPhib, and resPhib
!!  that each compute matrix-vector products of the form
!!         $$ \sum_{i=0}^n t^i \varphi_i(tL)b_i $$
!!  in addition to the function f_eval for compluting the nonlinear term N(t,y).
!!  The complete description of these three functions is contained below.
! =====================================================================================

module pf_mod_exp

use pf_mod_dtype
use pf_mod_utils

implicit none

! Exponential SDC sweeper type, extends abstract pf_sweeper_t
type, extends(pf_sweeper_t), abstract :: pf_exp_t

    real(pfdp),        allocatable :: w(:,:,:)   ! weights
    real(pfdp),        allocatable :: nodes(:)   ! nodes
    real(pfdp),        allocatable :: eta(:)     ! normalized substeps (on interval [0, 1])
    class(pf_encap_t), allocatable :: b(:)       ! scratch space for computing nonlinear derivatives
    class(pf_encap_t), allocatable :: f_old(:)   ! scratch space for storing nonlinear terms
    class(pf_encap_t), allocatable :: newF       ! scratch space for storing new function evaluations
    LOGICAL :: use_phib=.true.                   ! if TRUE calls phib otherwise calls swpPhib and resPhib (set in derived sweeper)

    contains

        ! specialized procedures for exponential integrator
        procedure(pf_f_eval_p), deferred :: f_eval    ! computes nonlinear term in equation
        procedure(pf_phib),     deferred :: phib      ! computes x(t) = \sum_{i=0}^n t^i \varphi_i(tL)b_i
        procedure(pf_swpPhib),  deferred :: swpPhib   ! computes x(t) = \sum_{i=0}^n t^i \varphi_i(tL)b_i where t = t_{n,j+1} - t_{n,j} for j = 1, ... q - 1
        procedure(pf_resPhib),  deferred :: resPhib   ! computes x(t) = \sum_{i=0}^n t^i \varphi_i(tL)b_i where t = t_{n,j+1} - t_{n} for j = 1, ... q - 1
        procedure, private :: weights
        procedure, private :: LocalDerivsAtNode
        ! generic functions
        procedure :: initialize   => exp_initialize
        procedure :: sweep        => exp_sweep
        procedure :: evaluate     => exp_evaluate
        procedure :: integrate    => exp_integrate
        procedure :: residual     => exp_residual
        procedure :: spreadq0     => exp_spreadq0
        procedure :: evaluate_all => exp_evaluate_all
        procedure :: destroy      => exp_destroy
        ! functions that can be accessed directly by types that inherit pf_exp_t
        procedure :: exp_destroy
        procedure :: exp_initialize

    end type pf_exp_t

    interface ! DESCRIPTION OF REQUIRED FUNCTIONS

        ! =================================================================================             REMARK: ADDING ax operation would simplify dealing with h for b(2:end)
        ! PHIB: Computes the product of vectors and phi functions
        !
        !       y(t) = exp(t h L) b_1 + h \sum_{i=1}^n t^i \varphi_i(t h L)b_{i+1}
        !
        !       for a user specified t, and h.
        !
        !       NOTE: The operator L is not passed in as a parameter, and must be
        !       implemented appropriately within the function.
        !
        ! Arguments
        !
        !   t   (input) DOUBLE
        !       evaluation time for expression
        !
        !   h   (input) DOUBLE
        !       scaling factor for the linear opeartor AND the vectors b_i for i = 2 ... n
        !
        !   b   (input) pf_encap_t(:)
        !       array that stores the vectors b_i
        !
        !   y   (output) pf_encap_t
        !       once subroutine terminate this stores the result y(t)
        ! =================================================================================

        subroutine pf_phib(this, t, h, b, y)
            import pf_exp_t, pf_encap_t, pfdp
            class(pf_exp_t),    intent(inout) :: this
            real(pfdp),         intent(in)    :: t
            real(pfdp),         intent(in)    :: h
            class(pf_encap_t),  intent(in)    :: b(:)
            class(pf_encap_t),  intent(inout) :: y
        end subroutine pf_phib

        ! =================================================================================
        ! SWPPHIB: Computes the product of vectors and phi functions
        !
        !           y(t) = \varphi_0(h L) b_1 + h \sum_{i=1}^n t^i \varphi_i(t h L)b_{i+1}
        !
        !       where the time t is
        !
        !           t_{n,j+1} - t_{n,j}         j = 1, ... q - 1
        !
        !       and h is a user specified. This procedure is used when computing 
        !       exponential correction sweeps.
        !
        !       NOTE: The operator L is not passed in as a parameter, and must be
        !       implemented appropriately within the function.
        !
        ! Arguments
        !
        !   j   (input) Integer
        !       substep index for determining t: t = t_{n,j+1} - t_{n,j}
        !
        !   h   (input) DOUBLE
        !       scaling factor for the linear opeartor AND the vectors b_i for i = 2 ... n
        !
        !   b   (input) pf_encap_t(:)
        !       array that stores the vectors b_i
        !
        !   y   (output) pf_encap_t
        !       once subroutine terminate this stores the result y(t)
        ! =================================================================================

        subroutine pf_swpPhib(this, j, h, b, y)
            import pf_exp_t, pf_encap_t, pfdp
            class(pf_exp_t),    intent(inout) :: this
            integer,            intent(in)    :: j
            real(pfdp),         intent(in)    :: h
            class(pf_encap_t),  intent(in)    :: b(:)
            class(pf_encap_t),  intent(inout) :: y
        end subroutine pf_swpPhib

        ! =================================================================================
        ! RESPHIB: Computes the product of vectors and phi functions
        !
        !           y(t) = \varphi_0(h L) b_1 + h \sum_{i=1}^n t^i \varphi_i(t h L)b_{i+1}
        !
        !       where the time t is
        !
        !           t_{n,j} - t_{n}
        !
        !       and h is a user specified. This procedure is used when computing 
        !       exponential correction sweeps.
        !
        !       NOTE: The operator L is not passed in as a parameter, and must be
        !       implemented appropriately within the function.
        !
        ! Arguments
        !
        !   j   (input) Integer
        !       substep index for determining t: t = t_{n,j+1} - t_{n}
        !
        !   h   (input) DOUBLE
        !       scaling factor for the linear opeartor AND the vectors b_i for i = 2 ... n
        !
        !   b   (input) pf_encap_t(:)
        !       array that stores the vectors b_i
        !
        !   y   (output) pf_encap_t
        !       once subroutine terminate this stores the result y(t)
        ! =================================================================================

        subroutine pf_resPhib(this, j, h, b, y)
            import pf_exp_t, pf_encap_t, pfdp
            class(pf_exp_t),    intent(inout) :: this
            integer,            intent(in)    :: j
            real(pfdp),         intent(in)    :: h
            class(pf_encap_t),  intent(in)    :: b(:)
            class(pf_encap_t),  intent(inout) :: y
        end subroutine pf_resPhib

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
        !   level (input) INTEGER
        !         current level index
        !
        !   f     (output) pf_encap_t
        !         N(t,y)
        ! =================================================================================

        subroutine pf_f_eval_p(this, y, t, level, n)
            import pf_exp_t, pf_encap_t, pfdp
            class(pf_exp_t),   intent(inout) :: this
            class(pf_encap_t), intent(in)    :: y
            real(pfdp),        intent(in)    :: t
            integer,           intent(in)    :: level
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

    subroutine exp_initialize(this, lev)

        ! arguments
        class(pf_exp_t),   intent(inout) :: this
        class(pf_level_t), intent(inout) :: lev

        ! local variables
        integer :: i, nnodes
        real(pfdp), allocatable :: q(:)

        nnodes = lev%nnodes
        allocate(this%eta(nnodes - 1))
        allocate(this%nodes(nnodes))
        allocate(q(nnodes))
        ! set nodes and substeps
        this%nodes = lev%sdcmats%qnodes
        this%eta = this%nodes(2 : nnodes) - this%nodes(1 : nnodes - 1) ! substeps
        ! compute weights
        allocate(this%w(nnodes - 1, nnodes, nnodes))
        do i = 1, nnodes - 1
            q = this%nodes - this%nodes(i);
            call weights(this, real(0.0, pfdp), q, nnodes - 1, this%W(i, :, :));
        end do
        ! set number of rhs components
        this%npieces = 1
        ! initialize temporary storage objects
        call lev%ulevel%factory%create_single(this%newF, lev%index, lev%shape)
        call lev%ulevel%factory%create_array(this%b, nnodes + 1, lev%index, lev%shape)
        call lev%ulevel%factory%create_array(this%f_old, nnodes, lev%index, lev%shape)

    end subroutine exp_initialize

    ! SWEEP: exponential sweep subroutine ===============================================
    subroutine exp_sweep(this, pf, level_index, t0, dt,nsweeps, flags)

        use pf_mod_timer
        use pf_mod_hooks

        ! arguments
        class(pf_exp_t),   intent(inout) :: this
        type(pf_pfasst_t), intent(inout), target :: pf   !!  PFASST structure
        integer,           intent(in)    :: level_index  !!  which level to sweep on
        real(pfdp),        intent(in)    :: t0           !!  Time at beginning of time step
        real(pfdp),        intent(in)    :: dt           !!  time step size
        integer,           intent(in)    :: nsweeps      !!  number of sweeps to do
        integer, optional, intent(in)    :: flags

        ! local variables
        class(pf_level_t), pointer :: lev
        integer                    :: m, nnodes, j, k,i
        real(pfdp)                 :: t

        lev => pf%levels(level_index)

        nnodes = lev%nnodes
        call lev%Q(1)%copy(lev%q0)  
        call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1))      ! compute F_j^{[k+1]}

        ! error sweeps
        do k = 1, nsweeps
           call call_hooks(pf, level_index, PF_PRE_SWEEP)      ! NOTE: ensure that lev%F has been properly initialized here
           t = t0 
           do j = 1, nnodes
              call this%f_old(j)%copy(lev%F(j,1))  ! Save old f
           end do
           do j = 1, nnodes - 1
              t = t0 + dt * this%eta(j)
              ! form b vectors
              call LocalDerivsAtNode(this, 1, nnodes, this%f_old(:), this%b(2:nnodes+1))  ! phi expansion for exponential picard integral
!              call LocalDerivsAtNode(this, 1, nnodes, lev%F(:,1), this%b(2:nnodes+1))  ! phi expansion for exponential picard integral
              call this%b(1)%copy(lev%Q(1))  ! add term \phi_0(tL) y_n

              do i=1,j
!                 call this%b(2)%axpy(-this%eta(i), this%f_old(i))         ! add -\phi_1(tL) F_j^{[k]}
!                 call this%b(2)%axpy(this%eta(i), lev%F(i,1))          ! add \phi_1(tL) F_j^{[k+1]}
                 call this%b(2)%axpy(-this%nodes(i+1), this%f_old(i))         ! add -\phi_1(tL) F_j^{[k]}
                 call this%b(2)%axpy(this%nodes(i+1), lev%F(i,1))          ! add \phi_1(tL) F_j^{[k+1]}
              end do
              
        call start_timer(pf, TLEVEL+lev%index-1)                      
              ! compute phi products
              if (this%use_phib) then
                 call this%phib(this%nodes(j+1), dt, this%b, lev%Q(j+1))
              else
!                 call this%swpPhib(j+1, dt, this%b, lev%Q(j+1))
                 call this%resPhib(j, dt, this%b, lev%Q(j+1))                 
              end if
    call end_timer(pf, TLEVEL+lev%index-1)        
!              !  Now we have to add in the tauQ
              if (allocated(lev%tauQ)) then
                 call lev%Q(j+1)%axpy(1.0_pfdp, lev%tauQ(j))
              end if

              call this%f_eval(lev%Q(j+1), t, lev%index, lev%F(j+1,1))                      ! eval last nonlinear term
!
           end do
           call lev%qend%copy(lev%Q(lev%nnodes))
           call pf_residual(pf, lev, dt)

           call call_hooks(pf, level_index, PF_POST_SWEEP)
           
        end do

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

    subroutine exp_integrate(this, lev, qSDC, fSDC, dt, fintsdc, flags)
        ! parameters
        class(pf_exp_t),   intent(inout) :: this
        class(pf_level_t), intent(in   ) :: lev          !!  Current level
        class(pf_encap_t), intent(in   ) :: qSDC(:)      !!  Solution values
        class(pf_encap_t), intent(in   ) :: fSDC(:, :)   !!  RHS Function values
        real(pfdp),        intent(in   ) :: dt           !!  Time step
        class(pf_encap_t), intent(inout) :: fintsdc(:)   !!  Integral from t_n to t_m
        integer, optional, intent(in   ) :: flags
        ! local variables
        integer :: i, nnodes

        nnodes = lev%nnodes
        call LocalDerivsAtNode(this, 1, nnodes, fSDC(:,1), this%b(2:nnodes+1)) ! compute derivatives
        call this%b(1)%setval(real(0.0, pfdp))
        call this%b(1)%axpy(real(1.0, pfdp), qSDC(1))
        do i = 1, nnodes - 1 ! loop over integrals : compute \int_{t_{n,i}}^{t_{n, i + 1}}
            if (this%use_phib) then
                call this%phib(this%nodes(i+1), dt, this%b, fintsdc(i))
            else
               call this%resPhib(i, dt, this%b, fintsdc(i))
             end if
        end do


    end subroutine exp_integrate

    ! RESIDUAL: compute  residual (generic) ====================================
    subroutine exp_residual(this, lev, dt, flags)

        class(pf_exp_t),  intent(inout)  :: this
        class(pf_level_t), intent(inout) :: lev  !!  Current level
        real(pfdp),        intent(in   ) :: dt   !!  Time step
        integer, intent(in), optional    :: flags

        integer :: m
        !>  Compute the integral of F from t_n to t_m at each node
        call lev%ulevel%sweeper%integrate(lev, lev%Q, lev%F, dt, lev%I, flags)
        
        !> add tau if it exists
        if (allocated(lev%tauQ)) then
           do m = 1, lev%nnodes-1
              call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), flags)
           end do
        end if
        
        !> subtract out the solution value
        do m = 1, lev%nnodes-1      
           call lev%R(m)%copy(lev%I(m))
           call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
        end do


    end subroutine exp_residual

    ! SPREADQ: spread solution (generic) ======================================
    subroutine exp_spreadq0(this, lev, t0, flags, step)

        class(pf_exp_t),  intent(inout)  :: this
        class(pf_level_t), intent(inout) :: lev
        real(pfdp),        intent(in   ) :: t0
        integer, optional,   intent(in)  :: flags, step

        call pf_generic_spreadq0(this, lev, t0)

    end subroutine exp_spreadq0

    ! EVALUATE: evaluate the nonlinear term at node m ========================
    subroutine exp_evaluate(this, lev, t, m, flags, step)
        ! arguments
        class(pf_exp_t),   intent(inout) :: this
        class(pf_level_t), intent(inout) :: lev  !!  Current level
        real(pfdp),        intent(in   ) :: t    !!  Time at which to evaluate
        integer,           intent(in   ) :: m    !!  Node at which to evaluate
        integer, intent(in), optional    :: flags, step

        call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1))

    end subroutine exp_evaluate

    ! EVALUATE_ALL: evaluate the nonlinear term at all nodes =================
    subroutine exp_evaluate_all(this, lev, t, flags, step)
        ! arguments
        class(pf_exp_t),  intent(inout)  :: this
        class(pf_level_t), intent(inout) :: lev   !!  Current level
        real(pfdp),        intent(in   ) :: t(:)     !!  Array of times at each node
        integer, intent(in), optional    :: flags, step

        call pf_generic_evaluate_all(this, lev, t)

    end subroutine exp_evaluate_all

    ! DEALLOCATE: deallocate sweeper variables
    subroutine exp_destroy(this, lev)
        ! arguments
        class(pf_exp_t),   intent(inout) :: this
        class(pf_level_t), intent(inout) :: lev   !!  Current level

        deallocate(this%w)
        deallocate(this%eta)
        deallocate(this%newF)
        call lev%ulevel%factory%destroy_array(this%b, lev%index, lev%nnodes,  lev%shape)
        call lev%ulevel%factory%destroy_array(this%f_old, lev%index, lev%nnodes,  lev%shape)        
    end subroutine exp_destroy

    ! =======================================================================
    ! LocalDerivsAtNode: approximate the local derivative vector 
    !           at the substep t_{n,i} for the nonlinear function N(y(t)) 
    !           using the terms N(y{n,i}). Local coordinates coorespond to 
    !           s = h t.
    ! Arguments
    !
    !   i       (input) INTEGER
    !           index of the substep for which we want to approximate derivatives
    !
    !   nnodes  (input) INTEGER
    !           number of nodes
    !
    !   N_eval  (input) pf_encap_t(:)
    !           Nonlinear function evaluations; N_eval(i) contains N(y_{n,i})
    !
    !   N_deriv (output) pf_encap_t(:)
    !           N_deriv(i) approximates \frac{d^{i-1}}{dt^{i-1}} N(y(t)) 
    ! =======================================================================

    subroutine LocalDerivsAtNode(this, i, nnodes, N_eval, N_deriv)
        ! arguments
        class(pf_exp_t),   intent(inout) :: this
        integer,           intent(in)    :: i
        integer,           intent(in)    :: nnodes
        class(pf_encap_t), intent(in)    :: N_eval(:)
        class(pf_encap_t), intent(inout) :: N_deriv(:)
        ! local variables
        integer :: j, k

        ! form nonlinear derivative vectors b
        do j = 1, nnodes                                                ! loop over derivatives j = 1 ... n
            call N_deriv(j)%setval(real(0.0, pfdp))
            do k = 1, nnodes                                            ! look over nodes k = 1 ... n
                call N_deriv(j)%axpy(this%w(i, k, j), N_eval(k))
            end do
        end do
    end

    ! =======================================================================
    ! WEIGHTS   Compute coefficients for finite difference approximation for
    !           the derivatives 1 to m at point z assuming data is known at
    !           points in array x. Based on the program "weights" in
    !           B. Fornberg, "Calculation of weights in finite difference
    !           formulas", SIAM Review 40 (1998), pp. 685-691.
    ! Arguments
    !
    !   z   (input) DOUBLE
    !       location where approximations are to be accurate
    !
    !   x   (input) DOUBLE Array
    !       array containing interpolation points
    !
    !   m   (input) INTEGER
    !       highest derivative for which weights are sought
    !
    !   W   (output) DOUBLE array, dimension(size(x),m+1)
    !       matrix that gives weights at grid locations x for
    !       derivative of order j<=m are found in c(:,j)
    ! =======================================================================

    subroutine weights(this, z, x, m, W)

        ! Arguments
        class(pf_exp_t),  intent(inout)  :: this
        real(pfdp), intent(in)    :: z
        real(pfdp), intent(in)    :: x(:)
        integer,    intent(in)    :: m
        real(pfdp), intent(out)   :: W(size(x),m+1)

        ! Variable Declarations
        real(pfdp) :: c1, c2, c3, c4, c5
        integer  :: i,j,k,n,mn

        c1 = 1.0_pfdp
        c4 = x(1) - z
        W  = 0.0_pfdp
        W(1,1) = 1.0_pfdp

        n = size(x)
        do i=2,n
            mn = min(i,m+1)
            c2 = 1.0_pfdp
            c5 = c4
            c4 = x(i) - z
            do j=1,i-1
                c3 = x(i) - x(j)
                c2 = c2*c3;
                if(j == i-1) then
                    do k=mn,2,-1
                        W(i,k) = c1*(real(k-1,pfdp)*W(i-1,k-1) - c5*W(i-1,k))/c2;
                    enddo
                    W(i,1) = -c1*c5*W(i-1,1)/c2;
                endif
                do k=mn,2,-1
                    W(j,k) = (c4*W(j,k) - real(k-1,pfdp)*W(j,k-1))/c3;
                enddo
                W(j,1) = c4*W(j,1)/c3;
            enddo
            c1 = c2;
        enddo

end subroutine weights

end module pf_mod_exp
