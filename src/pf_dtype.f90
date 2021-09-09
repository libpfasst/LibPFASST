!!  Data types and interfaces
!
! This file is part of LIBPFASST.
!
!>  Module to define the main parameters, data types, and interfaces in pfasst
module pf_mod_dtype
  use iso_c_binding
  implicit none

  
  !>  pfasst static  paramters
  integer, parameter :: pfdp = selected_real_kind(15, 307)  !!  Defines double precision type for all real and complex variables
!  integer, parameter :: pfdp = selected_real_kind(33, 4931)  !! For quad precision everywhere (use at your risk and see top of pf_mpi.f90)
  integer, parameter :: pfqp = selected_real_kind(33, 4931) !!  Defines quad precision type for all real and complex variableso
  real(pfdp), parameter :: ZERO  = 0.0_pfdp
  real(pfdp), parameter :: ONE   = 1.0_pfdp
  real(pfdp), parameter :: TWO   = 2.0_pfdp
  real(pfdp), parameter :: THREE  = 3.0_pfdp
  real(pfdp), parameter :: HALF  = 0.5_pfdp
  complex(pfdp), parameter :: ZI  = cmplx(0.0,1.0,pfdp)
  complex(pfdp), parameter :: Z0  = cmplx(0.0,0.0,pfdp)  
  real(pfqp),parameter ::  qpi = 3.1415926535897932384626433832795_pfqp
  real(pfdp),parameter ::  two_pi = 2.0_pfqp*qpi
  integer, parameter :: PF_MAXLEVS = 4
  integer, parameter :: PF_MAX_HOOKS = 32

  !> Quadrature node varieties
  integer, parameter :: SDC_GAUSS_LOBATTO   = 1
  integer, parameter :: SDC_GAUSS_RADAU  = 2
  integer, parameter :: SDC_CLENSHAW_CURTIS = 3
  integer, parameter :: SDC_UNIFORM         = 4
  integer, parameter :: SDC_GAUSS_LEGENDRE  = 5
  integer, parameter :: SDC_CHEBYSHEV       = 6

  !> States of operation
  integer, parameter :: PF_STATUS_ITERATING = 1
  integer, parameter :: PF_STATUS_CONVERGED = 2
  integer, parameter :: PF_STATUS_PREDICTOR = 3

  !>  Type for storing timings  later output
  integer, parameter :: PF_NUM_TIMERS = 22
  type :: pf_timer_t
     real(pfdp) :: timers(PF_NUM_TIMERS,PF_MAXLEVS)=0.0d0
     real(pfdp) :: runtimes(PF_NUM_TIMERS,PF_MAXLEVS)=0.0d0
  end type pf_timer_t
  
  !>  The type that holds the state of the system
  type, bind(c) :: pf_state_t
    real(pfdp) :: t0  !!  Time at beginning of this time step
    real(pfdp) :: dt  !!  Time step size
    integer :: nsteps   !! total number of time steps
    integer :: pfblock  !! pfasst block being worked on
    integer :: iter     !! current iteration number
    integer :: step     !! current time step number assigned to processor
    integer :: level_index  !! which level is currently being operated on
    integer :: finest_level !! the current finest level (for variable depth V cycles)
    integer :: hook     !! which hook
    integer :: proc     !! which processor
    integer :: sweep    !! sweep number
    integer :: status   !! status (iterating, converged etc)
    integer :: pstatus  !! previous rank's status
    logical :: pconverged  !! is previous rank converged
    integer :: itcnt    !! total iterations by this processor
    integer :: skippedy !! skipped sweeps for state (for mixed integration)
    integer :: mysteps  !! steps I did
  end type pf_state_t

  !>  Abstract hook type: hooks call diagnostic routines from various places in code
  type :: pf_hook_t
     procedure(pf_hook_p), pointer, nopass :: proc
  end type pf_hook_t

  !>  The abstract SDC sweeper type (must be extended)
  type, abstract :: pf_sweeper_t
     integer     :: npieces
     logical     :: use_LUq
   contains
     procedure(pf_sweep_p),        deferred :: sweep
     procedure(pf_initialize_p),   deferred :: initialize
     procedure(pf_evaluate_p),     deferred :: evaluate
     procedure(pf_integrate_p),    deferred :: integrate
     procedure(pf_evaluate_all_p), deferred :: evaluate_all
     procedure(pf_residual_p),     deferred :: residual
     procedure(pf_spreadq0_p),     deferred :: spreadq0
     procedure(pf_compute_dt_p),   deferred :: compute_dt
     procedure(pf_destroy_sweeper_p),      deferred :: destroy
  end type pf_sweeper_t

  !>  The abstract time stepper type (must be extended)
  type, abstract :: pf_stepper_t
     integer     :: npieces !  Number of pieces of rhs
     integer     :: order   !  Method order
     integer     :: nsteps  !  Number of steps per big time step
     
   contains
     procedure(pf_do_n_steps_p),           deferred :: do_n_steps
     procedure(pf_initialize_stepper_p),   deferred :: initialize
     procedure(pf_destroy_stepper_p),      deferred :: destroy
  end type pf_stepper_t

  !>  The abstract data type of the solution (must be extended)
  type, abstract :: pf_encap_t
   contains
     procedure(pf_encap_setval_p),  deferred :: setval
     procedure(pf_encap_copy_p),    deferred :: copy
     procedure(pf_encap_norm_p),    deferred :: norm
     procedure(pf_encap_pack_p),    deferred :: pack
     procedure(pf_encap_unpack_p),  deferred :: unpack
     procedure(pf_encap_axpy_p),    deferred :: axpy
     procedure(pf_encap_eprint_p),  deferred :: eprint
  end type pf_encap_t

  !>  Abstract type for creation and destruction of objects
  type, abstract :: pf_factory_t
   contains
     procedure(pf_encap_create_single_p),  deferred :: create_single
     procedure(pf_encap_create_array_p),   deferred :: create_array
     procedure(pf_encap_destroy_single_p), deferred :: destroy_single
     procedure(pf_encap_destroy_array_p),  deferred :: destroy_array
  end type pf_factory_t

  !>  The absract definition of level which is inherited  to include problem dependent stuff
  type, abstract :: pf_user_level_t
     class(pf_factory_t), allocatable :: factory
     class(pf_sweeper_t), allocatable :: sweeper
     class(pf_stepper_t), allocatable :: stepper
   contains
     procedure(pf_transfer_p), deferred :: restrict
     procedure(pf_transfer_p), deferred :: interpolate
  end type pf_user_level_t

  !>  The type to store quadrature matrices
  type :: pf_sdcmats_t
     integer :: nnodes                      !  Number of nodes
     integer :: qtype                   !  Type of nodes
     real(pfdp), allocatable :: qnodes(:)   !  The quadrature nodes
     real(pfdp), allocatable :: Qmat(:,:)   !  Collocation matrix
     real(pfdp), allocatable :: QmatFE(:,:) !  Forward Euler matrix
     real(pfdp), allocatable :: QmatBE(:,:) !  Backward Euler matrix
     real(pfdp), allocatable :: QmatTrap(:,:) ! Trapezoid rule matrix
     real(pfdp), allocatable :: QmatVer(:,:)  ! Verlet Matrix
     real(pfdp), allocatable :: QmatLU(:,:)   !  LU of Wmat
     real(pfdp), allocatable :: Smat(:,:)    !  The node to node version of Qmat

     logical :: use_proper_nodes =  .false. !  If true use gauss nodes in coarsening
     logical :: use_composite_nodes = .false. ! If true, finer nodes are composite
     logical :: use_no_left_q = .false.       ! If true don't use left endpoint in rule
  end type pf_sdcmats_t


  !>  Data type of a PFASST level
  type :: pf_level_t
     !  ===Mandatory level parameters===
     integer  :: mpibuflen    = -1   !! size of solution in pfdp units

     !  level parameters set by the pfasst_t values
     integer  :: index        = -1   !! level number (1 is the coarsest)
     integer  :: nnodes       = -1   !! number of sdc nodes
     integer  :: nsweeps      = -1   !! number of sdc sweeps to perform
     integer  :: nsweeps_pred = -1      !! number of coarse sdc sweeps to perform predictor in predictor
     logical     :: Finterp = .false.   !! interpolate functions instead of solutions


     !  Diagnostics
     real(pfdp)  :: error            !! holds the user defined error
     real(pfdp)  :: residual         !! holds the user defined residual
     real(pfdp)  :: residual_rel     !! holds the user defined relative residual (scaled by solution magnitude)

     class(pf_user_level_t), allocatable :: ulevel  !!  user customized level info

     !>  Simple data storage at each level
     real(pfdp), allocatable :: &
          send(:),    &                 !! send buffer
          recv(:),    &                 !! recv buffer
          nodes(:),   &                 !! list of SDC nodes
          rmat(:,:),  &                 !! time restriction matrix
          tmat(:,:)                     !! time interpolation matrix

     integer, allocatable :: &
          nflags(:)                     !! sdc node flags

     !>  Solution variable storage
     class(pf_encap_t), allocatable :: &
          Q(:),     &           !! solution at sdc nodes
          pQ(:),    &           !! unknowns at sdc nodes, previous sweep
          R(:),     &           !! full residuals
          I(:),     &           !! 0 to node integrals
          Fflt(:),  &           !! functions values at sdc nodes (flat)
          Frkflt(:), &          !!  Stage Function values
          tauQ(:),  &           !! fas correction in Q form
          pFflt(:), &           !! functions at sdc nodes, previous sweep (flat)
          q0,       &           !! initial condition
          delta_q0, &           !! Space for interpolating q0, qend
          qend                  !! solution at end time
     real(pfdp) :: max_delta_q0=0.0_pfdp
     !>  Function  storage
     class(pf_encap_t), pointer :: &
          F(:,:), &                     !! functions values at sdc nodes
          pF(:,:)                       !! functions at sdc nodes, previous sweep

     !>  Interpolation and restriction data structures
     logical :: interp_workspace_allocated = .false.
     logical :: restrict_workspace_allocated = .false.
     class(pf_encap_t), allocatable :: &
          cf_delta(:), &   ! delta fine in space and coarse in time
          c_delta(:), &    ! delta on the coarse level
          f_encap_array_c(:)  !  fine solution restricted in space only
     
     integer, allocatable :: lev_shape(:)   !! user defined shape array
     type(pf_sdcmats_t), allocatable :: sdcmats
     logical :: allocated = .false.
  end type pf_level_t

  !>  Data type to define the communicator
  type :: pf_comm_t
     integer :: nproc = -1              ! total number of processors

     integer :: comm = -1               ! communicator
     integer, pointer :: &
          recvreq(:), &                 ! receive requests (indexed by level)
          sendreq(:)                    ! send requests (indexed by level)
     integer :: statreq                 ! status send request

     ! fakie, needs modernization
     !type(c_ptr), pointer :: pfs(:)     ! pfasst objects (indexed by rank)
     !type(c_ptr), pointer :: pfpth(:,:)

     !> Procedure interfaces
     procedure(pf_post_p),        pointer, nopass :: post
     procedure(pf_recv_p),        pointer, nopass :: recv
     procedure(pf_recv_status_p), pointer, nopass :: recv_status
     procedure(pf_send_p),        pointer, nopass :: send
     procedure(pf_send_status_p), pointer, nopass :: send_status
     procedure(pf_wait_p),        pointer, nopass :: wait
     procedure(pf_broadcast_p),   pointer, nopass :: broadcast
  end type pf_comm_t

  type :: pf_results_t
     real(pfdp), allocatable ::    errors(:,:,:,:)
     real(pfdp), allocatable :: residuals(:,:,:,:)  !  (level,block,niter+1,sweep)
     real(pfdp), allocatable ::  delta_q0(:,:,:,:)  !  (level,block,niter+1,sweep)
     real(pfdp), allocatable ::  iters(:)  !           (block)
     integer :: nlevs
     integer :: nsteps
     integer :: niters  !  really the max niters
     integer :: nprocs  
     integer :: nblocks
     integer :: max_nsweeps  !  max nsweeps for allocation
     integer :: rank

     logical :: save_residuals
     logical :: save_errors
     logical :: save_delta_q0

     character(len=128) :: datpath
     procedure(pf_results_p), pointer, nopass :: destroy 
  end type pf_results_t

  !>  The main PFASST data type which includes pretty much everythingl
  type :: pf_pfasst_t
     !> === Mandatory pfasst parameters (must be set on command line or input file)  ===
     integer :: nlevels = -1             !! number of pfasst levels

     !>  ===  Optional pfasst parameters ====
     integer :: niters  = 5             !! number of PFASST iterations to do
     integer :: MINiters  = 0           !! MIN number of PFASST iterations to do
     integer :: qtype   = SDC_GAUSS_LOBATTO  !! type of nodes
     logical :: use_proper_nodes =  .false.
     logical :: use_composite_nodes = .false.
     logical :: use_no_left_q = .false.

     ! --  level dependent parameters
     integer :: nsweeps(PF_MAXLEVS) = 1       !!  number of sweeps at each levels
     integer :: nsweeps_pred(PF_MAXLEVS) =1   !!  number of sweeps during predictor
     integer :: nnodes(PF_MAXLEVS)=3          !! number of nodes

     ! --  tolerances
     real(pfdp) :: abs_res_tol = 0.0_pfdp   !!  absolute convergence tolerance
     real(pfdp) :: rel_res_tol = 0.0_pfdp   !!  relative convergence tolerance

     ! --  predictor options  (should be set before pfasst_run is called)
     logical :: PFASST_pred = .true.    !!  true if the PFASST type predictor is used
     logical :: pipeline_pred = .false. !!  true if coarse sweeps after burn in are pipelined  (if nsweeps_pred>1 on coarse level)
     integer :: nsweeps_burn =  1       !!  number of sdc sweeps to perform during coarse level burn in
     integer :: q0_style =  0           !!  q0 can take 3 values
                                        !!  0:  Only the q0 at t=0 is valid  (default)
                                        !!  1:  The q0 at each processor is valid
                                        !!  2:  q0 and all nodes at each processor is valid


     ! --  run options  (should be set before pfasst_run is called)
     logical :: Vcycle = .true.         !!  decides if Vcycles are done
     logical :: use_pysdc_V = .false.         !!  decides if Vcycles are done
     logical :: sweep_at_conv = .false. !!  decides if one final sweep after convergence is done
     logical :: Finterp = .false.    !!  True if transfer functions operate on rhs
     logical :: use_LUq = .true.     !!  True if LU type implicit matrix is used
     logical :: use_Sform = .false.  !!  True if Qmat type of stepping is used
     integer :: taui0 = -99          !! iteration cutoff for tau inclusion


     ! -- RK and Parareal options
     logical :: use_sdc_sweeper =.true.   !! decides if SDC sweeper is used 
     logical :: use_rk_stepper = .false.  !! decides if RK steps are used instead of the sweeps
     integer :: nsteps_rk(PF_MAXLEVS)=-1  !! number of runge-kutta steps per time step
     integer :: rk_order(PF_MAXLEVS)=-1   !! order of runge-kutta method per level
     integer :: rk_nstages(PF_MAXLEVS)=-1 !! number of runge-kutta stages per level
     logical :: RK_pred = .false.         !!  true if the coarse level is initialized with Runge-Kutta instead of PFASST

     ! -- misc
     logical :: debug = .false.         !!  If true, debug diagnostics are printed

     ! -- controller for the results 
     logical :: save_residuals = .true.  !!  Will save residuals every time they are set
     logical :: save_delta_q0 = .true.   !!  Will save change in initial condition
     logical :: save_errors  = .true.    !!  Will save errors, but set_error must be called externally
     logical :: save_json = .true.       !!  Will save a jason file of run parameters
     integer :: save_timings  = 2        !!  0=none, 1=total only, 2=all, 3=all and echo
     integer :: save_solutions  = 0      !!  0=none, 1=end, 2=all time steps, 3=all iterations

     integer :: rank    = -1            !! rank of current processor

     !> pf objects
     type(pf_state_t), allocatable :: state   !!  Describes where in the algorithm  is
     type(pf_level_t), allocatable :: levels(:) !! Holds the levels
     type(pf_comm_t),  pointer :: comm    !! Points to communicator
     type(pf_results_t) :: results   !!  Hold results for each level
 
     !> hooks variables
     type(pf_hook_t), allocatable :: hooks(:,:,:)  !!  Holds the hooks
     integer,  allocatable :: nhooks(:,:)   !!  Holds the number hooks

     !> timing variables
     type(pf_timer_t) :: pf_timers
     type(pf_timer_t) :: pf_runtimes
     !> output directory
     character(len=256) :: outdir

  end type pf_pfasst_t

  !> Interfaces for subroutines
  interface
    !> hooks subroutines
    subroutine pf_hook_p(pf, level_index)
       use iso_c_binding
       import pf_pfasst_t
       type(pf_pfasst_t), intent(inout) :: pf
       integer, intent(in) :: level_index
     end subroutine pf_hook_p

     !> SDC sweeper subroutines
     subroutine pf_sweep_p(this, pf, level_index, t0, dt, nsweeps, flags)
       import pf_sweeper_t,pf_pfasst_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
       real(pfdp),          intent(in)    :: t0
       real(pfdp),          intent(in)    :: dt
       integer,             intent(in)    :: nsweeps
       integer, optional,   intent(in)    :: flags
     end subroutine pf_sweep_p

     subroutine pf_evaluate_p(this, pf, level_index, t, m, flags, step)
       import pf_sweeper_t,pf_pfasst_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
       real(pfdp),          intent(in)    :: t
       integer,             intent(in)    :: m
       integer, optional,   intent(in)    :: flags, step
     end subroutine pf_evaluate_p

     subroutine pf_evaluate_all_p(this, pf, level_index, t, flags, step)
       import pf_sweeper_t, pf_pfasst_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
       real(pfdp),          intent(in)    :: t(:)
       integer, optional,   intent(in)    :: flags, step
     end subroutine pf_evaluate_all_p

     subroutine pf_initialize_p(this, pf, level_index)
       import pf_sweeper_t, pf_pfasst_t
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
     end subroutine pf_initialize_p

     subroutine pf_destroy_sweeper_p(this,pf,level_index)
       import pf_sweeper_t, pf_pfasst_t
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
     end subroutine pf_destroy_sweeper_p

     subroutine pf_integrate_p(this, pf, level_index, qSDC, fSDC, dt, fintSDC, flags)
       import pf_sweeper_t, pf_pfasst_t, pf_encap_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
       class(pf_encap_t),   intent(in)    :: qSDC(:), fSDC(:, :)
       real(pfdp),          intent(in)    :: dt !!  Time step size
       class(pf_encap_t),   intent(inout) :: fintSDC(:)
       integer, optional,   intent(in)    :: flags
     end subroutine pf_integrate_p

     subroutine pf_residual_p(this, pf, level_index, dt, flags)
       import  pf_sweeper_t,pf_pfasst_t,  pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),  target,  intent(inout) :: pf
       integer,              intent(in)    :: level_index
       real(pfdp),          intent(in)    :: dt !!  Time step size
       integer, optional,   intent(in)    :: flags
     end subroutine pf_residual_p

     subroutine pf_spreadq0_p(this, pf, level_index, t0, flags, step)
       import pf_sweeper_t, pf_pfasst_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),  target,  intent(inout) :: pf
       integer,              intent(in)    :: level_index
       real(pfdp),          intent(in)    :: t0 !!  Time at beginning of step; if flags == 2, time at end of step
       integer, optional,   intent(in)    :: flags, step
     end subroutine pf_spreadq0_p

     subroutine pf_compute_dt_p(this, pf, level_index, t0, dt, flags)
       import pf_sweeper_t,pf_pfasst_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
       real(pfdp),          intent(in)    :: t0
       real(pfdp),          intent(inout)    :: dt
       integer, optional,   intent(in)    :: flags
     end subroutine pf_compute_dt_p


     subroutine pf_destroy_p(this, pf,level_index)
       import pf_sweeper_t,pf_pfasst_t, pfdp
       class(pf_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),  target,  intent(inout) :: pf
       integer,              intent(in)    :: level_index
     end subroutine pf_destroy_p

     !>  time stepper interfaces
     subroutine pf_do_n_steps_p(this, pf, level_index, t0,y0,yend, big_dt,nsteps_rk)
       import pf_pfasst_t, pf_stepper_t, pf_level_t, pfdp, pf_encap_t
       class(pf_stepper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       real(pfdp),          intent(in)    :: big_dt !!  Time step size
       real(pfdp),          intent(in)    :: t0
       class(pf_encap_t), intent(in   )         :: y0           !!  Starting value
       class(pf_encap_t), intent(inout)         :: yend         !!  Final value
       integer,             intent(in)    :: level_index
       integer,             intent(in)    :: nsteps_rk
     end subroutine pf_do_n_steps_p

     subroutine pf_initialize_stepper_p(this,pf,level_index)
       import pf_stepper_t, pf_pfasst_t
       class(pf_stepper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
     end subroutine pf_initialize_stepper_p

     subroutine pf_destroy_stepper_p(this,pf,level_index)
       import pf_stepper_t, pf_pfasst_t
       class(pf_stepper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index
     end subroutine pf_destroy_stepper_p

     !> transfer interfaces used for restriction and interpolation
     subroutine pf_transfer_p(this, f_lev, c_lev, f_vec, c_vec, t, flags)
       import pf_user_level_t, pf_level_t, pf_encap_t, pfdp
       class(pf_user_level_t), intent(inout) :: this
       class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
       class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
       real(pfdp),          intent(in)       :: t
       integer, optional,   intent(in)       :: flags
     end subroutine pf_transfer_p

     !> encapsulation interfaces
     subroutine pf_encap_create_single_p(this, x, level_index, lev_shape)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x
       integer,      intent(in   )              :: level_index,  lev_shape(:)
     end subroutine pf_encap_create_single_p

     subroutine pf_encap_create_array_p(this, x, n, level_index, lev_shape)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x(:)
       integer,      intent(in   )              :: n, level_index, lev_shape(:)
     end subroutine pf_encap_create_array_p

     subroutine pf_encap_destroy_single_p(this, x)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x
     end subroutine pf_encap_destroy_single_p

     subroutine pf_encap_destroy_array_p(this, x)
       import pf_factory_t, pf_encap_t
       class(pf_factory_t), intent(inout)              :: this
       class(pf_encap_t),   intent(inout), allocatable :: x(:)
     end subroutine pf_encap_destroy_array_p

     subroutine pf_encap_setval_p(this, val, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout)        :: this
       real(pfdp),        intent(in)           :: val
       integer,    intent(in), optional :: flags
     end subroutine pf_encap_setval_p

     subroutine pf_encap_copy_p(this, src, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout)           :: this
       class(pf_encap_t), intent(in   )           :: src
       integer,    intent(in   ), optional :: flags
     end subroutine pf_encap_copy_p

     function pf_encap_norm_p(this, flags) result (norm)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(in   ) :: this
       integer,    intent(in   ), optional :: flags
       real(pfdp) :: norm
     end function pf_encap_norm_p

     subroutine pf_encap_pack_p(this, z, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(in   ) :: this
       real(pfdp),        intent(  out) :: z(:)
       integer, optional,   intent(in)  :: flags
     end subroutine pf_encap_pack_p

     subroutine pf_encap_unpack_p(this, z, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout) :: this
       real(pfdp),        intent(in   ) :: z(:)
       integer, optional,   intent(in)  :: flags
     end subroutine pf_encap_unpack_p

     subroutine pf_encap_axpy_p(this, a, x, flags)
       import pf_encap_t, pfdp
       class(pf_encap_t), intent(inout)  :: this
       class(pf_encap_t), intent(in   )  :: x
       real(pfdp),  intent(in)           :: a
       integer, intent(in), optional :: flags
     end subroutine pf_encap_axpy_p

     subroutine pf_encap_eprint_p(this,flags)
       import pf_encap_t
       class(pf_encap_t), intent(inout) :: this
       integer, intent(in), optional :: flags
     end subroutine pf_encap_eprint_p

     !> communicator interfaces
     subroutine pf_post_p(pf, level, tag, ierror, source)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(in)    :: pf
       class(pf_level_t), intent(inout) :: level
       integer,    intent(in)           :: tag
       integer,    intent(inout)        :: ierror
       integer,    intent(in)           :: source
     end subroutine pf_post_p

     subroutine pf_recv_p(pf, level, tag, blocking, ierror, source)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       class(pf_level_t), intent(inout) :: level
       integer,    intent(in)    :: tag
       logical,           intent(in)    :: blocking
       integer,    intent(inout)       :: ierror
       integer,          intent(in)    :: source
     end subroutine pf_recv_p

     subroutine pf_recv_status_p(pf, tag,istatus,ierror, source)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       integer,    intent(in)         :: tag
       integer,    intent(inout)      :: istatus
       integer,    intent(inout)      :: ierror
       integer,     intent(in)        :: source
     end subroutine pf_recv_status_p

     subroutine pf_send_p(pf, level, tag, blocking,ierror, dest)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       class(pf_level_t), intent(inout) :: level
       integer,    intent(in)    :: tag
       logical,           intent(in)    :: blocking
       integer,    intent(inout)       :: ierror
       integer,             intent(in)    :: dest
     end subroutine pf_send_p

     subroutine pf_send_status_p(pf, tag,istatus,ierror, dest)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       integer,    intent(in)        :: tag
       integer,    intent(in)        :: istatus
       integer,    intent(inout)     :: ierror
       integer,    intent(in)        :: dest
     end subroutine pf_send_status_p

     subroutine pf_wait_p(pf, level,ierror)
       import pf_pfasst_t
       type(pf_pfasst_t), intent(in) :: pf
       integer,    intent(in) :: level
       integer,    intent(inout)       :: ierror
     end subroutine pf_wait_p

     subroutine pf_broadcast_p(pf, y, nvar, root,ierror)
       import pf_pfasst_t, pfdp
       type(pf_pfasst_t), intent(inout) :: pf
       integer,    intent(in)    :: nvar, root
       real(pfdp)  ,      intent(in)    :: y(nvar)
       integer,    intent(inout)       :: ierror
     end subroutine pf_broadcast_p


     subroutine pf_results_p(this)
       import pf_results_t
       type(pf_results_t), intent(inout) :: this

     end subroutine pf_results_p
     subroutine pf_resultsp_p(this,pf)
       import pf_results_t,pf_pfasst_t
       type(pf_results_t), intent(inout) :: this
       type(pf_pfasst_t), intent(inout) :: pf
     end subroutine pf_resultsp_p



    end interface

end module pf_mod_dtype
