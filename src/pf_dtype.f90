!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module pf_mod_dtype
  use iso_c_binding
  use encap
  implicit none


  integer, parameter :: PF_MAX_HOOKS = 8

  integer, parameter :: SDC_GAUSS_LOBATTO   = 1
  integer, parameter :: SDC_GAUSS_RADAU     = 2
  integer, parameter :: SDC_CLENSHAW_CURTIS = 3
  integer, parameter :: SDC_UNIFORM         = 4
  integer, parameter :: SDC_GAUSS_LEGENDRE  = 5

  ! state type
  type :: pf_state_t
     real(pfdp) :: t0, dt
     integer    :: block, cycle, step, iter, nsteps
  end type pf_state_t


  ! hook type
  type :: pf_hook_t
     integer :: hook = -1               ! hook type (see pf_mod_hooks)
     procedure(hook_proc), pointer, nopass :: proc
  end type pf_hook_t


  ! level type
  type :: pf_level_t
     integer     :: nvars = -1          ! number of variables (dofs)
     integer     :: nnodes = -1         ! number of sdc nodes
     integer     :: nsweeps = 1         ! number of sdc sweeps to perform
     integer     :: level = -1          ! level number (1=finest)
     logical     :: Finterp = .false.   ! Interpolate functions instead of solution
     ! tolerances
     real(pfdp)  :: residual_tol = 0.0_pfdp

     ! arrays
     real(pfdp), pointer :: &
          q0(:), &                      ! initial condition (packed)
          send(:), &                    ! send buffer
          recv(:), &                    ! recv buffer
          nodes(:), &                   ! sdc nodes
          qmat(:,:), &                  ! integration matrix (0 to node)
          s0mat(:,:), &                 ! integration matrix (node to node)
          smat(:,:,:) => null(), &      ! sdc matrices (allocated by the sweeper)
          rmat(:,:) => null(), &        ! time restriction matrix
          tmat(:,:) => null()           ! time interpolation matrix

     logical, pointer :: nmask(:)       ! sdc node mask

     type(pf_encap_t) :: qend           ! end value
     type(pf_encap_t) :: qex            ! exact value for diagnostics

     type(pf_encap_t), pointer :: &
          qSDC(:) => null(), &          ! unknowns at sdc nodes
          pSDC(:) => null(), &          ! unknowns at sdc nodes, previous sweep
          fSDC(:,:) => null(), &        ! functions values at sdc nodes
          pfSDC(:,:) => null(), &       ! functions at sdc nodes, previous sweep
          tau(:) => null()              ! fas correction

     type(c_ptr) :: ctx  = c_null_ptr   ! user context
     type(c_ptr) :: dctx = c_null_ptr   ! dump context

     integer, pointer :: shape(:) => null() ! user shape

     ! rhs
     procedure(imp_rhs_proc), pointer, nopass  :: gen_imp_rhs => null()
     procedure(imex_rhs_proc), pointer, nopass :: gen_imex_rhs => null()

     logical :: allocated = .false.
  end type pf_level_t


  ! pfasst communicator
  type :: pf_comm_t
     integer :: nproc = -1              ! total number of processors

     ! mpi
     integer :: comm = -1               ! communicator
     integer, pointer :: &
          recvreq(:), &                 ! receive requests (indexed by level)
          sendreq(:)                    ! send requests (indexed by level)

     ! pthreads
     type(c_ptr), pointer :: pfs(:)     ! pfasst objects (indexed by rank)
     type(c_ptr), pointer :: pfpth(:,:) ! mutexes and conditions (indexed by rank, level)
  end type pf_comm_t


  ! pfasst type
  type :: pf_pfasst_t
     integer :: nlevels = -1            ! number of pfasst levels
     integer :: niters  = 5             ! number of iterations
     integer :: qtype   = 1             ! type of quadrature nodes
     integer :: rank    = -1            ! rank of current processor
     real(pfdp)        :: Htol          ! Hamiltonian tol in Verlet
     ! pf objects
     type(pf_state_t)          :: state
     type(pf_level_t), pointer :: levels(:)
     type(pf_comm_t), pointer  :: comm

     ! hooks
     type(pf_hook_t), pointer :: hooks(:,:)
     integer, pointer         :: nhooks(:)

     ! logging, timer
     logical :: echo_timings = .false.
     integer(8) :: timers(100)   = 0
     integer(8) :: runtimes(100) = 0
     integer :: log = -1 ! log file unit number
  end type pf_pfasst_t


  ! hook interface
  interface
     subroutine hook_proc(pf, level, state, ctx)
       use iso_c_binding
       import pf_pfasst_t
       import pf_level_t
       import pf_state_t
       type(pf_pfasst_t), intent(inout) :: pf
       type(pf_level_t),  intent(inout) :: level
       type(pf_state_t),  intent(in)    :: state
       type(c_ptr),       intent(in)    :: ctx
     end subroutine hook_proc
  end interface


  ! gen rhs interfaces
  interface
     subroutine imex_rhs_proc(rhs, q0, dt, f1, S, level, ctx)
       use iso_c_binding
       use encap
       type(pf_encap_t), intent(inout) :: rhs
       type(pf_encap_t), intent(in)    :: q0, f1, S
       real(pfdp),       intent(in)    :: dt
       integer,          intent(in)    :: level
       type(c_ptr),      intent(in)    :: ctx
     end subroutine imex_rhs_proc
  end interface

  interface
     subroutine imp_rhs_proc(rhs, q0, S, level, ctx)
       use iso_c_binding
       use encap
       type(pf_encap_t), intent(inout) :: rhs
       type(pf_encap_t), intent(in)    :: q0, S
       integer,          intent(in)    :: level
       type(c_ptr),      intent(in)    :: ctx
     end subroutine imp_rhs_proc
  end interface

end module pf_mod_dtype
