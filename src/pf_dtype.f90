!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
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

module pf_mod_dtype
  use iso_c_binding
  implicit none

  !  static pfasst paramters
!  integer, parameter :: pfdp = c_long_double
  integer, parameter :: pfdp = c_double


  real(pfdp), parameter :: ZERO  = 0.0_pfdp
  real(pfdp), parameter :: ONE   = 1.0_pfdp
  real(pfdp), parameter :: TWO   = 2.0_pfdp
  real(pfdp), parameter :: HALF  = 0.5_pfdp

  integer, parameter :: PF_MAX_HOOKS = 32

  integer, parameter :: SDC_GAUSS_LOBATTO   = 1
  integer, parameter :: SDC_GAUSS_RADAU     = 2
  integer, parameter :: SDC_CLENSHAW_CURTIS = 3
  integer, parameter :: SDC_UNIFORM         = 4
  integer, parameter :: SDC_GAUSS_LEGENDRE  = 5
  integer, parameter :: SDC_PROPER_NODES    = 2**8
  integer, parameter :: SDC_COMPOSITE_NODES = 2**9
  integer, parameter :: SDC_NO_LEFT         = 2**10

  integer, parameter :: SDC_KIND_SOL_FEVAL    = 1
  integer, parameter :: SDC_KIND_SOL_NO_FEVAL = 2
  integer, parameter :: SDC_KIND_FEVAL        = 3
  integer, parameter :: SDC_KIND_INTEGRAL     = 4
  integer, parameter :: SDC_KIND_CORRECTION   = 5

  integer, parameter :: PF_WINDOW_BLOCK = 1
  integer, parameter :: PF_WINDOW_RING  = 2
  integer, parameter :: PF_TAG_NMOVED   = 666

  integer, parameter :: PF_STATUS_ITERATING = 1
  integer, parameter :: PF_STATUS_CONVERGED = 2
  integer, parameter :: PF_STATUS_PREDICTOR = 3


  type, bind(c) :: pf_state_t
     real(pfdp) :: t0, dt
!     real(c_double) :: t0, dt
     integer(c_int) :: nsteps, block, cycle, step, iter, level, hook, proc
     integer(c_int) :: status       ! status (iterating, converged etc)
     integer(c_int) :: pstatus      ! previous rank's status
     integer(c_int) :: nmoved       ! how many processors behind me have moved
     integer(c_int) :: first        ! rank of first processor in time block
     integer(c_int) :: last         ! rank of last processor in time block
     integer(c_int) :: itcnt        ! iteration counter
     integer(c_int) :: mysteps      ! steps I did
!     real(c_double) :: res
     real(pfdp) :: res
  end type pf_state_t

  type :: pf_hook_t
     procedure(pf_hook_p), pointer, nopass :: proc
  end type pf_hook_t

  type :: pf_sweeper_t
     type(c_ptr) :: sweeperctx
     integer     :: npieces
     procedure(pf_sweep_p),        pointer, nopass :: sweep
     procedure(pf_initialize_p),   pointer, nopass :: initialize
     procedure(pf_evaluate_p),     pointer, nopass :: evaluate
     procedure(pf_evaluate_all_p), pointer, nopass :: evaluate_all
     procedure(pf_integrate_p),    pointer, nopass :: integrate
     procedure(pf_residual_p),     pointer, nopass :: residual
     procedure(pf_sweepdestroy_p), pointer, nopass :: destroy
  end type pf_sweeper_t

  type :: pf_encap_t
     procedure(pf_encap_create_p),  pointer, nopass :: create
     procedure(pf_encap_destroy_p), pointer, nopass :: destroy
     procedure(pf_encap_setval_p),  pointer, nopass :: setval
     procedure(pf_encap_printme_p),  pointer, nopass :: printme
     procedure(pf_encap_copy_p),    pointer, nopass :: copy
     procedure(pf_encap_norm_p),    pointer, nopass :: norm
     procedure(pf_encap_pack_p),    pointer, nopass :: pack
     procedure(pf_encap_unpack_p),  pointer, nopass :: unpack
     procedure(pf_encap_axpy_p),    pointer, nopass :: axpy
     procedure(pf_encap_eprint_p),    pointer, nopass :: eprint
  end type pf_encap_t

  type :: pf_level_t
     integer     :: nvars = -1          ! number of variables (dofs)
     integer     :: nnodes = -1         ! number of sdc nodes
     integer     :: nsweeps = 1         ! number of sdc sweeps to perform
     integer     :: nsweeps_pred = 1         ! number of sdc sweeps to perform
     integer     :: level = -1          ! level number (1 is the coarsest)
     logical     :: Finterp = .false.   ! interpolate functions instead of solutions


     real(pfdp)  :: residual

     type(pf_encap_t),         pointer         :: encap
     type(pf_sweeper_t),       pointer         :: sweeper
     procedure(pf_transfer_p), pointer, nopass :: interpolate, restrict

     real(pfdp), pointer :: &
          send(:), &                    ! send buffer
          recv(:), &                    ! recv buffer
          nodes(:), &                   ! sdc nodes
          qmat(:,:), &                  ! integration matrix (0 to node)
          s0mat(:,:), &                 ! integration matrix (node to node)
          rmat(:,:), &                  ! time restriction matrix
          tmat(:,:)                     ! time interpolation matrix

     integer(c_int), pointer :: &
          nflags(:)                     ! sdc node flags

     type(c_ptr), pointer :: &
          q0(:), &                      ! initial condition encap
          Q(:), &                       ! unknowns at sdc nodes
          pQ(:), &                      ! unknowns at sdc nodes, previous sweep
          F(:,:), &                     ! functions values at sdc nodes
          pF(:,:), &                    ! functions at sdc nodes, previous sweep
          R(:), &                       ! full residuals
          I(:), &                       ! 0 to node integrals
          S(:), &                       ! node to node integrals
          tau(:), &                     ! fas correction
          tauQ(:)                       ! fas correction in Q form

     type(c_ptr) :: qend                ! solution at last node

     type(c_ptr)      :: ctx       ! user context
     integer, pointer :: shape(:)       ! user shape

     logical :: allocated = .false.
  end type pf_level_t

  type :: pf_comm_t
     integer :: nproc = -1              ! total number of processors

     ! mpi
     integer :: comm = -1               ! communicator
     integer, pointer :: &
          recvreq(:), &                 ! receive requests (indexed by level)
          sendreq(:)                    ! send requests (indexed by level)
     integer :: statreq                 ! status send request

     ! fake
     type(c_ptr), pointer :: pfs(:)     ! pfasst objects (indexed by rank)
     type(c_ptr), pointer :: pfpth(:,:) ! mutexes and conditions (indexed by rank, level)

     procedure(pf_post_p),        pointer, nopass :: post
     procedure(pf_recv_p),        pointer, nopass :: recv
     procedure(pf_recv_status_p), pointer, nopass :: recv_status
     procedure(pf_recv_status_p), pointer, nopass :: recv_nmoved
     procedure(pf_send_p),        pointer, nopass :: send
     procedure(pf_send_status_p), pointer, nopass :: send_status
     procedure(pf_send_status_p), pointer, nopass :: send_nmoved
     procedure(pf_wait_p),        pointer, nopass :: wait
     procedure(pf_broadcast_p),   pointer, nopass :: broadcast
  end type pf_comm_t

  type :: pf_pfasst_t
     integer :: nlevels = -1            ! number of pfasst levels
     integer :: niters  = 5             ! number of iterations
     integer :: rank    = -1            ! rank of current processor
     integer :: qtype   = SDC_GAUSS_LOBATTO

     real(pfdp) :: abs_res_tol = 0.d0
     real(pfdp) :: rel_res_tol = 0.d0

     integer :: window = PF_WINDOW_BLOCK

     logical :: Pipeline_G =  .false.
     logical :: PFASST_pred = .false.

     integer     :: taui0 = -999999     ! Cutoff for tau inclusion

     ! pf objects
     type(pf_state_t), pointer :: state
     type(pf_level_t), pointer :: levels(:)
     type(pf_comm_t),  pointer :: comm

     ! hooks
     type(pf_hook_t), pointer :: hooks(:,:,:)
     integer,         pointer :: nhooks(:,:)

     ! timing
     logical    :: echo_timings  = .false.
     integer(8) :: timers(100)   = 0
     integer(8) :: runtimes(100) = 0

     ! misc
     character(512) :: outdir
  end type pf_pfasst_t

  interface
     ! hook interface
     subroutine pf_hook_p(pf, level, state, ctx)
       use iso_c_binding
       import pf_pfasst_t, pf_level_t, pf_state_t
       type(pf_pfasst_t), intent(inout) :: pf
       type(pf_level_t),  intent(inout) :: level
       type(pf_state_t),  intent(in)    :: state
       type(c_ptr),       intent(in)    :: ctx
     end subroutine pf_hook_p

     ! sweeper interfaces
     subroutine pf_sweep_p(pf, Lev, t0, dt)
       import pf_pfasst_t, pf_level_t, pfdp, c_ptr
       type(pf_pfasst_t), intent(inout) :: pf
       real(pfdp),        intent(in)    :: dt, t0
       type(pf_level_t),  intent(inout) :: Lev
     end subroutine pf_sweep_p

     subroutine pf_evaluate_p(Lev, t, m)
       import pf_level_t, pfdp, c_ptr
       type(pf_level_t), intent(inout) :: Lev
       real(pfdp),       intent(in)    :: t
       integer,          intent(in)    :: m
     end subroutine pf_evaluate_p

     subroutine pf_evaluate_all_p(Lev, t)
       import pf_level_t, pfdp, c_ptr
       type(pf_level_t), intent(inout) :: Lev
       real(pfdp),       intent(in)    :: t(:)
     end subroutine pf_evaluate_all_p

     subroutine pf_initialize_p(Lev)
       import pf_level_t
       type(pf_level_t), intent(inout) :: Lev
     end subroutine pf_initialize_p

     subroutine pf_sweepdestroy_p(sweeper)
       import pf_sweeper_t
       type(pf_sweeper_t), intent(inout)  :: sweeper
     end subroutine pf_sweepdestroy_p

     subroutine pf_integrate_p(Lev, qSDC, fSDC, dt, fintSDC)
       import pf_level_t, c_ptr, pfdp
       type(pf_level_t),  intent(in)    :: Lev
       type(c_ptr),       intent(in)    :: qSDC(:), fSDC(:, :)
       real(pfdp),        intent(in)    :: dt
       type(c_ptr),       intent(inout) :: fintSDC(:)
     end subroutine pf_integrate_p

     subroutine pf_residual_p(Lev, dt)
       import pf_level_t, c_ptr, pfdp
       type(pf_level_t),  intent(inout) :: Lev
       real(pfdp),        intent(in)    :: dt
     end subroutine pf_residual_p

     ! transfer interfaces
     subroutine pf_transfer_p(qF, qG, levelF, ctxF, levelG, ctxG, t)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value :: qF, qG, ctxF, ctxG
       integer,     intent(in)        :: levelF, levelG
       real(pfdp),  intent(in)        :: t
     end subroutine pf_transfer_p

     ! encapsulation interfaces
     subroutine pf_encap_create_p(sol, level, kind, nvars, shape, ctx)
       import c_ptr
       type(c_ptr),  intent(inout)     :: sol
       type(c_ptr),  intent(in), value :: ctx
       integer,      intent(in)        :: level, nvars, shape(:)
       integer,      intent(in)        :: kind
     end subroutine pf_encap_create_p

     subroutine pf_encap_destroy_p(sol)
       import c_ptr
       type(c_ptr), intent(in), value :: sol
     end subroutine pf_encap_destroy_p

     subroutine pf_encap_setval_p(sol, val, flags)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value    :: sol
       real(pfdp),  intent(in)           :: val
       integer,     intent(in), optional :: flags
     end subroutine pf_encap_setval_p

     subroutine pf_encap_printme_p(sol)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value    :: sol

     end subroutine pf_encap_printme_p

     subroutine pf_encap_copy_p(dst, src, flags)
       import c_ptr
       type(c_ptr), intent(in), value    :: dst, src
       integer,     intent(in), optional :: flags
     end subroutine pf_encap_copy_p

     function pf_encap_norm_p(sol) result (norm)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value    :: sol
       real(pfdp) :: norm
     end function pf_encap_norm_p

     subroutine pf_encap_pack_p(z, q)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value :: q
       real(pfdp),  intent(out)       :: z(:)
     end subroutine pf_encap_pack_p

     subroutine pf_encap_unpack_p(q, z)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value :: q
       real(pfdp),  intent(in)        :: z(:)
     end subroutine pf_encap_unpack_p

     subroutine pf_encap_axpy_p(y, a, x, flags)
       import c_ptr, pfdp
       real(pfdp),  intent(in)           :: a
       type(c_ptr), intent(in), value    :: x, y
       integer,     intent(in), optional :: flags
     end subroutine pf_encap_axpy_p

     subroutine pf_encap_eprint_p(x)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value    :: x
     end subroutine pf_encap_eprint_p

     ! communicator interfaces
     subroutine pf_post_p(pf, level, tag)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(in)    :: pf
       type(pf_level_t),  intent(inout) :: level
       integer,           intent(in)    :: tag
     end subroutine pf_post_p

     subroutine pf_recv_p(pf, level, tag, blocking)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       type(pf_level_t),  intent(inout) :: level
       integer,           intent(in)    :: tag
       logical,           intent(in)    :: blocking
     end subroutine pf_recv_p

     subroutine pf_recv_status_p(pf, tag)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       integer,           intent(in)    :: tag
     end subroutine pf_recv_status_p

     subroutine pf_send_p(pf, level, tag, blocking)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       type(pf_level_t),  intent(inout) :: level
       integer,           intent(in)    :: tag
       logical,           intent(in)    :: blocking
     end subroutine pf_send_p

     subroutine pf_send_status_p(pf, tag)
       import pf_pfasst_t, pf_level_t
       type(pf_pfasst_t), intent(inout) :: pf
       integer,           intent(in)    :: tag
     end subroutine pf_send_status_p

     subroutine pf_wait_p(pf, level)
       import pf_pfasst_t
       type(pf_pfasst_t), intent(in) :: pf
       integer,           intent(in) :: level
     end subroutine pf_wait_p

     subroutine pf_broadcast_p(pf, y, nvar, root)
       import pf_pfasst_t, pfdp
       type(pf_pfasst_t), intent(inout) :: pf
       integer,           intent(in)    :: nvar, root
       real(pfdp)  ,      intent(in)    :: y(nvar)
     end subroutine pf_broadcast_p

  end interface

end module pf_mod_dtype
