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

  integer, parameter :: pfdp = c_double

  real(pfdp), parameter :: ZERO  = 0.0_pfdp
  real(pfdp), parameter :: ONE   = 1.0_pfdp
  real(pfdp), parameter :: TWO   = 2.0_pfdp
  real(pfdp), parameter :: HALF  = 0.5_pfdp

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
     procedure(pf_hook_p), pointer, nopass :: proc
  end type pf_hook_t


  ! sweeper type
  type :: pf_sweeper_t
     type(c_ptr) :: ctx
     integer     :: npieces
     procedure(pf_sweep_p),      pointer, nopass :: sweep
     procedure(pf_initialize_p), pointer, nopass :: initialize
     procedure(pf_evaluate_p),   pointer, nopass :: evaluate
     procedure(pf_integrate_p),  pointer, nopass :: integrate
  end type pf_sweeper_t


  ! encap type
  type :: pf_encap_t
     type(c_ptr) :: ctx
     procedure(pf_encap_create_p),  pointer, nopass :: create
     procedure(pf_encap_destroy_p), pointer, nopass :: destroy
     procedure(pf_encap_setval_p),  pointer, nopass :: setval
     procedure(pf_encap_copy_p),    pointer, nopass :: copy
     procedure(pf_encap_pack_p),    pointer, nopass :: pack
     procedure(pf_encap_unpack_p),  pointer, nopass :: unpack
     procedure(pf_encap_axpy_p),    pointer, nopass :: axpy
  end type pf_encap_t


  ! level type
  type :: pf_level_t
     integer     :: nvars = -1          ! number of variables (dofs)
     integer     :: nnodes = -1         ! number of sdc nodes
     integer     :: nsweeps = 1         ! number of sdc sweeps to perform
     integer     :: level = -1          ! level number (1=finest)
     logical     :: Finterp = .false.   ! Interpolate functions instead of solution

     type(pf_encap_t) :: encap
     type(pf_sweeper_t) :: sweeper
     procedure(pf_transfer_p), pointer, nopass :: interpolate, restrict

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

     integer(c_int), pointer :: &
          nflags(:)                     ! sdc node flags

     type(c_ptr), pointer :: &
          qSDC(:), &                    ! unknowns at sdc nodes
          pSDC(:), &                    ! unknowns at sdc nodes, previous sweep
          fSDC(:,:), &                  ! functions values at sdc nodes
          pfSDC(:,:), &                 ! functions at sdc nodes, previous sweep
          tau(:), &                     ! fas correction
          qend                          ! solution at last node

     type(c_ptr) :: ctx  = c_null_ptr   ! user context
     integer, pointer :: shape(:)       ! user shape

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

     ! pf objects
     type(pf_state_t)          :: state
     type(pf_level_t), pointer :: levels(:)
     type(pf_comm_t), pointer  :: comm

     ! hooks
     type(pf_hook_t), pointer :: hooks(:,:)
     integer, pointer         :: nhooks(:)

     ! timing
     logical    :: echo_timings  = .false.
     integer(8) :: timers(100)   = 0
     integer(8) :: runtimes(100) = 0
  end type pf_pfasst_t


  ! hook interface
  interface
     subroutine pf_hook_p(pf, level, state, ctx)
       use iso_c_binding
       import pf_pfasst_t, pf_level_t, pf_state_t
       type(pf_pfasst_t), intent(inout) :: pf
       type(pf_level_t),  intent(inout) :: level
       type(pf_state_t),  intent(in)    :: state
       type(c_ptr),       intent(in)    :: ctx
     end subroutine pf_hook_p
  end interface


  ! sweeper interfaces
  interface
     subroutine pf_sweep_p(pf, F, t0, dt)
       import pf_pfasst_t, pf_level_t, pfdp, c_ptr
       type(pf_pfasst_t), intent(inout) :: pf
       real(pfdp),        intent(in)    :: dt, t0
       type(pf_level_t),  intent(inout) :: F
     end subroutine pf_sweep_p
  end interface

  interface
     subroutine pf_evaluate_p(F, t, m)
       import pf_level_t, pfdp, c_ptr
       type(pf_level_t), intent(inout) :: F
       real(pfdp),       intent(in)    :: t
       integer,          intent(in)    :: m
     end subroutine pf_evaluate_p
  end interface

  interface
     subroutine pf_initialize_p(F)
       import pf_level_t, c_ptr
       type(pf_level_t), intent(inout) :: F
     end subroutine pf_initialize_p
  end interface

  interface
     subroutine pf_integrate_p(F, fSDC, dt, fintSDC)
       import pf_level_t, c_ptr, pfdp
       type(pf_level_t),  intent(in) :: F
       type(c_ptr),       intent(in) :: fSDC(:, :), fintSDC(:)
       real(pfdp),        intent(in) :: dt
     end subroutine pf_integrate_p
  end interface


  ! transfer interfaces
  interface
     subroutine pf_transfer_p(qF, qG, levelF, ctxF, levelG, ctxG)
       import c_ptr
       type(c_ptr), intent(in), value :: qF, qG, ctxF, ctxG
       integer,     intent(in)        :: levelF, levelG
     end subroutine pf_transfer_p
  end interface


  ! encapsulation interfaces
  interface
     subroutine pf_encap_create_p(sol, level, feval, nvars, shape, ctx)
       import c_ptr
       type(c_ptr),  intent(out)       :: sol
       type(c_ptr),  intent(in), value :: ctx
       integer,      intent(in)        :: level, nvars, shape(:)
       logical,      intent(in)        :: feval
     end subroutine pf_encap_create_p
  end interface

  interface
     subroutine pf_encap_destroy_p(sol)
       import c_ptr
       type(c_ptr), intent(in), value :: sol
     end subroutine pf_encap_destroy_p
  end interface

  interface
     subroutine pf_encap_setval_p(sol, val)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value :: sol
       real(pfdp),  intent(in)        :: val
     end subroutine pf_encap_setval_p
  end interface

  interface
     subroutine pf_encap_copy_p(dst, src)
       import c_ptr
       type(c_ptr), intent(in), value :: dst, src
     end subroutine pf_encap_copy_p
  end interface

  interface
     subroutine pf_encap_pack_p(z, q)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value :: q
       real(pfdp),  intent(out)       :: z(:)
     end subroutine pf_encap_pack_p
  end interface

  interface
     subroutine pf_encap_unpack_p(q, z)
       import c_ptr, pfdp
       type(c_ptr), intent(in), value :: q
       real(pfdp),  intent(in)        :: z(:)
     end subroutine pf_encap_unpack_p
  end interface

  interface
     subroutine pf_encap_axpy_p(y, a, x)
       import c_ptr, pfdp
       real(pfdp),  intent(in)        :: a
       type(c_ptr), intent(in), value :: x, y
     end subroutine pf_encap_axpy_p
  end interface


  ! communicator interfaces
  ! XXX


end module pf_mod_dtype
