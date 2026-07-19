!!  IMEX Sweeper Module
!
! This file is part of LIBPFASST.
!
!>  Module for the IMEX Sweeper  of the  the derived sweeper class for doing IMEX sweeps for an equation of the form
!!         $$   y' = f_1(y) + f_2(y)  $$
!!  The \(f_1\) piece is treated explicitly and \(f_2\) implicitly
!!  Afer this sweeper is initialized (usually in main), the logical flags can be changed if desired
!! ---
!!     explicit:  Make false if there is no explicit piece
!!
!!     implicit:  Make false if there is no implicit piece
!! ---
!!  The user needs to supply the feval and fcomp routines for a given example
module pf_mod_imex_diag_sweeper
   use pf_mod_dtype
   use pf_mod_utils
   use pf_mod_timer
   use pf_mod_mpi
   use pf_mod_comm_mpi

   implicit none

   !>  IMEX SDC sweeper type, extends abstract sweeper
   type, extends(pf_sweeper_t), abstract :: pf_imex_diag_sweeper_t
      real(pfdp), allocatable :: QtilE(:,:)   !!  Approximate explicit quadrature rule
      real(pfdp), allocatable :: QtilI(:,:)   !!  Approximate implicit quadrature rule
      real(pfdp), allocatable :: dtsdc(:)     !!  SDC step sizes (single value here)
      real(pfdp), allocatable :: QdiffE(:,:)  !!  qmat-QtilE
      real(pfdp), allocatable :: QdiffI(:,:)  !!  qmat-QtilI

      logical    :: explicit  !!  True if there is an explicit piece (must set in derived sweeper)
      logical    :: implicit  !!  True if there an implicit piece (must set in derived sweeper)
      integer    :: m_sub     !!  Substep loop variable (useful in the function evaluation routines in derived sweepers)
      class(pf_encap_t), allocatable :: rhs        !! holds rhs for implicit solve
      class(pf_encap_t), allocatable :: local_sum  !! holds intermediate results on each rank

      real(pfdp), allocatable :: send_flat(:)  !!  flat array for sending encap data
      real(pfdp), allocatable :: recv_flat(:)  !!  flat array for receiving encap data

   contains
      procedure(pf_f_eval_p), deferred :: f_eval   !!  RHS function evaluations
      procedure(pf_f_comp_p), deferred :: f_comp   !!  Implicit solver
      !>  Set the generic functions
      procedure :: sweep      => imex_diag_sweep
      procedure :: initialize => imex_diag_initialize
      procedure :: evaluate   => imex_diag_evaluate
      procedure :: integrate  => imex_diag_integrate
      procedure :: residual   => imex_diag_residual
      procedure :: spreadq0   => imex_diag_spreadq0
      procedure :: compute_dt => imex_compute_dt
      procedure :: evaluate_all => imex_diag_evaluate_all
      procedure :: destroy   => imex_diag_destroy
      procedure :: imex_diag_destroy
      procedure :: imex_diag_initialize
   end type pf_imex_diag_sweeper_t

   interface
      !>  The interface to the routine to compute the RHS function values
      !>  Evaluate f_piece(y), where piece is one or two
      subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
         import pf_imex_diag_sweeper_t, pf_encap_t, pfdp
         class(pf_imex_diag_sweeper_t),  intent(inout) :: this
         class(pf_encap_t), intent(in   )  :: y           !!  Argument for evaluation
         real(pfdp),        intent(in   )  :: t           !!  Time at evaluation
         integer,    intent(in   )         :: level_index !!  Level index
         class(pf_encap_t), intent(inout)  :: f           !!  RHS function value
         integer,    intent(in   )         :: piece       !!  Which piece to evaluate
      end subroutine pf_f_eval_p

      !>  The interface to the routine to do implicit solve
      !>  i.e, solve the equation y - dtq*f_2(y) =rhs
      subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
         import pf_imex_diag_sweeper_t, pf_encap_t, pfdp
         class(pf_imex_diag_sweeper_t),  intent(inout) :: this
         class(pf_encap_t), intent(inout)  :: y           !!  Solution of implicit solve
         real(pfdp),        intent(in   )  :: t           !!  Time of solve
         real(pfdp),        intent(in   )  :: dtq         !!  dt*quadrature weight
         class(pf_encap_t), intent(in   )  :: rhs         !!  RHS for solve
         integer,    intent(in   )         :: level_index !!  Level index
         class(pf_encap_t), intent(inout) :: f            !!  f_2 of solution y
         integer,    intent(in   ) :: piece               !!  Which piece to evaluate
      end subroutine pf_f_comp_p
   end interface

contains

   !> Perform nsweeps SDC sweeps on level level_index and set qend appropriately.
   subroutine imex_diag_sweep(this, pf, level_index, t0, dt,nsweeps, flags)
      use pf_mod_hooks

      !>  Inputs
      class(pf_imex_diag_sweeper_t), intent(inout) :: this
      type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
      integer,           intent(in)    :: level_index  !!  level on which to sweep
      real(pfdp),        intent(in   ) :: t0           !!  time at beginning of time step
      real(pfdp),        intent(in   ) :: dt           !!  time step size
      integer,           intent(in)    :: nsweeps      !!  number of sweeps to do
      integer, optional, intent(in   ) :: flags        !!  sweep specific flags

      !>  Local variables
      type(pf_level_t), pointer :: lev    !!  points to current level

      integer     :: m,n,k      !!  Loop variables
      real(pfdp)  :: t          !!  Time at nodes

      !> Diagonal sweeper specific
      integer :: rank, fidx, ierror

      !>  Assign level pointer
      lev => pf%levels(level_index)   ! Assign level pointer
      rank = pf%rank_diag             ! Diag rank is index for nodes in diagonal sweeper (shifted by 2 or 1 since first node is not stored explicitly in PFASST)
      fidx = 1 + merge(1,0,pf%rank_diag .eq. 0)        ! We need this since rank 0 holds two F-arrays for node 1 and node 2, while all other ranks hold only one F-array for their node

      call this%local_sum%setval(0.0_pfdp)

      sweeps: do k = 1,nsweeps   !!  Loop over sweeps
         pf%state%sweep=k
         call call_hooks(pf, level_index, PF_PRE_SWEEP)
         call pf_start_timer(pf, T_SWEEP,level_index)

         !  Add terms from previous iteration  (not passing CI tests)
         !do m = 1, lev%nnodes-1
         !   call lev%I(m)%setval(0.0_pfdp)
         !end do

         ! compute integrals and add fas correction
         call lev%I(1)%setval(0.0_pfdp)
         ! initial value from left-node is omitted since QdiffE is zero in first column
         do m = 1, lev%nnodes-1    ! node m here is stored on rank m-1 since first node is not stored in PFASST
            ! add explicit part
            if (this%explicit) then
               ! loop over nodes replaced by MPI_Reduce with sum
               ! omit node 1 since QdiffE is zero in first column
               call this%local_sum%setval(0.0_pfdp)
               call this%local_sum%axpy(dt * this%QdiffE(m,rank+2),lev%F(fidx,1))
               call this%local_sum%pack(this%send_flat)       ! encap to array
               ! collect at rank m-1 since node m is stored on rank m-1
               if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce sweep explicit before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'nproc=', pf%diag_comm%nproc
               call mpi_reduce(this%send_flat,this%recv_flat,lev%mpibuflen,myMPI_Datatype,MPI_SUM,m-1,pf%diag_comm%comm,ierror)
               if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce sweep explicit after :', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'ierror=', ierror
               if (pf%rank_diag == m-1) then
                  call this%local_sum%unpack(this%recv_flat)     ! array to encap
                  call lev%I(1)%axpy(1.0_pfdp, this%local_sum)   ! add to integral
               end if
            end if
            ! add implicit part
            if (this%implicit) then
               ! loop over nodes replaced by MPI_Reduce with sum
               ! omit node 1 since QdiffI is zero in first column
               call this%local_sum%setval(0.0_pfdp)
               call this%local_sum%axpy(dt * this%QdiffI(m,rank+2),lev%F(fidx,2))
               call this%local_sum%pack(this%send_flat)       ! encap to array
               ! collect at rank m-1 since node m is stored on rank m-1
               if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce sweep implicit before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'nproc=', pf%diag_comm%nproc
               call mpi_reduce(this%send_flat,this%recv_flat,lev%mpibuflen,myMPI_Datatype,MPI_SUM,m-1,pf%diag_comm%comm,ierror)
               if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce sweep implicit after :', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'ierror=', ierror
               if (pf%rank_diag == m-1) then
                  call this%local_sum%unpack(this%recv_flat)     ! array to encap
                  call lev%I(1)%axpy(1.0_pfdp, this%local_sum)   ! add to integral
               end if
            end if
         end do   ! node loop

         !  Add the tau FAS correction
         ! Removed correction step for use_Sform since it is not compatible with diagonal sweeper
         if (level_index < pf%state%finest_level) then
            ! no outer loop required since each rank holds single I & tauQ
            call lev%I(1)%axpy(1.0_pfdp, lev%tauQ(1))
         end if

         !  Recompute the first function value if this is first sweep
         ! We only need to call this on rank 0
         ! Rank 0 holds 2 F-arrays F(1,:) for node 1 and F(2,:) for node 2
         if (k .eq. 1 .and. rank == 0) then
            ! explicit part
            if (this%explicit) then
               if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'f_eval sweep explicit before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'nproc=', pf%diag_comm%nproc
               call pf_start_timer(pf,T_FEVAL,level_index)
               call this%f_eval(lev%q0, t0, level_index, lev%F(1,1),1)
               call pf_stop_timer(pf,T_FEVAL,level_index)
               ! broadcast to all ranks ?
               !call lev%F(1,1)%pack(this%send_flat)    ! encap to array
               !call mpi_bcast(this%send_flat,lev%mpibuflen,myMPI_Datatype,0,pf%diag_comm, ierror)
               !call lev%F(1,1)%unpack(this%send_flat)  ! array to encap
            end if
            ! implicit part
            if (this%implicit) then
               if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'feval sweep implicit before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'nproc=', pf%diag_comm%nproc
               call pf_start_timer(pf,T_FEVAL,level_index)
               call this%f_eval(lev%q0, t0, level_index, lev%F(1,2),2)
               call pf_stop_timer(pf,T_FEVAL,level_index)
               ! broadcast to all ranks ?
               !call lev%F(1,2)%pack(this%send_flat)    ! encap to array
               !call mpi_bcast(this%send_flat,lev%mpibuflen,myMPI_Datatype,0,pf%diag_comm, ierror)
               !call lev%F(1,2)%unpack(this%send_flat)  ! array to encap
            end if
         end if

         !t = t0

         !> loop over substeps is parallel loop over ranks here
         ! boils down to a single implicit step due to diagonal QtilI and zero QtilE
         t = t0 + dt*this%dtsdc(1)

         this%m_sub=rank+1
         !>  Accumulate rhs
         call this%rhs%setval(0.0_pfdp)
         ! only need to add the implicit part since QtilE is zero
         ! omitting the first node since QtilI is zero in first column
         !if (this%implicit) then
         ! isn't this also not 0, since org. sweeper goes n=1,m ?? which is just lower-triangular matrix
         !   call this%rhs%axpy(dt*this%QtilI(rank+1,rank+2), lev%F(fidx,2))
         !end if
         !>  Add the integral term
         call this%rhs%axpy(1.0_pfdp, lev%I(1))

         !>  Add the starting value - changed from lev%Q(1) to lev%q0
         call this%rhs%axpy(1.0_pfdp, lev%q0)

         !>  Solve for the implicit piece
         if (this%implicit) then
            if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'fcomp sweep implicit before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'nproc=', pf%diag_comm%nproc
            call pf_start_timer(pf,T_FCOMP,level_index)
            call this%f_comp(lev%Q(1), t, dt*this%QtilI(rank+1,rank+2), this%rhs, level_index,lev%F(fidx,2),2)
            call pf_stop_timer(pf,T_FCOMP,level_index)
         else
            call lev%Q(1)%copy(this%rhs)
         end if

         !>  Compute explicit function on new value
         if (this%explicit) then
            call pf_start_timer(pf,T_FEVAL,level_index)
            call this%f_eval(lev%Q(1), t, level_index, lev%F(fidx,1),1)
            call pf_stop_timer(pf,T_FEVAL,level_index)
         end if

         ! all substeps computed in parallel

         ! compute residual at every node/rank
         call pf_residual(pf, level_index, dt)
         !> Identify the rank holding the last node, broadcast and copy into qend.
         !> MPI_Bcast must be called by all ranks in the communicator, so the
         !> non-root ranks need to enter the collective as well.
         if (rank+2 .eq. lev%nnodes) then
            !> Copy the last node to qend on the root rank.
            call lev%qend%pack(this%send_flat)    ! encap to array
         end if
         if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'bcast sweep end before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', lev%nnodes-2, 'nproc=', pf%diag_comm%nproc
         call mpi_bcast(this%send_flat,lev%mpibuflen,myMPI_Datatype,lev%nnodes-2,pf%diag_comm%comm, ierror)
         call lev%qend%unpack(this%send_flat)  ! array to encap
         call pf_stop_timer(pf, T_SWEEP,level_index)

         call call_hooks(pf, level_index, PF_POST_SWEEP)
      end do sweeps  !  End loop on sweeps

   end subroutine imex_diag_sweep

   !> Subroutine to initialize matrices and space for sweeper
   subroutine imex_diag_initialize(this, pf, level_index)
      use pf_mod_quadrature
      class(pf_imex_diag_sweeper_t), intent(inout) :: this
      type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
      integer,           intent(in)    :: level_index  !!  level on which to initialize

      type(pf_level_t), pointer  :: lev    !  Current level

      integer    ::  nnodes,ierr,i,j
      lev => pf%levels(level_index)   !  Assign level pointer
      this%npieces = 2

      !  The default is to use both pieces, but can be overriddent in local sweeper
      this%explicit=.TRUE.
      this%implicit=.TRUE.

      !
      nnodes = lev%nnodes
      if (pf%diag_comm%nproc /= nnodes - 1) then
         ! diagonal sweeper only doable if diagsize
         call pf_stop(__FILE__,__LINE__,"initialization error - number of diagonal processes need to be nnodes-1")
      end if

      !
      allocate(this%QdiffE(nnodes-1,nnodes),stat=ierr)
      if (ierr /=0) stop "allocate fail in imex_diag_initialize for QdiffE"
      allocate(this%QdiffI(nnodes-1,nnodes),stat=ierr)
      if (ierr /=0) stop "allocate fail in imex_diag_initialize for QdiffI"
      allocate(this%QtilE(nnodes-1,nnodes),stat=ierr)
      if (ierr /=0) stop "allocate fail in imex_diag_initialize for QtilE"
      allocate(this%QtilI(nnodes-1,nnodes),stat=ierr)
      if (ierr /=0) stop "allocate fail in imex_diag_initialize for QtilI"
      allocate(this%dtsdc(1),stat=ierr)      ! only single value here due to diagonal sweeper
      if (ierr /=0) stop "allocate fail in imex_diag_initialize for dtsdc"
      ! allocate additional space for flat arrays for sending and receiving encap data
      allocate(this%send_flat(lev%mpibuflen),stat=ierr)
      if (ierr /=0) stop "allocate fail in imex_diag_initialize for send_flat"
      allocate(this%recv_flat(lev%mpibuflen),stat=ierr)
      if (ierr /=0) stop "allocate fail in imex_diag_initialize for recv_flat"

      this%QtilE = 0.0_pfdp
      this%QtilI = 0.0_pfdp

      !>  Substep sizes (zero-to-node) for diagonal sweeper
      ! Use
      this%dtsdc(1) = lev%sdcmats%qnodes(pf%rank_diag+2) - lev%sdcmats%qnodes(1)

      ! Implicit matrix
      ! Use row sum to compute diagonal version of QtilI
      do i = 1, nnodes - 1
         do j = 1, i+1
            if (this%use_LUq) then
               this%QtilI(i,i+1) = this%QtilI(i,i+1) + lev%sdcmats%qmatLU(i,j)
            else
               this%QtilI(i,i+1) = this%QtilI(i,i+1) +lev%sdcmats%qmatBE(i,j)
            end if
         end do
      end do

      ! Explicit matrix
      ! Explicit case diagonal Qtil is Zero matrix
      ! this%QtilE =  lev%sdcmats%qmatFE

      this%QdiffE = lev%sdcmats%qmat-this%QtilE
      this%QdiffI = lev%sdcmats%qmat-this%QtilI

      if (pf%use_Sform) then
         ! node-to-node not possible with diagonal sweeper
         call pf_stop(__FILE__,__LINE__,"initialization error - pf%useSform=true for diagonal sweeper - how did we came here?")
      end if
      !>  Make space for rhs & local_sum
      call lev%ulevel%factory%create_single(this%rhs, level_index,   lev%lev_shape)
      call lev%ulevel%factory%create_single(this%local_sum, level_index,   lev%lev_shape)


   end subroutine imex_diag_initialize

   !>  Subroutine to deallocate sweeper
   subroutine imex_diag_destroy(this, pf,level_index)
      class(pf_imex_diag_sweeper_t),  intent(inout) :: this
      type(pf_pfasst_t),  target,  intent(inout) :: pf
      integer,              intent(in)    :: level_index

      type(pf_level_t), pointer  :: lev        !  Current level
      lev => pf%levels(level_index)   !  Assign level pointer

      deallocate(this%QdiffE)
      deallocate(this%QdiffI)
      deallocate(this%QtilE)
      deallocate(this%QtilI)
      deallocate(this%dtsdc)
      ! destroy additional space for flat arrays for sending and receiving encap data
      deallocate(this%send_flat)
      deallocate(this%recv_flat)

      call lev%ulevel%factory%destroy_single(this%rhs)
      call lev%ulevel%factory%destroy_single(this%local_sum)

   end subroutine imex_diag_destroy


   !> Subroutine to compute  Picard integral of function values
   subroutine imex_diag_integrate(this,pf,level_index, qSDC, fSDC, dt, fintSDC, flags)
      class(pf_imex_diag_sweeper_t), intent(inout) :: this
      type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
      integer,           intent(in)    :: level_index  !!  level on which to initialize
      class(pf_encap_t), intent(in   ) :: qSDC(:)      !!  Solution values
      class(pf_encap_t), intent(in   ) :: fSDC(:, :)   !!  RHS Function values
      real(pfdp),        intent(in   ) :: dt           !!  Time step
      class(pf_encap_t), intent(inout) :: fintSDC(:)   !!  Integral from t_1 to t_m
      integer, optional, intent(in   ) :: flags

      integer :: m, rank, fidx, ierror
      type(pf_level_t), pointer :: lev        !  Current level
      lev => pf%levels(level_index)   !  Assign level pointer
      rank = pf%rank_diag
      fidx = 1 + merge(1,0,rank .eq. 0)        ! We need this since rank 0 holds two F-arrays for node 1 and node 2, while all other ranks hold only one F-array for their node

      !
      call fintSDC(1)%setval(0.0_pfdp)
      !> again ignoring first node since QdiffE and QdiffI are zero in first column
      do m = 1, lev%nnodes-1
         ! explicit part
         if (this%explicit) then
            call this%local_sum%setval(0.0_pfdp)
            call this%local_sum%axpy(dt*lev%sdcmats%qmat(m,rank+2), fSDC(fidx,1))
            call this%local_sum%pack(this%send_flat)       ! encap to array
            if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce explicit before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'nproc=', pf%diag_comm%nproc
            call mpi_reduce(this%send_flat,this%recv_flat,lev%mpibuflen,myMPI_Datatype,MPI_SUM,m-1,pf%diag_comm%comm,ierror)
            if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce explicit after :', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'ierror=', ierror
            if (pf%rank_diag == m-1) then
               call this%local_sum%unpack(this%recv_flat)     ! array to encap
               call fintSDC(1)%axpy(1.0_pfdp, this%local_sum) ! add to integral
            end if
         end if
         ! implicit part
         if (this%implicit) then
            call this%local_sum%setval(0.0_pfdp)
            call this%local_sum%axpy(dt*lev%sdcmats%qmat(m,rank+2), fSDC(fidx,2))
            call this%local_sum%pack(this%send_flat)       ! encap to array
            if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce implicit before:', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'nproc=', pf%diag_comm%nproc
            call mpi_reduce(this%send_flat,this%recv_flat,lev%mpibuflen,myMPI_Datatype,MPI_SUM,m-1,pf%diag_comm%comm,ierror)
            if (pf%debug) print*, 'DEBUG --', pf%rank_global, 'diag_reduce implicit after :', 'rank=', pf%rank, 'rank_diag=', pf%rank_diag, 'm=', m, 'root=', m-1, 'ierror=', ierror
            if (pf%rank_diag == m-1) then
               call this%local_sum%unpack(this%recv_flat)     ! array to encap
               call fintSDC(1)%axpy(1.0_pfdp, this%local_sum) ! add to integral
            end if
         end if

      end do

      !if (this%explicit) call pf_apply_mat(fintSDC, dt, lev%sdcmats%Qmat, fSDC(:,1), .false.)
      !if (this%implicit) call pf_apply_mat(fintSDC, dt, lev%sdcmats%Qmat, fSDC(:,2), .false.)
   end subroutine imex_diag_integrate


   !> Subroutine to compute  Residual
   subroutine imex_diag_residual(this, pf, level_index, dt, flags)
      class(pf_imex_diag_sweeper_t),  intent(inout) :: this
      type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
      integer,           intent(in)    :: level_index  !!  level on which to initialize
      real(pfdp),        intent(in   ) :: dt           !!  Time step
      integer, intent(in), optional   :: flags

      type(pf_level_t), pointer :: lev        !  Current level
      lev => pf%levels(level_index)   !  Assign level pointer

      call imex_diag_integrate(this,pf,level_index, lev%Q, lev%F, dt, lev%I, flags)
      !> add tau if it exists
      if (lev%index < pf%state%finest_level) then
         !> again ignoring first node
         ! single addition per rank
         call lev%I(1)%axpy(1.0_pfdp, lev%tauQ(1), flags)
      end if
      !> compute residuals at each node
      call lev%R(1)%copy(lev%I(1))
      call lev%R(1)%axpy(-1.0_pfdp, lev%Q(1))
      !> in general no flags should be used
      if (present(flags)) then
         if (flags .eq. 0) then
            call lev%R(1)%axpy(1.0_pfdp, lev%q0)
         else
            call lev%R(1)%axpy(1.0_pfdp, lev%Q(1))
         end if
      else
         !> changes default case to use q0 instead of Q(1) for residual computation
         call lev%R(1)%axpy(1.0_pfdp, lev%q0)
      end if

      !call pf_generic_residual(this, pf, level_index, dt)
   end subroutine imex_diag_residual


   subroutine imex_diag_spreadq0(this, pf,level_index, t0, flags, step)
      class(pf_imex_diag_sweeper_t),  intent(inout) :: this
      type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
      integer,           intent(in)    :: level_index  !!  level on which to initialize
      real(pfdp),        intent(in   ) :: t0
      integer, optional,   intent(in)    :: flags, step

      integer :: p
      class(pf_level_t), pointer :: lev  !!  Level on which to spread
      lev => pf%levels(level_index)   !!  Assign level pointer

      !  Stick initial condition into node slot - happens on all ranks
      call lev%Q(1)%copy(lev%q0)

      !  Evaluate F at all ranks for t0 - happens on all ranks
      call lev%ulevel%sweeper%evaluate(pf,level_index, t0, 1)

      ! No need to spread F and Q
      ! Copy F on rank 0 into second spot (holds F for node 1 and 2)
      if (pf%rank_diag .eq. 0) then
         do p = 1, this%npieces
            call lev%F(2,p)%copy(lev%F(1,p))
         end do
      end if
   end subroutine imex_diag_spreadq0

   subroutine imex_compute_dt(this,pf,level_index,  t0, dt,flags)
      class(pf_imex_diag_sweeper_t),  intent(inout) :: this
      type(pf_pfasst_t), target, intent(inout) :: pf
      integer,              intent(in)    :: level_index
      real(pfdp),        intent(in   ) :: t0
      real(pfdp),        intent(inout) :: dt
      integer, optional,   intent(in)    :: flags

      type(pf_level_t),    pointer :: lev
      lev => pf%levels(level_index)   !!  Assign level pointer
      !  Do nothing now
      return
   end subroutine imex_compute_dt

   !> Subroutine to evaluate function value at node m
   ! only used in imex_diag_spreadq0()
   subroutine imex_diag_evaluate(this, pf,level_index, t, m, flags, step)
      class(pf_imex_diag_sweeper_t),  intent(inout) :: this
      type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
      integer,           intent(in)    :: level_index  !!  level on which to initialize
      real(pfdp),        intent(in   ) :: t    !!  Time at which to evaluate
      integer,           intent(in   ) :: m    !!  Node at which to evaluate
      integer, intent(in), optional   :: flags, step

      integer :: fidx
      type(pf_level_t), pointer :: lev        !  Current level
      lev => pf%levels(level_index)   !  Assign level pointer
      fidx = 1 + merge(1,0,pf%rank_diag .eq. 0)    !!  We need this since rank 0 holds two F-arrays for node 1 and node 2, while all other ranks hold only one F-array for their node

      ! no need for m here since each rank holds only one node
      ! we keep it for consistency with generic sweeper
      ! do some sanity check here - should only be called for m=1
      if (m > 1) then
         call pf_stop(__FILE__,__LINE__,"imex_diag_evaluate: m>1 - evaluate only for intial condition")
      end if
      if (this%explicit) then
         call pf_start_timer(pf,T_FEVAL,level_index)
         call this%f_eval(lev%Q(1), t, level_index, lev%F(fidx,1),1)
         call pf_stop_timer(pf,T_FEVAL,level_index)
      end if

      if (this%implicit) then
         call pf_start_timer(pf,T_FEVAL,level_index)
         call this%f_eval(lev%Q(1), t, level_index, lev%F(fidx,2),2)
         call pf_stop_timer(pf,T_FEVAL,level_index)
      end if

   end subroutine imex_diag_evaluate

   !> Subroutine to evaluate the function values at all nodes
   subroutine imex_diag_evaluate_all(this, pf,level_index, t, flags, step)
      class(pf_imex_diag_sweeper_t),  intent(inout) :: this
      type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
      integer,           intent(in)    :: level_index  !!  level on which to initialize
      real(pfdp),        intent(in   ) :: t(:)  !!  Array of times at each node
      integer, intent(in), optional   :: flags, step

      integer :: fidx
      class(pf_level_t), pointer :: lev   !!  points to current level
      lev => pf%levels(level_index)       !!  Assign level pointer
      fidx = 1 + merge(1,0,pf%rank_diag .eq. 0)    !!  We need this since rank 0 holds two F-arrays for node 1 and node 2, while all other ranks hold only one F-array for their node

      if (this%explicit) then
         call pf_start_timer(pf,T_FEVAL,level_index)
         call this%f_eval(lev%Q(1), t(pf%rank_diag+2), level_index, lev%F(fidx,1),1)
         call pf_stop_timer(pf,T_FEVAL,level_index)
      end if

      if (this%implicit) then
         call pf_start_timer(pf,T_FEVAL,level_index)
         call this%f_eval(lev%Q(1), t(pf%rank_diag+2), level_index, lev%F(fidx,2),2)
         call pf_stop_timer(pf,T_FEVAL,level_index)
      end if

   end subroutine imex_diag_evaluate_all

end module pf_mod_imex_diag_sweeper
