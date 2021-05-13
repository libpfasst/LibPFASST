!!  Routines that run the MGRIT algorithm
!
! This file is part of LIBPFASST.
!

!> Module of routines to run MGRIT
module pf_mod_MGRIT
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_pfasst
  use pf_mod_comm
  use pf_mod_results
  use pf_mod_parareal
  use mpi
  
  implicit none

  type :: int_vector
     integer, allocatable :: val(:)
  end type int_vector

  !> This data structure holds important information for each MGRIT level for this proc
  type :: mgrit_level_data
     integer :: Nt !> Number of time points local to this proc
     integer :: Nt_glob !> global number of time points
     real(pfdp) :: dt !> Time step size
     real(pfdp) :: T0_glob !> Global starting time step
     real(pfdp) :: Tfin_glob !> Global ending time step
     real(pfdp) :: t0 !> Starting time step local to this proc
     real(pfdp) :: tfin !> Ending time step local to this proc
     real(pfdp) :: res_norm_loc !> Residual norm local to this proc
     integer :: res_norm_index 
     integer :: cycle_phase !> current phase of the V-cycle (down-cycle, exact solve, up-cycle)
     logical :: FCF_flag = .true. !> If false, use only F-relaxation (set to false for Parareal)
     logical :: FAS_flag = .false. !> If .true., use the Full Approximation Scheme (FAS)
     logical :: setup_start_coarse_flag = .false. !> If true, set up the hierarchy of time grids starting from the coarsest
     type(int_vector), allocatable :: f_blocks(:) !> Blocks of F-points. Each block contains the indices of the F-points
     integer, allocatable :: c_pts(:) !> Array of C-point indices
     class(pf_encap_t), allocatable :: g(:) !> Right-hand side
     class(pf_encap_t), allocatable :: uc(:) !> Values at C-points
     class(pf_encap_t), allocatable :: uc_prev(:) !> Temp variable for values at C-points
     class(pf_encap_t), allocatable :: uc_fas(:) !> FAS approximation at C-points
     class(pf_encap_t), allocatable :: r !> Residual
     class(pf_encap_t), allocatable :: uc_ghost !> Ghost points (data received by neighboring procs)
     class(pf_encap_t), allocatable :: uc_fas_ghost !> FAS ghost points
     class(pf_encap_t), allocatable :: u_temp !> Temp encap
     class(pf_encap_t), allocatable :: Q0 !> Global initial condition
     class(pf_encap_t), allocatable :: g_temp !> Temp encap
     class(pf_encap_t), allocatable :: e(:) !> Correction to be interpolated 
     integer :: send_to_rank !> Send neighbor
     integer :: recv_from_rank !> Receive neighbor
     integer :: rank_shifted !> Shifted rank.  Needed for coarser grids when the global number of time points < num procs
     integer :: coarsest_level !> Index of coarsest level
     logical :: f_pts_flag = .true. !> True if this proc has F-points
     logical :: c_pts_flag = .true. !> True if this proc has C-points
     logical :: c_converge
     integer :: c_converge_index
     integer, allocatable :: recv_procs(:)
     integer, allocatable :: send_procs(:)
     integer, allocatable :: recv_converge_status(:)
     integer :: num_recvs
     integer :: num_sends
     integer :: num_recv_converged
  end type 
  
contains

  !> Set up MGRIT data
  subroutine mgrit_initialize(pf, mg_ld, T0_glob, Tfin_glob, Nt_start, coarsen_factor, FAS_flag, FCF_flag, start_coarse_flag)
     type(pf_pfasst_t), intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: Nt_start, coarsen_factor
     logical, intent(in) :: FAS_flag, FCF_flag, start_coarse_flag
     real(pfdp) :: T0_glob, Tfin_glob, x, y
     integer :: level_index, nlevels, N, i, kk, p, Nt, p_rank, k
     type(pf_level_t) :: pf_lev, pf_f_lev, pf_c_lev
     type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev
     integer :: count_levels
     logical :: start_dropping_procs
     integer :: my_rank, rank_send, rank_recv, Nt_temp
     integer :: Nt_coarsest, coarsest_level
     integer :: num_sends, num_recvs

     if (start_coarse_flag .eqv. .false.) then !> Generate starting from coarsest grid
        x = real(pf%comm%nproc, pfdp)
        y = log(x) / log(real(coarsen_factor, pfdp))
        if (ceiling(y) .ne. floor(y)) then
           print *,'Error in MGRIT setup: number of procs must be a power of the coarsening factor'
           stop
        end if
        x = real(Nt_start * pf%comm%nproc, pfdp)
        y = log(x) / log(real(coarsen_factor, pfdp))
        if (ceiling(y) .ne. floor(y)) then
           print *,'Error in MGRIT setup: number of finest-grid points must be a power of the coarsening factor'
           stop
        end if
     end if

     Nt_coarsest = coarsen_factor

     nlevels = pf%nlevels
     n = coarsen_factor

     allocate(mg_ld(nlevels))
     mg_ld(nlevels)%setup_start_coarse_flag = start_coarse_flag

     do level_index = 1,nlevels
        mg_lev => mg_ld(level_index)
        mg_lev%T0_glob = T0_glob
        mg_lev%Tfin_glob = Tfin_glob
     end do
     mg_ld(nlevels)%num_sends = 1
     mg_ld(nlevels)%num_recvs = 1
     !> Generate the hierarchy of grids
     if (mg_ld(nlevels)%setup_start_coarse_flag .eqv. .true.) then !> Generate starting from the coarsest grid
        coarsest_level = 1
        mg_ld(nlevels)%coarsest_level = 1
        N = Nt_start
        do level_index = 1,nlevels
           mg_lev => mg_ld(level_index)
           mg_lev%Nt = N
           mg_lev%Nt_glob = N * pf%comm%nproc
           mg_lev%dt = (Tfin_glob - T0_glob)/real(mg_lev%Nt_glob, pfdp)
           mg_lev%t0 = T0_glob + real(N, pfdp) * mg_lev%dt * real(pf%rank, pfdp) + mg_lev%dt
           mg_lev%tfin = mg_lev%t0 + N * mg_lev%dt - mg_lev%dt
           pf%state%t0 = mg_lev%t0
           N = N * coarsen_factor
           
           mg_lev%rank_shifted = pf%rank
           mg_lev%send_to_rank = pf%rank+1
           mg_lev%recv_from_rank = pf%rank-1
        end do
     else !> Generate starting from finest grid
        Nt_temp = Nt_start * pf%comm%nproc
        level_index = nlevels
        mg_ld(level_index)%Nt_glob = Nt_temp
        !> Loop until coarsest grid or maximum number of levels is reached
        do while (level_index .gt. 1)
           Nt_temp = Nt_temp / coarsen_factor
           !> If the grid cannot be coarser, stop
           if (Nt_temp < Nt_coarsest) then
              exit
           end if
           level_index = level_index - 1
           mg_ld(level_index)%Nt_glob = Nt_temp
        end do
        coarsest_level = level_index
        mg_ld(nlevels)%coarsest_level = coarsest_level
        mg_ld(nlevels)%Nt = Nt_start
        mg_ld(nlevels)%rank_shifted = pf%rank
        mg_ld(nlevels)%send_to_rank = pf%rank+1
        mg_ld(nlevels)%recv_from_rank = pf%rank-1
        start_dropping_procs = .false.
        i = 0
        !> Loop over levels to determine number of time points, C-points, and F-points
        do level_index = (nlevels-1),coarsest_level,-1
           mg_lev => mg_ld(level_index)
           mg_f_lev => mg_ld(level_index+1)
           mg_lev%send_to_rank = -1
           mg_lev%recv_from_rank = -1
           mg_lev%rank_shifted = -1
           p = pf%rank+1
           !> Several cases need to be handled if the fine grid has only one point
           if (mg_f_lev%Nt .eq. 1) then
              if (start_dropping_procs .eqv. .false.) then !> If there is one point per proc, we need to start "dropping" procs
                 start_dropping_procs = .true.
              end if
              if (mod(p, coarsen_factor**(i+1)) .eq. 0) then !> If my rank divides coarsen_factor^(i+1), the fine grid has a C-point
                 mg_lev%Nt = 1
                 mg_lev%rank_shifted = (mg_f_lev%rank_shifted+1)/ coarsen_factor - 1
                 mg_lev%send_to_rank = (mg_lev%rank_shifted+2) * coarsen_factor**(i+1) - 1
                 mg_lev%recv_from_rank = mg_lev%rank_shifted * coarsen_factor**(i+1) - 1
                 mg_f_lev%f_pts_flag = .false.
                 mg_f_lev%c_pts_flag = .true.
              else !> Otherwise, the fine grid has an F-point
                 mg_f_lev%f_pts_flag = .true.
                 mg_f_lev%c_pts_flag = .false.
                 mg_lev%c_pts_flag = .false.
                 mg_lev%f_pts_flag = .false.
                 mg_lev%Nt = 0
              end if
           else if (mg_f_lev%Nt .gt. 1) then
              mg_lev%Nt = mg_f_lev%Nt / coarsen_factor
              mg_lev%rank_shifted = pf%rank
              mg_lev%send_to_rank = pf%rank+1
              mg_lev%recv_from_rank = pf%rank-1
           else
              mg_lev%Nt = 0
              mg_lev%c_pts_flag = .false.
              mg_lev%f_pts_flag = .false.
           end if
           if (start_dropping_procs .eqv. .true.) then
              i = i + 1;
              if ((mg_lev%send_to_rank .gt. -1) .and. (mg_lev%send_to_rank .lt. pf%comm%nproc)) then
                 mg_ld(nlevels)%num_sends = mg_ld(nlevels)%num_sends + 1
              end if
              if ((mg_lev%recv_from_rank .gt. -1) .and. (mg_lev%recv_from_rank .lt. pf%comm%nproc)) then
                 mg_ld(nlevels)%num_recvs = mg_ld(nlevels)%num_recvs + 1
              end if
           end if
        end do

        do level_index = coarsest_level,nlevels
           mg_lev => mg_ld(level_index)
           mg_lev%dt = (Tfin_glob - T0_glob)/mg_lev%Nt_glob
           mg_lev%t0 = T0_glob + mg_lev%dt * pf%rank * mg_lev%Nt + mg_lev%dt 
           mg_lev%tfin = mg_lev%t0 + mg_lev%Nt * mg_lev%dt
           pf%state%t0 = mg_lev%t0
        end do

        !do level_index = nlevels,coarsest_level,-1
        !   mg_lev => mg_ld(level_index)
        !   print *,level_index,pf%rank,mg_lev%send_to_rank,mg_lev%recv_from_rank
        !end do
     end if

     allocate(mg_ld(nlevels)%send_procs(mg_ld(nlevels)%num_sends))
     allocate(mg_ld(nlevels)%recv_procs(mg_ld(nlevels)%num_recvs))
     allocate(mg_ld(nlevels)%recv_converge_status(mg_ld(nlevels)%num_recvs))

     mg_ld(nlevels)%recv_procs(1) = -1
     i = 1;
     do level_index = coarsest_level,nlevels
        mg_lev => mg_ld(level_index)
        mg_c_lev => mg_ld(level_index-1)
        if (mg_ld(nlevels)%recv_procs(i) .eq. -1) then
           mg_ld(nlevels)%recv_procs(i) = mg_lev%recv_from_rank
        else
           if ((mg_lev%recv_from_rank .gt. -1) .and. (mg_lev%recv_from_rank .lt. pf%comm%nproc)) then
              if (mg_c_lev%recv_from_rank .ne. mg_lev%recv_from_rank) then
                 i = i + 1
                 mg_ld(nlevels)%recv_procs(i) = mg_lev%recv_from_rank
              end if
           end if
        end if
     end do

     mg_ld(nlevels)%send_procs(1) = mg_ld(nlevels)%send_to_rank
     i = 2
     do level_index = nlevels,coarsest_level+1,-1
        mg_lev => mg_ld(level_index)
        mg_c_lev => mg_ld(level_index-1)
        if ((mg_c_lev%send_to_rank .gt. -1) .and. (mg_c_lev%send_to_rank .lt. pf%comm%nproc)) then
           if (mg_c_lev%send_to_rank .ne. mg_lev%send_to_rank) then
              mg_ld(nlevels)%send_procs(i) = mg_c_lev%send_to_rank
              i = i + 1
           end if
        end if 
     end do

     !print *,pf%rank,mg_ld(nlevels)%num_recvs

     !if (pf%rank .eq. 1) then
     !   do i = 1,mg_ld(nlevels)%num_sends
     !      print *,mg_ld(nlevels)%send_procs(i)
     !   end do
     !end if
     
     pf%state%t0 = mg_ld(nlevels)%tfin - mg_ld(nlevels)%dt
     pf%state%dt = mg_ld(nlevels)%dt
     pf%state%step = (pf%rank+1) * mg_ld(nlevels)%Nt - 1
     pf%state%nsteps = pf%comm%nproc * mg_ld(nlevels)%Nt

     !> Set up the F and C points
     do level_index = nlevels,coarsest_level+1,-1
        mg_lev => mg_ld(level_index)
        mg_c_lev => mg_ld(level_index-1)

        if (mg_ld(nlevels)%setup_start_coarse_flag .eqv. .true.) then
           call FC_Setup(mg_c_lev%Nt, mg_c_lev%Nt, coarsen_factor, mg_lev%f_blocks, mg_lev%c_pts)
        else
           call FC_Setup(mg_c_lev%Nt, mg_lev%Nt, coarsen_factor, mg_lev%f_blocks, mg_lev%c_pts)
        end if
     end do

     !> Allocate memory for all levels
     do level_index = nlevels,coarsest_level,-1
        mg_lev => mg_ld(level_index)
        pf_lev = pf%levels(level_index)
        !if ((level_index .lt. nlevels) .and. (level_index .gt. 1) .and. (mg_lev%Nt .gt. 0)) then
        if (level_index .lt. nlevels) then
           call pf_lev%ulevel%factory%create_array(mg_lev%g, max(1, mg_lev%Nt), level_index, pf_lev%lev_shape)
           do i = 1,mg_lev%Nt
              call mg_lev%g(i)%setval(0.0_pfdp)
           end do
        end if
        if (FAS_flag .eqv. .true.) then
           mg_ld(level_index)%FAS_flag = .true.
           !if ((level_index .lt. nlevels) .and. (mg_lev%Nt .gt. 0)) then
           if (level_index .lt. nlevels) then
              call pf_lev%ulevel%factory%create_array(mg_lev%uc_fas, max(1, mg_lev%Nt), level_index, pf_lev%lev_shape)
              do i = 1,mg_lev%Nt
                 call mg_lev%uc_fas(i)%setval(0.0_pfdp)
              end do
           end if
        end if
        if (FCF_flag .eqv. .false.) then
           mg_ld(level_index)%FCF_flag = .false.
        end if
        if (level_index .gt. coarsest_level) then
           mg_c_lev => mg_ld(level_index-1)
           !if ((mg_c_lev%Nt .eq. 0) .and. (mg_lev%Nt .eq. 1)) then
           !   call pf_lev%ulevel%factory%create_array(mg_lev%uc, 1, level_index, pf_lev%lev_shape)
           !   call pf_lev%ulevel%factory%create_array(mg_lev%uc_prev, 1, level_index, pf_lev%lev_shape)
           !   call mg_lev%uc(1)%setval(0.0_pfdp)
           !   call mg_lev%uc_prev(1)%setval(0.0_pfdp)
           !else if (mg_c_lev%Nt .gt. 0) then
              call pf_lev%ulevel%factory%create_array(mg_lev%uc, max(1, mg_c_lev%Nt), level_index, pf_lev%lev_shape)
              call pf_lev%ulevel%factory%create_array(mg_lev%uc_prev, max(1, mg_c_lev%Nt), level_index, pf_lev%lev_shape)
              do i = 1,mg_c_lev%Nt
                 call mg_lev%uc(i)%setval(0.0_pfdp)
                 call mg_lev%uc_prev(i)%setval(0.0_pfdp)
              end do
           !end if
        else
           call pf_lev%ulevel%factory%create_array(mg_lev%uc, max(1, mg_lev%Nt), level_index, pf_lev%lev_shape)
           call pf_lev%ulevel%factory%create_array(mg_lev%uc_prev, max(1, mg_lev%Nt), level_index, pf_lev%lev_shape)
           do i = 1,mg_lev%Nt
              call mg_lev%uc(i)%setval(0.0_pfdp)
              call mg_lev%uc_prev(i)%setval(0.0_pfdp)
           end do
        end if
        call pf_lev%ulevel%factory%create_single(mg_lev%uc_ghost, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_single(mg_lev%u_temp, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_single(mg_lev%g_temp, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_single(mg_lev%r, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_single(mg_lev%Q0, level_index, pf_lev%lev_shape)
        if (FAS_flag .eqv. .true.) then
           call pf_lev%ulevel%factory%create_single(mg_lev%uc_fas_ghost, level_index, pf_lev%lev_shape)
        end if
     end do

     level_index = nlevels
     mg_lev => mg_ld(level_index)
     pf_lev = pf%levels(level_index)
     call pf_lev%ulevel%factory%create_array(mg_lev%e, max(1, mg_lev%Nt), level_index, pf_lev%lev_shape)

     if (mg_ld(nlevels)%coarsest_level .eq. nlevels) then
        print *,'Error in MGRIT setup: number of levels must be greater than 1'
        stop
     end if
  end subroutine mgrit_initialize

  !> Set up F and C points 
  subroutine FC_Setup(Nc, Nf, coarsen_factor, f_blocks, c_pts)
     type(int_vector), allocatable, intent(inout) :: f_blocks(:)
     integer, allocatable, intent(inout) :: c_pts(:)
     integer, intent(in) :: Nc, Nf, coarsen_factor
     integer, allocatable :: f_pts(:)
     integer :: i, j, n
     integer :: j_start, j_end, jj
     integer :: f_size

     n = coarsen_factor

     !> Find indices of C-points if this proc has any C-points
     if (Nc .gt. 0) then
        allocate(c_pts(Nc))
        j = min(Nf, n)
        do i = 1,Nc
           c_pts(i) = j;
           j = j + n;
        end do
     end if

     !> Find indices of F-points if this proc has any F-points
     if (Nf .gt. 0) then
        if (Nc .gt. 0) then
           allocate(f_blocks(Nc))
           do i = 1,Nc
              j_start = (i-1)*n+1;
              j_end = i*n-1;
              f_pts = [(jj, jj = j_start,j_end)]
              f_size = size(f_pts)
              allocate(f_blocks(i)%val(f_size))
              j = 1
              do jj = j_start,j_end
                 f_blocks(i)%val(j) = jj
                 j = j + 1
              end do
           end do
        else
           allocate(f_blocks(1))
           allocate(f_blocks(1)%val(1))
           f_blocks(1)%val(1) = 1
        end if
     end if
  end subroutine FC_Setup

  !>  Run the MGRIT solver
  subroutine pf_MGRIT_run(pf, mg_ld, Q0, qend)
    type(pf_pfasst_t), target, intent(inout)   :: pf   !!  The complete PFASST structure
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    class(pf_encap_t), intent(in   )           :: Q0   !!  The initial condition
    class(pf_encap_t), intent(inout), optional :: qend    !!  The computed solution at tend
    class(pf_encap_t), pointer :: uc(:)
    type(mgrit_level_data), pointer :: mg_lev
    type(pf_level_t), pointer :: pf_lev
    logical :: zero_rhs_flag
    integer :: iter, k, i, level_index, l
    integer :: ierror, source, tag
    integer :: Nt, Nt_glob
    integer :: uc_len
    real(pfdp) :: t0, dt, big_dt, t0_glob
    integer :: nlevels, coarsest_level
    double precision :: wtime_start
    integer :: stat(MPI_STATUS_SIZE)

    k = 1
    nlevels = pf%nlevels
    coarsest_level = mg_ld(nlevels)%coarsest_level
    pf%state%pfblock = 1
    pf%state%sweep = 1
    !pf%state%status = PF_STATUS_PREDICTOR
    !pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING

    mg_lev => mg_ld(nlevels)

    if (pf%niters .lt. 0) then
       uc => mg_ld(nlevels)%uc
       call qend%copy(uc(size(uc)))
       return
    end if

    !>  Allocate stuff for holding results 
    call initialize_results(pf)

    !> Restrict initial condition
    level_index = coarsest_level
    call mg_ld(nlevels)%Q0%copy(Q0)
    do l = nlevels,level_index+1,-1
       call pf%levels(l)%ulevel%restrict(pf%levels(l), pf%levels(l-1), mg_ld(l)%Q0, mg_ld(l-1)%Q0, mg_ld(l)%t0)
    end do

    !> Do predictor
    wtime_start = MPI_Wtime()
    call pf_start_timer(pf, T_TOTAL)
    call pf_start_timer(pf, T_PREDICTOR)
    mg_ld(nlevels)%cycle_phase = 0
    call PredictExactSolve(pf, mg_ld, Q0, level_index)
    zero_rhs_flag = .true.
    call IdealInterp(pf, mg_ld, 0, zero_rhs_flag)
    call pf_stop_timer(pf, T_PREDICTOR)

    if (pf%niters .lt. 1) then
       uc => mg_ld(nlevels)%uc
       call qend%copy(uc(size(uc)))
       return
    end if

    mg_ld(nlevels)%num_recv_converged = 0
    do i = 1,mg_ld(nlevels)%num_recvs
       mg_ld(nlevels)%recv_converge_status(i) = PF_STATUS_ITERATING
    end do

    mg_ld(nlevels)%c_converge_index = 1
    mg_ld(nlevels)%c_converge = .false.

    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING

    !> Iterate until convergence is achieved
    do iter = 1, pf%niters
       pf%state%iter = iter

       do level_index = coarsest_level,nlevels
          mg_lev => mg_ld(level_index)
          pf_lev => pf%levels(level_index)
          uc_len = size(mg_lev%uc)
          call pf_lev%delta_q0%copy(mg_lev%uc(uc_len))
       end do       

       !> Do a V-cycle
       call pf_start_timer(pf, T_ITERATION)
       call pf_MGRIT_v_cycle(pf, mg_ld, iter)
       call pf_set_iter(pf, iter)
       call pf_stop_timer(pf, T_ITERATION)

       do level_index = coarsest_level,nlevels
          mg_lev => mg_ld(level_index)
          pf_lev => pf%levels(level_index)
          uc_len = size(mg_lev%uc)
          call pf_lev%delta_q0%axpy(-1.0_pfdp, mg_lev%uc(uc_len))
          pf_lev%max_delta_q0 = pf_lev%delta_q0%norm()
          call pf_set_delta_q0(pf, level_index, pf_lev%delta_q0%norm())
          call pf_set_resid(pf, level_index, mg_lev%res_norm_loc)
          call pf_lev%qend%copy(mg_lev%uc(uc_len))
       end do

       call call_hooks(pf, nlevels, PF_POST_ITERATION)
       if (pf%state%status .eq. PF_STATUS_CONVERGED) then
          exit
       end if
    end do

    call pf_stop_timer(pf, T_TOTAL)
    !if (pf%rank .eq. pf%comm%nproc-1) print *,"solve time ",MPI_Wtime()-wtime_start

    call pf_dump_stats(pf)

    uc => mg_ld(nlevels)%uc
    call qend%copy(uc(size(uc)))
  end subroutine pf_MGRIT_run

  !> Execute a MGRIT V-cycle (iteration)
  subroutine pf_MGRIT_v_cycle(pf, mg_ld, iteration)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: iteration
    type(pf_level_t), pointer :: pf_lev, pf_f_lev, pf_c_lev
    type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev
    logical :: zero_rhs_flag, zero_c_pts_flag, interp_flag
    integer :: uc_len
    integer :: i, j, k, n, ierr
    integer :: nlevels, level_index, level_index_f
    integer :: ierror
    real(pfdp) :: res_norm_glob
    class(pf_encap_t), pointer :: uc(:)
    real(pfdp), pointer :: qend(:)
    real(pfdp) :: qexact, maxerr
    integer :: coarsest_level

    nlevels = pf%nlevels;
    coarsest_level = mg_ld(nlevels)%coarsest_level


    mg_ld(nlevels)%cycle_phase = 1
    level_index = nlevels
    mg_lev => mg_ld(level_index)

    !> FCF-relaxation on finest grid
    call pf_start_timer(pf, T_SWEEP, level_index)
    call FCF_Relax(pf, mg_ld, level_index, iteration)
    call pf_stop_timer(pf, T_SWEEP, level_index)

    do level_index = nlevels-1,coarsest_level+1,-1
        level_index_f = level_index+1
        mg_lev => mg_ld(level_index)
        mg_f_lev => mg_ld(level_index_f)

        !> FCF-relaxation on intermediate grids
        call pf_start_timer(pf, T_SWEEP, level_index)
        call FCF_Relax(pf, mg_ld, level_index, iteration)
        call pf_stop_timer(pf, T_SWEEP, level_index)
        
        if (pf%state%status .eq. PF_STATUS_CONVERGED) then
           return
        end if
    end do

    mg_ld(nlevels)%cycle_phase = 2
    level_index = coarsest_level
    level_index_f = level_index+1
    call Restrict(pf, mg_ld, level_index, level_index_f)

    if (pf%state%status .eq. PF_STATUS_CONVERGED) then
       return
    end if

    !> Coarsest grid solve
    call pf_start_timer(pf, T_SWEEP, level_index)
    call ExactSolve(pf, mg_ld, level_index)
    call pf_stop_timer(pf, T_SWEEP, level_index)

    !> Ideal interpolation
    zero_rhs_flag = .false.
    mg_ld(nlevels)%cycle_phase = 3
    call pf_start_timer(pf, T_SWEEP, level_index)
    call IdealInterp(pf, mg_ld, iteration, zero_rhs_flag)
    call pf_stop_timer(pf, T_SWEEP, level_index)
  end subroutine pf_MGRIT_v_cycle

  !> FCF-relaxation
  subroutine FCF_Relax(pf, mg_ld, level_index, iteration)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: iteration
     integer, intent(in) :: level_index
     type(pf_level_t), pointer :: pf_lev
     type(mgrit_level_data), pointer :: mg_lev, mg_c_lev
     logical :: zero_rhs_flag, zero_c_pts_flag, interp_flag, restrict_flag, send_flag, recv_flag, FC_relax_flag
     integer :: uc_len, i, nlevels
     integer :: ierror
     class(pf_encap_t), pointer :: uc(:)

     nlevels = pf%nlevels
     mg_lev => mg_ld(level_index)
     mg_c_lev => mg_ld(level_index-1)
     pf_lev => pf%levels(level_index)

     !> If not at finest level, initial values at C-points are zero
     if (level_index .lt. nlevels) then
        zero_c_pts_flag = .true.
     else
        zero_c_pts_flag = .false. 
     end if

     !> If not at finest level, right-hand side is non-zero
     if (level_index .lt. nlevels) then
        zero_rhs_flag = .false.
     else
        zero_rhs_flag = .true.
     end if

     !> Check if we are doing FCF-relaxations (just F-relaxations for Parareal)
     if (mg_lev%FCF_flag .eqv. .true.) then
        FC_relax_flag = .true.
     else
        FC_relax_flag = .false.
     end if

     !> Check if we are restricting
     if (level_index .eq. nlevels) then
        restrict_flag = .false.
     !> For the finest level, restrict at C-points first so that
     !> we can compute residuals and check for convergence
     else if (level_index .eq. nlevels-1) then
        restrict_flag = .false.
        call Restrict(pf, mg_ld, level_index, level_index+1)
     else
        restrict_flag = .true.
     end if
     if (pf%state%status .eq. PF_STATUS_CONVERGED) then
        return
     end if
     
     interp_flag = .false.

     uc_len = size(mg_lev%uc)
     do i = 1,uc_len
        call mg_lev%uc_prev(i)%copy(mg_lev%uc(i))
     end do

     !> F- or FC-relaxation
     call RelaxTransfer(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, restrict_flag, zero_c_pts_flag, FC_relax_flag) 

     !> F-relaxations if we are doing FCF-relaxations
     if (mg_lev%FCF_flag .eqv. .true.) then
        restrict_flag = .false.
        zero_c_pts_flag = .false.
        FC_relax_flag = .false.
        do i = 1,uc_len
           call mg_lev%uc_prev(i)%copy(mg_lev%uc(i))
        end do
        call RelaxTransfer(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, restrict_flag, zero_c_pts_flag, FC_relax_flag)
     end if
  end subroutine FCF_Relax

  !> Ideal interpolation, i.e., interpolation followed by F-relaxations
  subroutine IdealInterp(pf, mg_ld, iteration, zero_rhs_flag)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: iteration
    logical, intent(in) :: zero_rhs_flag
    type(pf_level_t), pointer :: pf_lev, pf_f_lev, pf_c_lev
    type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev
    logical :: zero_c_pts_flag, interp_flag, restrict_flag, send_flag, recv_flag, FC_relax_flag
    integer :: uc_len
    integer :: i, j, k, n, ierr
    integer :: nlevels, level_index, level_index_f, level_index_c
    class(pf_encap_t), pointer :: uc(:)
    integer :: coarsest_level

    nlevels = pf%nlevels;
    coarsest_level = mg_ld(nlevels)%coarsest_level
    zero_c_pts_flag = .false.
    interp_flag = .true.
    restrict_flag = .false.
    FC_relax_flag = .false.
    do level_index = coarsest_level+1,(nlevels-1)
       level_index_f = level_index+1
       mg_lev => mg_ld(level_index)
       pf_lev => pf%levels(level_index)
       mg_f_lev => mg_ld(level_index_f)
       pf_f_lev => pf%levels(level_index_f)
       uc_len = size(mg_lev%uc)
       do i = 1,uc_len
          call mg_lev%uc_prev(i)%copy(mg_lev%uc(i))
       end do
       !> F-relaxation followed by interpolation
       call RelaxTransfer(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, restrict_flag, zero_c_pts_flag, FC_relax_flag) 
    end do
  end subroutine IdealInterp

  !> Do a series of relxations followed by an optional grid transfer (restriction or interpolation)
  subroutine RelaxTransfer(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, restrict_flag, zero_c_pts_flag, FC_relax_flag)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     logical, intent(in) :: zero_rhs_flag, zero_c_pts_flag, interp_flag, restrict_flag, FC_relax_flag
     type(mgrit_level_data), pointer :: mg_f_lev, mg_lev
     type(pf_level_t), pointer :: pf_lev, pf_f_lev
     integer :: i, j, n, ii, nlevels, level_index_f, i_upper, uc_len, i_f

     nlevels = pf%nlevels
     level_index_f = level_index + 1
     mg_lev => mg_ld(level_index)
     mg_f_lev => mg_ld(level_index_f)
     pf_lev => pf%levels(level_index)
     pf_f_lev => pf%levels(level_index_f)

     !> In this case of one time point on the fine level, values at C-points need to be sent in advance for restriction
     if ((restrict_flag .eqv. .true.) .and. (mg_f_lev%Nt .eq. 1)) then
        if (mg_f_lev%f_pts_flag .eqv. .true.) then
           call mgrit_send(pf, mg_ld, mg_f_lev%uc(1), level_index_f, 1, .false.)
        else
           call mgrit_post(pf, mg_ld, level_index_f, 1)
           call mgrit_recv(pf, mg_ld, mg_f_lev%uc_ghost, level_index_f, 1, .false.)
           call mg_f_lev%uc(1)%copy(mg_f_lev%uc_ghost)
        end if
     end if

     !> Return if no time points
     if (mg_lev%Nt .lt. 1) then
        return
     end if

     !> Post non-blocking receives
     if (FC_relax_flag .eqv. .true.) then
        if ((zero_c_pts_flag .eqv. .false.) .or. ((mg_lev%Nt .eq. 1) .and. (mg_lev%c_pts_flag .eqv. .true.))) then
           call mgrit_post(pf, mg_ld, level_index, 1)
        end if
     else
        if ((zero_c_pts_flag .eqv. .false.) .and. (mg_lev%f_pts_flag .eqv. .true.)) then
           call mgrit_post(pf, mg_ld, level_index, 1)
        end if
     end if

     if (mg_lev%f_pts_flag .eqv. .true.) then
        i_upper = size(mg_lev%f_blocks)
     else
        i_upper = 1
     end if

     uc_len = size(mg_lev%uc)
     !> Start at right-most F-block such that computation and communication can be overlapped
     do i = i_upper,1,-1
     !do i = 1,i_upper
        !> If this proc has C-points, send the right-most C-point
        if ((i .eq. i_upper) .and. (zero_c_pts_flag .eqv. .false.) .and. (mg_lev%c_pts_flag .eqv. .true.)) then
           call mgrit_send(pf, mg_ld, mg_lev%uc(uc_len), level_index, 1, .false.)
        end if
        !> If this proc has F-points, do F-relaxations
        if (mg_lev%f_pts_flag .eqv. .true.) then
           !> Set initial condition for this block
           if ((i .eq. 1) .and. (mg_lev%rank_shifted .eq. 0)) then
              if ((level_index .eq. nlevels) .or. (mg_ld(nlevels)%cycle_phase .eq. 0)) then
                 call pf_lev%q0%copy(mg_lev%Q0)
              else if (mg_ld(nlevels)%FAS_flag .eqv. .true.) then
                 if (zero_c_pts_flag .eqv. .true.) then
                    call pf_lev%q0%setval(0.0_pfdp)
                 else
                    call pf_lev%q0%copy(mg_lev%Q0)
                 end if
              else
                 call pf_lev%q0%setval(0.0_pfdp)
              end if
           else
              if (zero_c_pts_flag .eqv. .true.) then
                 call pf_lev%q0%setval(0.0_pfdp)
              else if (i .eq. 1) then
                 !> Receive off-proc C-point
                 call mgrit_recv(pf, mg_ld, mg_lev%uc_ghost, level_index, 1, .false.)
                 call pf_lev%q0%copy(mg_lev%uc_ghost)
              else
                 call pf_lev%q0%copy(mg_lev%uc_prev(i-1))
              end if
           end if
       
           do n = 1,size(mg_lev%f_blocks(i)%val)
              j = mg_lev%f_blocks(i)%val(n)
              !> Restrict at this point
              if (restrict_flag .eqv. .true.) then
                 i_f = mg_f_lev%c_pts(j)
                 call InjectRestrictPoint(pf, mg_ld, mg_lev%g(j), level_index, level_index_f, j, i_f, zero_rhs_flag);
                 call mg_f_lev%uc(j)%copy(mg_f_lev%uc_prev(j))
              end if
              !> Do a single step
              call PointRelax(pf, mg_ld, level_index, j, pf_lev%q0, pf_lev%qend)
              !> Add right-hand side to the result
              if (zero_rhs_flag .eqv. .false.) then
                 call pf_lev%qend%axpy(1.0_pfdp, mg_lev%g(j))
              end if

              !> Interpolate at this F-point
              if (interp_flag .eqv. .true.) then
                call InjectInterpPoint(pf, mg_ld, pf_lev%qend, level_index, level_index_f, j, j)
              end if
              call pf_lev%q0%copy(pf_lev%qend)
           end do
           !> If this proc has an F-point but no C-point, send F-point forward
           if ((FC_relax_flag .eqv. .true.) .and. (mg_lev%c_pts_flag .eqv. .false.)) then
              call mgrit_send(pf, mg_ld, pf_lev%qend, level_index, 1, .false.)
           end if
        end if

        !> Do C-relaxation
        if (FC_relax_flag .eqv. .true.) then
           !> If this proc has C-points do C-relaxations
           if (mg_lev%c_pts_flag .eqv. .true.) then
              j = mg_lev%c_pts(i)
              !> Restrict at this C-point
              if (restrict_flag .eqv. .true.) then
                 i_f = mg_f_lev%c_pts(j)
                 call InjectRestrictPoint(pf, mg_ld, mg_lev%g(j), level_index, level_index_f, j, i_f, zero_rhs_flag);
                 call mg_f_lev%uc(j)%copy(mg_f_lev%uc_prev(j))
              end if
              !> If there is only one time point, receive F-point
              if (mg_lev%f_pts_flag .eqv. .false.) then
                 call mgrit_recv(pf, mg_ld, mg_lev%uc_ghost, level_index, 1, .false.)
                 call pf_lev%q0%copy(mg_lev%uc_ghost)
              end if
              call PointRelax(pf, mg_ld, level_index, j, pf_lev%q0, pf_lev%qend)
              call mg_lev%uc(i)%copy(pf_lev%qend)
              if (zero_rhs_flag .eqv. .false.) then
                 call mg_lev%uc(i)%axpy(1.0_pfdp, mg_lev%g(j))
              end if
           end if
        else if (mg_lev%f_pts_flag .eqv. .true.) then
           call mg_lev%uc(i)%copy(pf_lev%qend)
        end if

        !> Interpolate at this C-point
        if ((interp_flag .eqv. .true.) .and. (mg_lev%c_pts_flag .eqv. .true.)) then
           if (mg_lev%f_pts_flag .eqv. .false.) then
              j = 1
           else
              j = j + 1
           end if
           call InjectInterpPoint(pf, mg_ld, mg_lev%uc_prev(i), level_index, level_index_f, j, j)
           call mg_lev%uc(i)%copy(mg_lev%uc_prev(i))
        end if
     end do
  end subroutine RelaxTransfer

  !> Blocking send-receive
  subroutine mgrit_send_recv(pf, mg_ld, level_index, send_flag, recv_flag)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: level_index
    logical, intent(in) :: send_flag, recv_flag
    type(pf_level_t), pointer :: pf_lev
    integer :: ierror
    type(mgrit_level_data), pointer :: mg_lev

    pf_lev => pf%levels(level_index)
    mg_lev => mg_ld(level_index)
    if (mod(mg_lev%rank_shifted, 2) .eq. 0) then
       if (recv_flag .eqv. .true.) then
          call mgrit_recv(pf, mg_ld, mg_lev%uc_ghost, level_index, 1, .true.)
       end if
       if (send_flag .eqv. .true.) then
          call mgrit_send(pf, mg_ld, pf_lev%qend, level_index, 1, .true.)
       end if 
    else
       if (send_flag .eqv. .true.) then
          call mgrit_send(pf, mg_ld, pf_lev%qend, level_index, 1, .true.)
       end if
       if (recv_flag .eqv. .true.) then
          call mgrit_recv(pf, mg_ld, mg_lev%uc_ghost, level_index, 1, .true.)
       end if
    end if
  end subroutine mgrit_send_recv

  !> Exact solve for coarsest grid
  subroutine ExactSolve(pf, mg_ld, level_index)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     integer :: level_index_f
     type(pf_level_t), pointer :: lev
     class(pf_encap_t), allocatable :: gi
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: zero_rhs_flag, send_flag, recv_flag
     type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_finest_lev
     type(pf_level_t), pointer :: pf_lev, pf_f_lev

     mg_lev => mg_ld(level_index)
     nlevels = pf%nlevels
     level_index_f = level_index + 1

     mg_f_lev => mg_ld(level_index_f)
     mg_finest_lev => mg_ld(nlevels)

     pf_lev => pf%levels(level_index)
     pf_f_lev => pf%levels(level_index_f)

     if (mg_lev%Nt .lt. 1) then
        return
     end if

     if ((mg_lev%rank_shifted .gt. 0) .and. (pf%comm%nproc .gt. 1)) then
        call mgrit_post(pf, mg_ld, level_index, 1)
        call mgrit_recv(pf, mg_ld, mg_lev%uc_ghost, level_index, 1, .false.)
        call pf_lev%q0%copy(mg_lev%uc_ghost)
     else
        if (mg_ld(nlevels)%FAS_flag .eqv. .true.) then
           call pf_lev%q0%copy(mg_lev%Q0)
        else
           call pf_lev%q0%setval(0.0_pfdp);
        end if
     end if

     do i = 1,mg_lev%Nt
        call PointRelax(pf, mg_ld, level_index, i, pf_lev%q0, pf_lev%qend)
        call pf_lev%qend%axpy(1.0_pfdp, mg_lev%g(i))
        call pf_lev%q0%copy(pf_lev%qend)
        call mg_lev%uc(i)%copy(pf_lev%qend)
     end do

     if ((pf%rank .lt. pf%comm%nproc-1) .and. (pf%comm%nproc .gt. 1)) then
        call mgrit_send(pf, mg_ld, pf_lev%qend, level_index, 1, .false.)
     end if

     do i = 1,mg_lev%Nt
        call InjectInterpPoint(pf, mg_ld, mg_lev%uc(i), level_index, level_index_f, i, i)
     end do
  end subroutine ExactSolve

  !> Coarsest grid exact solver used by predictor
  subroutine PredictExactSolve(pf, mg_ld, Q0, level_index)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     class(pf_encap_t), intent(in) :: Q0
     integer, intent(in) :: level_index
     integer :: level_index_f
     type(pf_level_t), pointer :: lev
     integer :: i, j, k, n, ii, l
     integer :: nlevels
     type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_finest_lev
     type(pf_level_t), pointer :: pf_lev, pf_f_lev


     mg_lev => mg_ld(level_index)
     nlevels = pf%nlevels
     level_index_f = level_index + 1
     mg_f_lev => mg_ld(level_index_f)
     mg_finest_lev => mg_ld(nlevels)

     pf_lev => pf%levels(level_index)
     pf_f_lev => pf%levels(level_index_f)

     if (mg_lev%Nt .lt. 1) then
        return
     end if

     if ((mg_lev%rank_shifted .gt. 0) .and. (pf%comm%nproc .gt. 1)) then
        call mgrit_post(pf, mg_ld, level_index, 3)
        call mgrit_recv(pf, mg_ld, pf_lev%q0, level_index, 3, .false.)
     else
        call pf_lev%q0%copy(mg_lev%Q0)
     end if

     do i = 1,mg_lev%Nt
        call PointRelax(pf, mg_ld, level_index, i, pf_lev%q0, pf_lev%qend)
        call pf_lev%q0%copy(pf_lev%qend)
        call mg_lev%uc(i)%copy(pf_lev%qend)
     end do

     call mgrit_send(pf, mg_ld, pf_lev%qend, level_index, 3, .false.)

     do i = 1,mg_lev%Nt
        call pf_start_timer(pf, T_INTERPOLATE, level_index_f)
        call pf_lev%ulevel%interpolate(pf_f_lev, pf_lev, mg_f_lev%uc(i), mg_lev%uc(i), mg_f_lev%t0)
        call pf_stop_timer(pf, T_INTERPOLATE, level_index_f)
     end do
  end subroutine PredictExactSolve

  !> Restrict residuals corresponding to C-points
  subroutine Restrict(pf, mg_ld, level_index_c, level_index_f)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index_c, level_index_f
     class(pf_encap_t), pointer :: gi
     integer :: i_c, i_f, k, n
     integer :: nlevels, level_index
     logical :: zero_rhs_flag, send_flag, recv_flag
     real(pfdp) :: r_norm
     type(mgrit_level_data), pointer :: mg_lev, mg_c_lev, mg_f_lev
     type(pf_level_t), pointer :: pf_lev, pf_f_lev
     real(pfdp) :: res_norm_glob
     integer :: ierr, i_upper

     nlevels = pf%nlevels

     mg_f_lev => mg_ld(level_index_f)
     if (mg_f_lev%Nt .lt. 1) then
        return
     end if

     mg_c_lev => mg_ld(level_index_c)
     pf_f_lev => pf%levels(level_index_f)

     !> In this case of one time point on the fine level, values at C-points need to be sent in advance for restriction
     if (mg_f_lev%Nt .eq. 1) then
        if (mg_f_lev%f_pts_flag .eqv. .true.) then
           call mgrit_send(pf, mg_ld, mg_f_lev%uc(1), level_index_f, 1, .false.) 
        else
           call mgrit_post(pf, mg_ld, level_index_f, 1)
           call mgrit_recv(pf, mg_ld, mg_f_lev%uc_ghost, level_index_f, 1, .false.)
           call mg_f_lev%uc(1)%copy(mg_f_lev%uc_ghost)
        end if
     end if

     mg_f_lev%res_norm_loc = 0.0_pfdp
     i_upper = mg_c_lev%Nt

     !> Restrict at C-points
     !do i_c = i_upper,1,-1
     do i_c = 1,i_upper
        i_f = mg_f_lev%c_pts(i_c)
        if (level_index_f .eq. nlevels) then
           zero_rhs_flag = .true.
        else
           zero_rhs_flag = .false.
        end if
        call InjectRestrictPoint(pf, mg_ld, mg_c_lev%g(i_c), level_index_c, level_index_f, i_c, i_f, zero_rhs_flag);
     end do

     do i_c = 1,i_upper
        call mg_f_lev%uc(i_c)%copy(mg_f_lev%uc_prev(i_c))
     end do

     !> If the finest level is finest, check residual norm
     if (level_index_f .eq. nlevels) then
!        if () then
!           call mpi_allreduce(mg_f_lev%res_norm_loc, res_norm_glob, 1, myMPI_Datatype, MPI_SUM, pf%comm%comm, ierr)
!           res_norm_glob = sqrt(res_norm_glob)
!           pf%levels(nlevels)%residual = res_norm_glob
!           if (res_norm_glob .lt. pf%abs_res_tol) then
!              pf%state%pstatus = PF_STATUS_CONVERGED
!           end if
!        else 
           pf%levels(nlevels)%residual = sqrt(mg_f_lev%res_norm_loc)
           call mgrit_check_convergence(mg_ld, pf, nlevels, 0)
           do level_index = mg_ld(nlevels)%coarsest_level,(nlevels-1)
              mg_lev => mg_ld(level_index)
              call mg_lev%uc_ghost%setval(0.0_pfdp)
           end do
!        end if
     end if
  end subroutine Restrict

  !> Interpolate correction using injection
  subroutine InjectInterpPoint(pf, mg_ld, eci, level_index_c, level_index_f, i_c, i_f)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     class(pf_encap_t), intent(inout) :: eci
     integer, intent(in) :: level_index_c, level_index_f, i_c, i_f
     type(mgrit_level_data), pointer :: mg_c_lev, mg_f_lev
     type(pf_level_t), pointer :: pf_f_lev, pf_c_lev
     integer :: nlevels

     mg_c_lev => mg_ld(level_index_c)
     mg_f_lev => mg_ld(level_index_f)

     pf_f_lev => pf%levels(level_index_f)
     pf_c_lev => pf%levels(level_index_c)
     nlevels = pf%nlevels

     call pf_start_timer(pf, T_INTERPOLATE, level_index_f)
     call mg_c_lev%u_temp%copy(eci)
     if ((mg_ld(nlevels)%FAS_flag .eqv. .true.) .and. (mg_ld(nlevels)%cycle_phase .gt. 0)) then
        call mg_c_lev%u_temp%axpy(-1.0_pfdp, mg_c_lev%uc_fas(i_f))
     end if
     !> Spatial interpolation
     call pf_f_lev%ulevel%interpolate(pf_f_lev, pf_c_lev, mg_f_lev%u_temp, mg_c_lev%u_temp, mg_f_lev%t0)
     call mg_f_lev%uc(i_f)%axpy(1.0_pfdp, mg_f_lev%u_temp)
     call pf_stop_timer(pf, T_INTERPOLATE, level_index_f)
  end subroutine InjectInterpPoint

  !> Restrict residual using injection
  subroutine InjectRestrictPoint(pf, mg_ld, gci, level_index_c, level_index_f, i_c, i_f, zero_rhs_flag)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     class(pf_encap_t), intent(inout) :: gci
     integer, intent(in) :: level_index_c, level_index_f, i_c, i_f
     logical, intent(in) :: zero_rhs_flag
     type(mgrit_level_data), pointer :: mg_c_lev, mg_f_lev
     type(pf_level_t), pointer :: pf_f_lev, pf_c_lev
     integer :: nlevels
     logical :: send_flag, recv_flag

     mg_c_lev => mg_ld(level_index_c)
     mg_f_lev => mg_ld(level_index_f)

     pf_f_lev => pf%levels(level_index_f)
     pf_c_lev => pf%levels(level_index_c)
     nlevels = pf%nlevels

     call mg_f_lev%u_temp%copy(mg_f_lev%uc(i_c))
     !> Take step
     call PointRelax(pf, mg_ld, level_index_f, i_f, mg_f_lev%u_temp, mg_f_lev%r)
     !> Finish computing residual
     call mg_f_lev%r%axpy(-1.0_pfdp, mg_f_lev%uc_prev(i_c))
     if (zero_rhs_flag .eqv. .false.) then
        call mg_f_lev%r%axpy(1.0_pfdp, mg_f_lev%g(i_f));
     end if
     !> If on the finest level, update res norm
     if (level_index_f .eq. nlevels) then
        call SumResNorm(pf, mg_ld, level_index_f, mg_f_lev%res_norm_loc, mg_f_lev%r, i_c, 2)
     end if
     call pf_start_timer(pf, T_RESTRICT, level_index_f)
     !> Spatial restriction
     call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%r, mg_c_lev%r, mg_f_lev%t0)
     !> FAS restriction
     if (mg_ld(nlevels)%FAS_flag .eqv. .true.) then
        !> FAS spatial restriction  
        if (i_c .eq. 1) then !> Case of left-most C-point
           !> Restrict off-proc C-point
           if (mg_c_lev%rank_shifted .eq. 0) then
              call mg_c_lev%u_temp%copy(mg_c_lev%Q0)
           else
              if (mg_c_lev%Nt .eq. 1) then
                 call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%uc_prev(i_c), mg_c_lev%uc_fas(i_c), mg_f_lev%t0)
              end if
              call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%uc_ghost, mg_c_lev%u_temp, mg_f_lev%t0)
           end if
        else if (i_c .eq. mg_c_lev%Nt) then !> Case of right-most C-point
           !> Since we visit C-points from right to left, restrict at current and left C-point 
           call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%uc_prev(i_c), mg_c_lev%uc_fas(i_c), mg_f_lev%t0)
           call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%uc_prev(i_c-1), mg_c_lev%uc_fas(i_c-1), mg_f_lev%t0)
           call mg_c_lev%u_temp%copy(mg_c_lev%uc_fas(i_c-1))
        else
           !> Restrict at C-point to the left
           call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%uc_prev(i_c-1), mg_c_lev%uc_fas(i_c-1), mg_f_lev%t0)
           call mg_c_lev%u_temp%copy(mg_c_lev%uc_fas(i_c-1))
        end if
        !> Update RHS using approximation at C-point
        call mg_c_lev%r%axpy(1.0_pfdp, mg_c_lev%uc_fas(i_c))
        call PointRelax(pf, mg_ld, level_index_c, i_c, mg_c_lev%u_temp, pf_c_lev%qend)
        call mg_c_lev%r%axpy(-1.0_pfdp, pf_c_lev%qend)
     end if
     call pf_stop_timer(pf, T_RESTRICT, level_index_f)
     call gci%copy(mg_c_lev%r)
  end subroutine InjectRestrictPoint

  !> Take a step
  subroutine PointRelax(pf, mg_ld, level_index, n, q0, qend)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index, n
     class(pf_encap_t), intent(inout) :: q0, qend
     type(pf_level_t), pointer :: pf_lev
     type(mgrit_level_data), pointer :: mg_lev
     real(pfdp) :: t0n

     mg_lev => mg_ld(level_index)
     pf_lev => pf%levels(level_index)

     t0n = mg_lev%t0 + mg_lev%dt * real(n-1,pfdp)
     call pf_lev%ulevel%stepper%do_n_steps(pf, level_index, t0n, q0, qend, mg_lev%dt, 1)
  end subroutine PointRelax

  !> Compute norm of r and add it to sum_norms (sum of norms)
  subroutine SumResNorm(pf, mg_ld, level_index, sum_norms, r, i, norm_type)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index, i
     real(pfdp), intent(inout) :: sum_norms
     class(pf_encap_t), intent(in) :: r
     integer, intent(in) :: norm_type
     real(pfdp) :: r_norm
     type(mgrit_level_data), pointer :: mg_lev, mg_c_lev

     mg_lev => mg_ld(level_index)
     mg_c_lev => mg_ld(level_index-1)

     r_norm = r%norm()

     mg_lev%c_converge = .false.
     !if (pf%rank .eq. 1) print *,i,mg_lev%c_converge_index,mg_c_lev%Nt,r_norm
     if ((i .eq. mg_lev%c_converge_index) .and. (r_norm .lt. pf%abs_res_tol))  then
        if (mg_lev%c_converge_index .eq. mg_c_lev%Nt) then
           mg_lev%c_converge = .true.
        else
           mg_lev%c_converge_index = mg_lev%c_converge_index + 1
        end if
     end if

     if (norm_type .eq. 0) then
        sum_norms = max(sum_norms, r_norm)
     else if (norm_type .eq. 1) then
        sum_norms = sum_norms + r_norm
     else if (norm_type == 2) then
        sum_norms = sum_norms + r_norm*r_norm
     end if
  end subroutine SumResNorm

  subroutine mgrit_check_convergence(mg_ld, pf, level_index, tag)
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index
    integer, intent(in) :: tag  !! identifier for status send and receive
    logical :: residual_converged, converged, send_flag
    type(mgrit_level_data), pointer :: mg_lev
    integer :: ierr, pstatus, istatus, source, dest

    mg_lev => mg_ld(level_index)

    !> Check to see if tolerances are met
    residual_converged = .false.
    if (mg_lev%c_converge)  then
       residual_converged = .true.
    end if    

    !>  Until I hear the previous processor is done, recieve it's status
    if (pf%state%pstatus /= PF_STATUS_CONVERGED) then
       call mgrit_recv_status(pf, mg_ld, level_index, tag)
    end if

    !>  Check to see if I am converged
    converged = .false.
    if (residual_converged) then
       if (pf%rank == 0) then
          converged = .true.
       else  !  I am not the first processor, so I need to check the previous one
          if (pf%state%pstatus == PF_STATUS_CONVERGED) converged = .true.
       end if
    end if ! (residual_converged)

    if ((mg_lev%Nt * (pf%rank+1)) .lt. pf%state%iter) then
       converged = .true.
    end if

    send_flag = .false.
    !> Assign status and send it forward
    if (converged) then
       if (pf%state%status == PF_STATUS_ITERATING) then
          !  If I am converged for the first time
          !  then flip my flag and send the last status update
          pf%state%status = PF_STATUS_CONVERGED
          send_flag = .true.
       end if
    else
       !  I am not converged, send the news
       pf%state%status = PF_STATUS_ITERATING
       send_flag = .true.
    end if

    if (send_flag) call mgrit_send_status(pf, mg_ld, level_index, tag)
  end subroutine mgrit_check_convergence

  !> Non-blocking MPI send
  subroutine mgrit_send(pf, mg_ld, send_data, level_index, tag, blocking)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    class(pf_encap_t), intent(in) :: send_data
    integer, intent(in) :: level_index, tag
    logical, intent(in) :: blocking
    type(pf_level_t), pointer :: pf_lev
    integer :: ierror, dir, which, stride
    type(mgrit_level_data), pointer :: mg_lev

    if ((pf%rank .lt. pf%comm%nproc-1) .and. (pf%comm%nproc .gt. 1)) then
       pf_lev => pf%levels(level_index)

       dir = 1
       which = 1
       stride = mg_ld(level_index)%send_to_rank - pf%rank
       call pf_lev%qend%copy(send_data)
       call pf_send(pf, pf_lev, tag, blocking, dir, which, stride)
    end if
  end subroutine mgrit_send

  !> Non-blocking MPI receive
  subroutine mgrit_recv(pf, mg_ld, recv_data, level_index, tag, blocking)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    class(pf_encap_t), intent(inout) :: recv_data
    integer, intent(in) :: level_index, tag
    logical, intent(in) :: blocking
    type(pf_level_t), pointer :: pf_lev
    integer :: ierror, dir, which, stride
    type(mgrit_level_data), pointer :: mg_lev

    if (pf%state%pstatus /= PF_STATUS_CONVERGED) then
       mg_lev => mg_ld(level_index)
       if ((mg_lev%rank_shifted .gt. 0) .and. (pf%comm%nproc .gt. 1)) then
          pf_lev => pf%levels(level_index)

          dir = 1
          which = 1
          stride = pf%rank - mg_lev%recv_from_rank
          call pf_recv(pf, pf_lev, tag, blocking, dir, which, stride)
          call recv_data%copy(pf_lev%q0)
       end if
    end if
  end subroutine mgrit_recv

  !> MPI post used in non-blocking communication
  subroutine mgrit_post(pf, mg_ld, level_index, tag)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: level_index, tag
    type(pf_level_t), pointer :: pf_lev
    integer :: ierror, dir, stride
    type(mgrit_level_data), pointer :: mg_lev

    if (pf%state%pstatus /= PF_STATUS_CONVERGED) then
       mg_lev => mg_ld(level_index)
       if ((mg_lev%rank_shifted .gt. 0) .and. (pf%comm%nproc .gt. 1)) then
          pf_lev => pf%levels(level_index)

          dir = 1
          stride = pf%rank - mg_lev%recv_from_rank
          call pf_post(pf, pf_lev, tag, dir, stride)
       end if
    end if
  end subroutine mgrit_post


  !> Non-blocking MPI send
  subroutine mgrit_send_status(pf, mg_ld, level_index, tag)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: level_index, tag
    type(mgrit_level_data), pointer :: mg_lev
    integer :: ierr, istatus, dest, i

    if (pf%rank .lt. pf%comm%nproc-1) then
       mg_lev => mg_ld(pf%nlevels)
       istatus = pf%state%status
       do i = 1,mg_lev%num_sends
          dest = mg_lev%send_procs(i)
          call pf%comm%send_status(pf, tag, istatus, ierr, dest)
       end do
    end if
  end subroutine mgrit_send_status

  !> Non-blocking MPI receive
  subroutine mgrit_recv_status(pf, mg_ld, level_index, tag)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: level_index, tag
    type(mgrit_level_data), pointer :: mg_lev
    integer :: i, ierr, pstatus, source

    if (pf%rank .gt. 0) then
       mg_lev => mg_ld(pf%nlevels)
       do i = 1,mg_lev%num_recvs
          if (mg_lev%recv_converge_status(i) .ne. PF_STATUS_CONVERGED) then
             source = mg_lev%recv_procs(i)
             call pf%comm%recv_status(pf, tag, pstatus, ierr, source)
             mg_lev%recv_converge_status(i) = pstatus
             if (pstatus .eq. PF_STATUS_CONVERGED) then
                 mg_lev%num_recv_converged = mg_lev%num_recv_converged + 1
             end if
          end if
       end do
       if (mg_lev%num_recv_converged .eq. mg_lev%num_recvs) then
          pf%state%pstatus = PF_STATUS_CONVERGED
       end if
    end if
  end subroutine mgrit_recv_status

end module pf_mod_MGRIT
