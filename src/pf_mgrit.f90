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
  use mpi
  
  implicit none

  type :: int_vector
     integer, allocatable :: val(:)
  end type int_vector

  type :: mgrit_level_data
     integer :: Nt
     real(pfdp) :: dt
     real(pfdp) :: glob_t0
     real(pfdp) :: glob_tfin
     real(pfdp) :: t0
     real(pfdp) :: tfin
     real(pfdp) :: res_norm_loc(1)
     integer :: res_norm_index
     integer :: cycle_phase
     logical :: FCF_flag = .true.
     logical :: FAS_flag = .false.
     logical :: setup_start_coarse = .true.
     type(int_vector), allocatable :: f_blocks(:)
     integer, allocatable :: c_pts(:)
     class(pf_encap_t), allocatable :: g(:)
     class(pf_encap_t), allocatable :: qc(:)
     class(pf_encap_t), allocatable :: qc_prev(:)
     class(pf_encap_t), allocatable :: qc_fas(:)
     class(pf_encap_t), allocatable :: r
     class(pf_encap_t), allocatable :: qc_boundary
     class(pf_encap_t), allocatable :: q_temp
     class(pf_encap_t), allocatable :: Q0
     integer, allocatable :: send_procs(:)
     integer, allocatable :: recv_procs(:)
     integer, allocatable :: level_ranks_shifted(:)
     integer, allocatable :: level_ranks(:)
     integer :: nprocs_level
  end type 
  
contains

  subroutine mgrit_initialize(pf, mg_ld, glob_t0, glob_tfin, Nt_start, coarsen_factor, FAS_flag, FCF_flag)
     type(pf_pfasst_t), intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: Nt_start, coarsen_factor
     logical, intent(in) :: FAS_flag
     logical, intent(in) :: FCF_flag
     real(pfdp) :: glob_t0, glob_tfin
     integer :: level_index, nlevels, N, i, kk, p, Nt, p_rank, k
     type(pf_level_t) :: pf_lev, pf_f_lev, pf_c_lev
     type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev
     class(pf_encap_t), allocatable :: qtemp
     integer :: count_levels
     logical :: start_dropping_procs
     integer :: my_rank, rank_send, rank_recv, Nt_temp
     integer :: Nt_coarsest

     Nt_coarsest = 2

     nlevels = pf%nlevels
     n = coarsen_factor

     allocate(mg_ld(nlevels))
     do level_index = 1,nlevels
        allocate(mg_ld(level_index)%send_procs(pf%comm%nproc))
        allocate(mg_ld(level_index)%recv_procs(pf%comm%nproc))
     end do
     if (mg_ld(nlevels)%setup_start_coarse .eqv. .true.) then
        N = Nt_start;
        do level_index = 1,nlevels
           mg_lev => mg_ld(level_index)
           mg_lev%res_norm_loc(1) = 0.0_pfdp

           mg_lev%glob_t0 = glob_t0
           mg_lev%glob_tfin = glob_tfin
           mg_lev%Nt = N
           mg_lev%dt = (glob_tfin - glob_t0)/(N * pf%comm%nproc)
           mg_lev%t0 = glob_t0 + N * mg_lev%dt * pf%rank + mg_lev%dt
           pf%state%t0 = mg_lev%t0
           mg_lev%tfin = mg_lev%t0 + N * mg_lev%dt - mg_lev%dt
           N = N * coarsen_factor

           if (level_index .eq. nlevels) then
              pf%state%t0 = mg_lev%tfin - mg_lev%dt
              pf%state%dt = mg_lev%dt
              pf%state%step = (pf%rank+1) * mg_lev%Nt - 1
           end if
           do p = 1,pf%comm%nproc
              if (p .lt. pf%comm%nproc) then
                 mg_lev%send_procs(p) = p
              end if
              if (p .gt. 1) then
                 mg_lev%recv_procs(p) = p-2
              end if
           end do
        end do
     else
        Nt_temp = Nt_start * pf%comm%nproc
        count_levels = 1
        do while (count_levels < pf%nlevels)
           Nt_temp = Nt_temp / coarsen_factor;
           if (Nt_temp < Nt_coarsest) then
              exit
           end if
           count_levels = count_levels + 1;
        end do
        pf%nlevels = count_levels
        nlevels = pf%nlevels
        mg_ld(nlevels)%Nt = Nt_start
        start_dropping_procs = .false.
        i = 1
        do level_index = (nlevels-1),1,-1
           mg_lev => mg_ld(level_index)
           mg_f_lev => mg_ld(level_index+1)
           mg_lev%nprocs_level = 0
           do p = 1,pf%comm%nproc
              p_rank = p-1
              if (mg_f_lev%Nt .eq. 1) then
                 if (start_dropping_procs .eqv. .false.) then
                    start_dropping_procs = .true.
                 end if
                 if (mod(p, coarsen_factor**(i+1)) .eq. 0) then 
                    mg_lev%Nt = mg_f_lev%Nt
                    mg_lev%nprocs_level = mg_lev%nprocs_level + 1
                 else
                    mg_lev%Nt = 0
                 end if
              else if (mg_f_lev%Nt .gt. 1) then
                 mg_lev%Nt = mg_f_lev%Nt / coarsen_factor
                 mg_lev%nprocs_level = mg_lev%nprocs_level + 1
              else
                 mg_lev%Nt = 0
              end if
           end do
           if (start_dropping_procs .eqv. .true.) then
              i = i + 1;
           end if
        end do


        do level_index = nlevels,1,-1
           mg_lev => mg_ld(level_index)
           allocate(mg_lev%level_ranks_shifted(pf%comm%nproc))
           allocate(mg_lev%level_ranks(mg_lev%nprocs_level))
           k = 1
           do p = 1,pf%comm%nproc
              p_rank = p-1
              mg_lev%level_ranks_shifted(p) = -1
              if (mg_lev%Nt .gt. 0) then
                 mg_lev%level_ranks_shifted(p) = k-1
                 mg_lev%level_ranks(k) = p_rank
                 k = k + 1
              end if 
           end do
           my_rank = mg_lev%level_ranks(1);
           rank_recv = mg_lev%level_ranks(2);
           mg_lev%send_procs(my_rank+1) = rank_recv;
           do i = 2,(mg_lev%nprocs_level-1)
               my_rank = mg_lev%level_ranks(i);
               rank_send = mg_lev%level_ranks(i+1);
               rank_recv = mg_lev%level_ranks(i-1);
               mg_lev%send_procs(my_rank+1) = rank_send;
               mg_lev%recv_procs(my_rank+1) = rank_recv;
           end do
           rank_send = mg_lev%level_ranks(mg_lev%nprocs_level-1);
           my_rank = mg_lev%level_ranks(mg_lev%nprocs_level);
           mg_lev%recv_procs(my_rank+1) = rank_send;
        end do
     end if

     do level_index = nlevels,2,-1
        mg_lev => mg_ld(level_index)
        mg_c_lev => mg_ld(level_index-1)

        call FC_Setup(mg_c_lev%Nt, coarsen_factor, mg_lev%f_blocks, mg_lev%c_pts)
     end do

     do level_index = nlevels,1,-1
        mg_lev => mg_ld(level_index)
        pf_lev = pf%levels(level_index)
        if ((level_index .lt. nlevels) .and. (level_index .gt. 1) .and. (mg_lev%Nt .gt. 0)) then
           call pf_lev%ulevel%factory%create_array(mg_lev%g, mg_lev%Nt, level_index, pf_lev%lev_shape)
           do i = 1,mg_lev%Nt
              call mg_lev%g(i)%setval(0.0_pfdp)
           end do
        end if
        if (FAS_flag .eqv. .true.) then
           mg_ld(level_index)%FAS_flag = .true.
           if ((level_index .lt. nlevels) .and. (mg_lev%Nt .gt. 0)) then
              call pf_lev%ulevel%factory%create_array(mg_lev%qc_fas, mg_lev%Nt, level_index, pf_lev%lev_shape)
              do i = 1,mg_lev%Nt
                 call mg_lev%qc_fas(i)%setval(0.0_pfdp)
              end do
           end if
        end if
        if (FCF_flag .eqv. .false.) then
           mg_ld(level_index)%FCF_flag = .false.
        end if
        if (level_index .gt. 1) then
           mg_c_lev => mg_ld(level_index-1)
           if (mg_c_lev%Nt .gt. 0) then
              call pf_lev%ulevel%factory%create_array(mg_lev%qc, mg_c_lev%Nt, level_index, pf_lev%lev_shape)
              call pf_lev%ulevel%factory%create_array(mg_lev%qc_prev, mg_c_lev%Nt, level_index, pf_lev%lev_shape)
              do i = 1,mg_c_lev%Nt
                 call mg_lev%qc(i)%setval(0.0_pfdp)
                 call mg_lev%qc_prev(i)%setval(0.0_pfdp)
              end do
           end if
        end if
        call pf_lev%ulevel%factory%create_single(mg_lev%qc_boundary, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_single(mg_lev%q_temp, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_single(mg_lev%r, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_single(mg_lev%Q0, level_index, pf_lev%lev_shape)
     end do
  end subroutine mgrit_initialize

  subroutine FC_Setup(Nc, coarsen_factor, f_blocks, c_pts)
     type(int_vector), allocatable, intent(inout) :: f_blocks(:)
     integer, allocatable, intent(inout) :: c_pts(:)
     integer, intent(in) :: Nc, coarsen_factor
     integer, allocatable :: f_pts(:)
     integer :: i, j, n
     integer :: j_start, j_end, jj
     integer :: f_size

     n = coarsen_factor

     allocate(f_blocks(Nc))
     allocate(c_pts(Nc))

     j = n
     do i = 1,Nc
        c_pts(i) = j;
        j = j + n;
     end do
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
  end subroutine FC_Setup

  !>  Do the MGRIT algorithm
  subroutine pf_MGRIT_run(pf, mg_ld, Q0, qend)
    type(pf_pfasst_t), target, intent(inout)   :: pf   !!  The complete PFASST structure
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    class(pf_encap_t), intent(in   )           :: Q0   !!  The initial condition
    class(pf_encap_t), intent(inout), optional :: qend    !!  The computed solution at tend
    class(pf_encap_t), pointer :: qc(:)
    type(mgrit_level_data), pointer :: mg_lev
    type(pf_level_t), pointer :: pf_lev
    logical :: zero_rhs_flag

    !  Local variables
    integer :: nproc  !!  Total number of processors
    integer :: nsteps_loc  !!  local number of time steps
    real(pfdp) :: tend_loc !!  The final time of run
    integer :: ierror, iter, k, i, level_index
    integer :: Nt, Nt_glob
    integer :: qc_len
    real(pfdp) :: t0, dt, big_dt, t0_glob
    integer :: nlevels

    k = 1
    nlevels = pf%nlevels
    pf%state%pfblock = 1
    pf%state%nsteps = pf%comm%nproc
    pf%state%sweep = 1
    pf%state%pstatus = PF_STATUS_ITERATING

    !>  Allocate stuff for holding results 
    call initialize_results(pf)


    mg_ld(nlevels)%cycle_phase = 0
    level_index = 1
    call InitExactSolve(pf, mg_ld, Q0, level_index)
    zero_rhs_flag = .true.
    call IdealInterp(pf, mg_ld, 0, zero_rhs_flag)

    do level_index = 1,nlevels
       pf_lev => pf%levels(level_index)
       call pf_lev%qend%setval(0.0_pfdp)
       call pf_lev%q0%setval(0.0_pfdp)
    end do
    
    !>  Try to sync everyone
    call mpi_barrier(pf%comm%comm, ierror)

    !> Start timer
    if (pf%save_timings > 0) call pf_start_timer(pf, T_TOTAL)
    do iter = 1, pf%niters
       pf%state%iter = iter

       do level_index = 1,nlevels
          mg_lev => mg_ld(level_index)
          pf_lev => pf%levels(level_index)
          qc_len = size(mg_lev%qc)
          if (level_index .gt. 1) then
             call pf_lev%delta_q0%copy(mg_lev%qc(qc_len))
          else
             call pf_lev%delta_q0%copy(pf_lev%qend) 
          end if
       end do       

       !  Do a v_cycle
       if (pf%save_timings > 1) call pf_start_timer(pf, T_ITERATION)
       call pf_MGRIT_v_cycle(pf, mg_ld, iter)
       call pf_set_iter(pf, iter)
       if (pf%save_timings > 1) call pf_stop_timer(pf, T_ITERATION)

       do level_index = 1,nlevels
          mg_lev => mg_ld(level_index)
          pf_lev => pf%levels(level_index)
          qc_len = size(mg_lev%qc)
          if (level_index .gt. 1) then
             call pf_lev%delta_q0%axpy(-1.0_pfdp, mg_lev%qc(qc_len))
          else
             call pf_lev%delta_q0%axpy(-1.0_pfdp, pf_lev%qend)
          end if
          pf_lev%max_delta_q0 = pf_lev%delta_q0%norm()
          call pf_set_delta_q0(pf, level_index, pf_lev%delta_q0%norm())
          call pf_set_resid(pf, level_index, mg_lev%res_norm_loc(1))
       end do

       call call_hooks(pf, -1, PF_POST_ITERATION)
       if (pf%state%pstatus .eq. PF_STATUS_CONVERGED) then
          exit
       end if
    end do
    if (pf%save_timings > 0) call pf_stop_timer(pf, T_TOTAL)

    call pf_dump_stats(pf)

    qc => mg_ld(nlevels)%qc
    call qend%copy(qc(size(qc)))

   ! qc => mg_ld(pf%nlevels)%qc

   ! if (pf%rank .eq. 0) then
   !    do i = 1,size(qc)
   !       print *,i,qc(i)%norm()
   !    end do
   ! end if
   ! call mpi_barrier(pf%comm%comm, ierror)
   ! call mpi_barrier(pf%comm%comm, ierror)
   ! call mpi_barrier(pf%comm%comm, ierror)
   ! if (pf%rank .eq. 1) then
   !    do i = 1,size(qc)
   !       print *,i,qc(i)%norm()
   !    end do
   ! end if


   ! if (pf%rank .eq. pf%comm%nproc-1) then
   !    call qend%eprint()
   !
   !    t0 = mg_ld(pf%nlevels)%glob_t0
   !    Nt = mg_ld(pf%nlevels)%Nt * pf%comm%nproc
   !    dt = mg_ld(pf%nlevels)%dt
   !    big_dt = dt * real(Nt, pfdp)
   !    call pf%levels(pf%nlevels)%q0%copy(q0)
   !    call pf%levels(pf%nlevels)%ulevel%stepper%do_n_steps(pf, pf%nlevels, t0, pf%levels(pf%nlevels)%q0, &
   !                               pf%levels(pf%nlevels)%qend, big_dt, Nt)
   !    call qend%copy(pf%levels(pf%nlevels)%qend)
   ! end if
  end subroutine pf_MGRIT_run

  !> Execute a MGRIT V-cycle (iteration)
  !!  It is assumed that we have two levels and two nodes here
  !!  When this is called the previous coarse integrator result should be stored in Q(1)
  !!  and the MGRIT iteration in qend (both on coarse level).  If this is called
  !!  directly after the predictor, these will be the same thing
  subroutine pf_MGRIT_v_cycle(pf, mg_ld, iteration)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: iteration
    type(pf_level_t), pointer :: pf_lev, pf_f_lev, pf_c_lev
    type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev
    logical :: zero_rhs_flag, zero_c_pts_flag, interp_flag
    integer :: qc_len
    integer :: i, j, k, n, ierr
    integer :: nlevels, level_index, level_index_f
    integer :: ierror
    real(pfdp) :: res_norm_glob(1)
    class(pf_encap_t), pointer :: qc(:)
    real(pfdp), pointer :: qend(:)
    real(pfdp) :: qexact, maxerr

    nlevels = pf%nlevels;

    !> FCF-relaxation on finest grid
    mg_ld(nlevels)%cycle_phase = 1
    level_index = nlevels
    mg_lev => mg_ld(level_index)
    call FCF_Relax(pf, mg_ld, level_index, iteration)
    if (pf%state%pstatus .eq. PF_STATUS_CONVERGED) then
        return
    end if

    do level_index = nlevels-1,2,-1
        level_index_f = level_index+1
        mg_lev => mg_ld(level_index)
        mg_f_lev => mg_ld(level_index_f)

        !> Restrict residual
        call Restrict(pf, mg_ld, level_index, level_index_f)
        do i = 1,size(mg_f_lev%qc)
           call mg_f_lev%qc(i)%copy(mg_f_lev%qc_prev(i))
        end do
        if ((level_index_f .eq. nlevels) .and. (iteration .gt. 1)) then
           call mpi_allreduce(mg_f_lev%res_norm_loc(1), res_norm_glob(1), 1, myMPI_Datatype, MPI_MAX, pf%comm%comm, ierr)
           if (res_norm_glob(1) .lt. pf%abs_res_tol) then
              pf%state%pstatus = PF_STATUS_CONVERGED
              return
           end if
        end if

        !> FCF-relaxation on intermediate grids
        call FCF_Relax(pf, mg_ld, level_index, iteration)
    end do

    !> Coarsest grid solve
    mg_ld(nlevels)%cycle_phase = 2
    level_index = 1
    call ExactSolve(pf, mg_ld, level_index)
    if ((level_index+1 .eq. nlevels) .and. (iteration .gt. 1)) then
       mg_f_lev => mg_ld(nlevels)
       call mpi_allreduce(mg_f_lev%res_norm_loc(1), res_norm_glob(1), 1, myMPI_Datatype, MPI_MAX, pf%comm%comm, ierr)
       if (res_norm_glob(1) .lt. pf%abs_res_tol) then
          pf%state%pstatus = PF_STATUS_CONVERGED
          return
       end if
    end if
    zero_rhs_flag = .false.
    mg_ld(nlevels)%cycle_phase = 3
    call IdealInterp(pf, mg_ld, iteration, zero_rhs_flag)

    mg_lev => mg_ld(nlevels)
    mg_c_lev => mg_ld(nlevels-1)
  end subroutine pf_MGRIT_v_cycle

  subroutine FCF_Relax(pf, mg_ld, level_index, iteration)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: iteration
     integer, intent(in) :: level_index
     type(pf_level_t), pointer :: pf_lev
     type(mgrit_level_data), pointer :: mg_lev, mg_c_lev
     logical :: zero_rhs_flag, zero_c_pts_flag, interp_flag
     integer :: qc_len, i, nlevels
     integer :: ierror
     class(pf_encap_t), pointer :: qc(:)

     nlevels = pf%nlevels
     mg_lev => mg_ld(level_index)
     mg_c_lev => mg_ld(level_index-1)
     pf_lev => pf%levels(level_index)

     if (level_index .lt. nlevels) then
        zero_c_pts_flag = .true.
     else
        zero_c_pts_flag = .false. 
     end if

     if (level_index .lt. nlevels) then
        zero_rhs_flag = .false.
     else
        zero_rhs_flag = .true.
     end if
     
     interp_flag = .false.

     qc_len = size(mg_lev%qc)
     do i = 1,qc_len
        call mg_lev%qc_prev(i)%copy(mg_lev%qc(i))
     end do

     if (zero_c_pts_flag .eqv. .false.) then
        if (pf%rank .lt. pf%comm%nproc-1) then
           call pf_lev%qend%copy(mg_lev%qc(qc_len))
        end if
        call send_recv_C_points(pf, mg_ld, level_index)
     else
        call pf_lev%q0%setval(0.0_pfdp)
     end if
     call F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag)

     if (mg_lev%FCF_flag .eqv. .true.) then
        zero_c_pts_flag = .false.
        call C_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag)
        !> Send and receive boundary C-points
        if (pf%rank .lt. pf%comm%nproc-1) then
           call pf_lev%qend%copy(mg_lev%qc(qc_len))
        end if
        call send_recv_C_points(pf, mg_ld, level_index)
        if (pf%rank .gt. 0) then
           call mg_lev%qc_boundary%copy(pf_lev%q0)
        end if
        do i = 1,qc_len
           call mg_lev%qc_prev(i)%copy(mg_lev%qc(i))
        end do
        call F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag)
     end if

  end subroutine FCF_Relax

  subroutine IdealInterp(pf, mg_ld, iteration, zero_rhs_flag)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: iteration
    logical, intent(in) :: zero_rhs_flag
    type(pf_level_t), pointer :: pf_lev, pf_f_lev, pf_c_lev
    type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev
    logical :: zero_c_pts_flag, interp_flag
    integer :: qc_len
    integer :: i, j, k, n, ierr
    integer :: nlevels, level_index
    class(pf_encap_t), pointer :: qc(:)

    nlevels = pf%nlevels;
    zero_c_pts_flag = .false.
    interp_flag = .true.

    do level_index = 2,(nlevels-1)
       mg_lev => mg_ld(level_index)
       pf_lev => pf%levels(level_index)
       qc_len = size(mg_lev%qc)
       do i = 1,qc_len
          call mg_lev%qc_prev(i)%copy(mg_lev%qc(i))
       end do
       if (pf%comm%nproc .gt. 1) then
          if (pf%rank .lt. pf%comm%nproc-1) then
             call pf_lev%qend%copy(mg_lev%qc(qc_len))
          end if
          call send_recv_C_points(pf, mg_ld, level_index)
       end if
       call F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag)
       do i = 1,qc_len
          call mg_lev%qc(i)%copy(mg_lev%qc_prev(i))
       end do
    end do
  end subroutine IdealInterp

  subroutine F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     logical, intent(in) :: zero_rhs_flag, zero_c_pts_flag, interp_flag
     type(mgrit_level_data), pointer :: mg_f_lev, mg_lev
     type(pf_level_t), pointer :: pf_lev, pf_f_lev
     integer :: i, j, n, ii
     integer :: nlevels
     logical :: q_zero_flag
     integer :: level_index_f

     nlevels = pf%nlevels
     level_index_f = level_index + 1;
     mg_lev => mg_ld(level_index)
     pf_lev => pf%levels(level_index)
     if (interp_flag .eqv. .true.) then
        mg_f_lev => mg_ld(level_index_f)
        pf_f_lev => pf%levels(level_index_f)
     end if

     do i = 1, size(mg_lev%f_blocks)
        if ((i .eq. 1) .and. (pf%rank .eq. 0)) then
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
           else if (i .gt. 1) then
              call pf_lev%q0%copy(mg_lev%qc_prev(i-1))
           end if
        end if
        
        do n = 1,size(mg_lev%f_blocks(i)%val)
           j = mg_lev%f_blocks(i)%val(n)
           !> Do a single step
           call PointRelax(pf, mg_ld, level_index, j, pf_lev%q0, pf_lev%qend)
           !> Add g to the result
           if (zero_rhs_flag .eqv. .false.) then
              call pf_lev%qend%axpy(1.0_pfdp, mg_lev%g(j))
           end if

           !> Interpolate to fine level
           if (interp_flag .eqv. .true.) then
              call mg_lev%q_temp%copy(pf_lev%qend)
              if ((mg_ld(nlevels)%FAS_flag .eqv. .true.) .and. (mg_ld(nlevels)%cycle_phase .gt. 0)) then
                 call mg_lev%q_temp%axpy(-1.0_pfdp, mg_lev%qc_fas(j))
              end if
              call pf_f_lev%ulevel%interpolate(pf_f_lev, pf_lev, mg_f_lev%q_temp, mg_lev%q_temp, mg_f_lev%t0)
              call mg_f_lev%qc(j)%axpy(1.0_pfdp, mg_f_lev%q_temp)
           end if
           call pf_lev%q0%copy(pf_lev%qend)
        end do
        call mg_lev%qc(i)%copy(pf_lev%qend)
        if (interp_flag .eqv. .true.) then
           call mg_lev%q_temp%copy(mg_lev%qc_prev(i))
           if ((mg_ld(nlevels)%FAS_flag .eqv. .true.) .and. (mg_ld(nlevels)%cycle_phase .gt. 0)) then
              call mg_lev%q_temp%axpy(-1.0_pfdp, mg_lev%qc_fas(j+1))
           end if
           call pf_f_lev%ulevel%restrict(pf_f_lev, pf_lev, mg_f_lev%q_temp, mg_lev%q_temp, mg_lev%t0)
           call mg_f_lev%qc(j+1)%axpy(1.0_pfdp, mg_f_lev%q_temp)
        end if
     end do
  end subroutine F_Relax

  subroutine C_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     logical, intent(in) :: zero_rhs_flag, interp_flag
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: q_zero_flag
     type(mgrit_level_data), pointer :: mg_lev, mg_finest_lev
     type(pf_level_t), pointer :: pf_lev

     nlevels = pf%nlevels

     mg_lev => mg_ld(level_index)

     pf_lev => pf%levels(level_index)

     do i = 1,size(mg_lev%c_pts)
        j = mg_lev%c_pts(i)
        call pf_lev%q0%copy(mg_lev%qc(i))
        call PointRelax(pf, mg_ld, level_index, j, pf_lev%q0, pf_lev%qend)
        call mg_lev%qc(i)%copy(pf_lev%qend)
        if (zero_rhs_flag .eqv. .false.) then
           call mg_lev%qc(i)%axpy(1.0_pfdp, mg_lev%g(j))
        end if
    end do
  end subroutine C_Relax

  subroutine send_recv_C_points(pf, mg_ld, level_index)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: level_index
    type(pf_level_t), pointer :: pf_lev
    integer :: ierror

    if (pf%comm%nproc .gt. 1) then
       pf_lev => pf%levels(level_index)
       if (mod(pf%rank, 2) .eq. 0) then
          if (pf%rank .gt. 0) then
             !call pf_recv(pf, pf_lev, 1, .true.)
             call mgrit_recv(pf, mg_ld, level_index, 1, .true.)
          end if
          if (pf%rank .lt. pf%comm%nproc-1) then
             !call pf_send(pf, pf_lev, 1, .true.)
             call mgrit_send(pf, mg_ld, level_index, 1, .true.)
          end if
       else
          if (pf%rank .lt. pf%comm%nproc-1) then
             !call pf_send(pf, pf_lev, 1, .true.)
             call mgrit_send(pf, mg_ld, level_index, 1, .true.)
          end if
          if (pf%rank .gt. 0) then
             !call pf_recv(pf, pf_lev, 1, .true.)
             call mgrit_recv(pf, mg_ld, level_index, 1, .true.)
          end if
       end if
    end if
  end subroutine send_recv_C_points

  subroutine mgrit_send(pf, mg_ld, level_index, tag, blocking)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: level_index, tag
    logical, intent(in) :: blocking
    type(pf_level_t), pointer :: pf_lev
    integer :: ierror

    pf_lev => pf%levels(level_index)
    call pf_lev%qend%pack(pf_lev%send)
    call pf%comm%send(pf, pf_lev, tag, blocking, ierror, mg_ld(level_index)%send_procs(pf%rank+1))
  end subroutine mgrit_send

  subroutine mgrit_recv(pf, mg_ld, level_index, tag, blocking)
    type(pf_pfasst_t), target, intent(inout) :: pf
    type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
    integer, intent(in) :: level_index, tag
    logical, intent(in) :: blocking
    type(pf_level_t), pointer :: pf_lev
    integer :: ierror

    pf_lev => pf%levels(level_index)
    call pf%comm%recv(pf, pf_lev, tag, blocking, ierror, mg_ld(level_index)%recv_procs(pf%rank+1))
    call pf_lev%q0%unpack(pf_lev%recv)
  end subroutine mgrit_recv

  subroutine ExactSolve(pf, mg_ld, level_index)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     integer :: level_index_f
     type(pf_level_t), pointer :: lev
     class(pf_encap_t), allocatable :: gi
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: zero_rhs_flag
     type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_finest_lev
     type(pf_level_t), pointer :: pf_lev, pf_f_lev

     nlevels = pf%nlevels
     level_index_f = level_index + 1

     mg_lev => mg_ld(level_index)
     mg_f_lev => mg_ld(level_index_f)
     mg_finest_lev => mg_ld(nlevels)

     pf_lev => pf%levels(level_index)
     pf_f_lev => pf%levels(level_index_f)

     mg_f_lev%res_norm_loc(1) = 0.0_pfdp

     if ((pf%rank .gt. 0) .and. (pf%comm%nproc .gt. 1)) then
        call pf_recv(pf, pf_lev, 2, .true.) 
     else
        if (mg_ld(nlevels)%FAS_flag .eqv. .true.) then
           call pf_lev%q0%copy(mg_lev%Q0)
        else
           call pf_lev%q0%setval(0.0_pfdp);
        end if
     end if

     call pf_lev%ulevel%factory%create_single(gi, level_index, pf_lev%lev_shape)
     do i = 1,mg_lev%Nt
        ii = mg_f_lev%c_pts(i)
        if (level_index_f .eq. nlevels) then
            zero_rhs_flag = .true.
        else
            zero_rhs_flag = .false.
        end if
        call InjectRestrictPoint(pf, mg_ld, gi, level_index, level_index_f, i, ii, zero_rhs_flag)
        call mg_f_lev%qc(i)%copy(mg_f_lev%qc_prev(i))

        call PointRelax(pf, mg_ld, level_index, i, pf_lev%q0, pf_lev%qend)
        call pf_lev%qend%axpy(1.0_pfdp, gi)
        call pf_lev%q0%copy(pf_lev%qend)

        call mg_lev%q_temp%copy(pf_lev%qend)
        if ((mg_ld(nlevels)%FAS_flag .eqv. .true.) .and. (mg_ld(nlevels)%cycle_phase .gt. 0)) then
           call mg_lev%q_temp%axpy(-1.0_pfdp, mg_lev%qc_fas(i))
        end if
        call pf_f_lev%ulevel%interpolate(pf_f_lev, pf_lev, mg_f_lev%q_temp, mg_lev%q_temp, mg_f_lev%t0)
        call mg_f_lev%qc(i)%axpy(1.0_pfdp, mg_f_lev%q_temp)
     end do
     call pf_lev%ulevel%factory%destroy_single(gi)

     pf_f_lev%residual = mg_f_lev%res_norm_loc(1)
     pf_lev%residual = 0.0_pfdp

     if ((pf%rank .lt. pf%comm%nproc-1) .and. (pf%comm%nproc .gt. 1)) then
        call pf_send(pf, pf_lev, 2, .true.)
     end if
  end subroutine ExactSolve
 
  subroutine InitExactSolve(pf, mg_ld, Q0, level_index)
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

     nlevels = pf%nlevels
     level_index_f = level_index + 1

     mg_lev => mg_ld(level_index)
     mg_f_lev => mg_ld(level_index_f)
     mg_finest_lev => mg_ld(nlevels)

     pf_lev => pf%levels(level_index)
     pf_f_lev => pf%levels(level_index_f)

     if ((pf%rank .gt. 0) .and. (pf%comm%nproc .gt. 1)) then
        call pf_recv(pf, pf_lev, 3, .true.)
     else
        call mg_finest_lev%Q0%copy(Q0)
        do l = nlevels,2,-1
           call pf%levels(l)%ulevel%restrict(pf%levels(l), pf%levels(l-1), mg_ld(l)%Q0, mg_ld(l-1)%Q0, mg_ld(l)%t0)
        end do
        call pf_lev%q0%copy(mg_lev%Q0)
     end if

     do i = 1,mg_lev%Nt
        call PointRelax(pf, mg_ld, level_index, i, pf_lev%q0, pf_lev%qend)
        call pf_lev%ulevel%interpolate(pf_f_lev, pf_lev, mg_f_lev%qc(i), pf_lev%qend, mg_f_lev%t0)
        call pf_lev%q0%copy(pf_lev%qend)
     end do

     if ((pf%rank .lt. pf%comm%nproc-1) .and. (pf%comm%nproc .gt. 1)) then
        call pf_send(pf, pf_lev, 3, .true.)
     end if
  end subroutine InitExactSolve

  subroutine Restrict(pf, mg_ld, level_index_c, level_index_f)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index_c, level_index_f
     class(pf_encap_t), pointer :: gi
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: zero_rhs_flag
     real(pfdp) :: r_norm
     type(mgrit_level_data), pointer :: mg_c_lev, mg_f_lev
     type(pf_level_t), pointer :: pf_lev, pf_f_lev

     nlevels = pf%nlevels

     mg_c_lev => mg_ld(level_index_c)
     mg_f_lev => mg_ld(level_index_f)

     pf_f_lev => pf%levels(level_index_f)

     mg_f_lev%res_norm_loc(1) = 0.0_pfdp
     mg_f_lev%res_norm_index = 0

     do i = 1,mg_c_lev%Nt
        ii = mg_f_lev%c_pts(i)
        if (level_index_f .eq. nlevels) then
           zero_rhs_flag = .true.;
        else
           zero_rhs_flag = .false.;
        end if
        call InjectRestrictPoint(pf, mg_ld, mg_c_lev%g(i), level_index_c, level_index_f, i, ii, zero_rhs_flag);
     end do

     pf_f_lev%residual = mg_f_lev%res_norm_loc(1)
  end subroutine Restrict

  subroutine InjectRestrictPoint(pf, mg_ld, gci, level_index_c, level_index_f, i_c, i_f, zero_rhs_flag)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     class(pf_encap_t), intent(inout) :: gci
     integer, intent(in) :: level_index_c, level_index_f, i_c, i_f
     logical, intent(in) :: zero_rhs_flag
     type(mgrit_level_data), pointer :: mg_c_lev, mg_f_lev
     type(pf_level_t), pointer :: pf_f_lev, pf_c_lev
     integer :: nlevels

     mg_c_lev => mg_ld(level_index_c)
     mg_f_lev => mg_ld(level_index_f)

     pf_f_lev => pf%levels(level_index_f)
     pf_c_lev => pf%levels(level_index_c)
     nlevels = pf%nlevels

     call mg_f_lev%q_temp%copy(mg_f_lev%qc(i_c))
     call PointRelax(pf, mg_ld, level_index_f, i_f, mg_f_lev%q_temp, mg_f_lev%r)
     call mg_f_lev%r%axpy(-1.0_pfdp, mg_f_lev%qc_prev(i_c))
     if (zero_rhs_flag .eqv. .false.) then
         call mg_f_lev%r%axpy(1.0_pfdp, mg_f_lev%g(i_f));
     end if
     if (i_c .eq. mg_c_lev%Nt) then
        call ResNorm(pf, mg_ld, level_index_f, mg_f_lev%r, i_c)
     end if
     call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%r, mg_c_lev%r, mg_f_lev%t0)
     if (mg_ld(nlevels)%FAS_flag .eqv. .true.) then
        call pf_f_lev%ulevel%restrict(pf_f_lev, pf_c_lev, mg_f_lev%qc_prev(i_c), mg_c_lev%qc_fas(i_c), mg_f_lev%t0)
        call mg_c_lev%r%axpy(1.0_pfdp, mg_c_lev%qc_fas(i_c))
        if (i_c .eq. 1) then
           if (pf%rank .eq. 0) then
              call mg_c_lev%q_temp%copy(mg_c_lev%Q0)
           else
              call mg_c_lev%q_temp%copy(mg_f_lev%qc_boundary)
           end if
        else
           call mg_c_lev%q_temp%copy(mg_c_lev%qc_fas(i_c-1))
        end if
        call PointRelax(pf, mg_ld, level_index_c, i_c, mg_c_lev%q_temp, pf_c_lev%qend)
        call mg_c_lev%r%axpy(-1.0_pfdp, pf_c_lev%qend)
     end if
     call gci%copy(mg_c_lev%r)
  end subroutine InjectRestrictPoint

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

  subroutine ResNorm(pf, mg_ld, level_index, r, i)
     type(pf_pfasst_t), target, intent(inout) :: pf
     type(mgrit_level_data), allocatable, target, intent(inout) :: mg_ld(:)
     class(pf_encap_t), intent(in) :: r
     integer, intent(in) :: level_index, i
     type(mgrit_level_data), pointer :: mg_lev
     type(pf_level_t), pointer :: pf_lev
     integer :: nlevels
     real(pfdp) :: r_norm

     nlevels = pf%nlevels

     mg_lev => mg_ld(level_index)
     pf_lev => pf%levels(level_index)
     r_norm = r%norm()
     if (r_norm .gt. mg_lev%res_norm_loc(1)) then
        mg_lev%res_norm_loc(1) = r_norm
        mg_lev%res_norm_index = i
     end if
  end subroutine ResNorm

end module pf_mod_MGRIT
