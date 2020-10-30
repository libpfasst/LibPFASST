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
  
  implicit none

  type :: int_vector
     integer, allocatable :: val(:)
  end type int_vector

  type :: mgrit_level_data
     integer :: Nt
     real(pfdp) :: dt
     real(pfdp) :: T0
     real(pfdp) :: Tfin
     type(int_vector), allocatable :: f_blocks(:)
     integer, allocatable :: c_pts(:)
     class(pf_encap_t), allocatable :: g(:)
     class(pf_encap_t), allocatable :: qc(:)
     class(pf_encap_t), allocatable :: qc_prev(:)
     integer, allocatable :: interp_map(:)
  end type 
  
contains

  subroutine mgrit_initialize(pf, mg_ld, T0, Tfin, n_coarse, refine_factor)
     type(pf_pfasst_t), intent(inout) :: pf
     type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
     integer, intent(in) :: n_coarse, refine_factor
     real(pfdp) :: T0, Tfin
     integer :: level_index, nlevels, N, i, kk
     type(pf_level_t) :: pf_lev, pf_f_lev, pf_c_lev
     type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev

     nlevels = pf%nlevels
     n = refine_factor

     allocate(mg_ld(nlevels))

     N = n_coarse;
     do level_index = 1,nlevels
        mg_lev => mg_ld(level_index)

        mg_lev%T0 = T0
        mg_lev%Tfin = Tfin
        mg_lev%Nt = N
        mg_lev%dt = (Tfin-T0)/(mg_lev%Nt)
        N = N * refine_factor
     end do

     do level_index = nlevels,2,-1
        mg_lev => mg_ld(level_index)
        mg_c_lev => mg_ld(level_index-1)

        call FC_Setup(mg_c_lev%Nt, refine_factor, mg_lev%f_blocks, mg_lev%c_pts)
     end do

     kk = 0;
     do level_index = (nlevels-1),1,-1
        mg_lev => mg_ld(level_index)

        allocate(mg_lev%interp_map(mg_lev%Nt))
        do i = 1,mg_lev%Nt
            mg_lev%interp_map(i) = i*refine_factor**kk
        end do
        kk = kk + 1;
     end do
     do level_index = nlevels,2,-1
        mg_lev => mg_ld(level_index)
        mg_c_lev => mg_ld(level_index-1)

        pf_lev = pf%levels(level_index)
        if (level_index < nlevels .and. level_index > 1) then
           call pf_lev%ulevel%factory%create_array(mg_lev%g, mg_lev%Nt, level_index, pf_lev%lev_shape)
           do i = 1,mg_lev%Nt
               call mg_lev%g(i)%setval(0.0_pfdp)
           end do
        end if
        call pf_lev%ulevel%factory%create_array(mg_lev%qc, mg_c_lev%Nt, level_index, pf_lev%lev_shape)
        call pf_lev%ulevel%factory%create_array(mg_lev%qc_prev, mg_c_lev%Nt, level_index, pf_lev%lev_shape)
        do i = 1,mg_c_lev%Nt
           call mg_lev%qc(i)%setval(0.0_pfdp)
           call mg_lev%qc_prev(i)%setval(0.0_pfdp)
        end do
     end do
  end subroutine mgrit_initialize

  subroutine FC_Setup(Nc, refine_factor, f_blocks, c_pts)
     type(int_vector), allocatable, intent(inout) :: f_blocks(:)
     integer, allocatable, intent(inout) :: c_pts(:)
     integer, intent(in) :: Nc, refine_factor
     integer, allocatable :: f_pts(:)
     integer :: i, j, n
     integer :: j_start, j_end, jj
     integer :: f_size

     n = refine_factor

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
  subroutine pf_MGRIT_run(pf, mg_ld, q0, qend)
    type(pf_pfasst_t), intent(inout), pointer   :: pf   !!  The complete PFASST structure
    type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
    class(pf_encap_t), intent(in   )           :: q0   !!  The initial condition
    class(pf_encap_t), intent(inout), optional :: qend    !!  The computed solution at tend
    class(pf_encap_t), allocatable :: qc(:)

    !  Local variables
    integer :: nproc  !!  Total number of processors
    integer :: nsteps_loc  !!  local number of time steps
    real(pfdp) :: tend_loc !!  The final time of run
    integer :: ierr, iter, k

    k = 1


    !>  Allocate stuff for holding results 
    call initialize_results(pf)
    
    !>  Try to sync everyone
    call mpi_barrier(pf%comm%comm, ierr)

    !> Start timer
    if (pf%save_timings > 0) call pf_start_timer(pf, T_TOTAL)

    do iter = 1, pf%niters
       call call_hooks(pf, -1, PF_PRE_ITERATION)
       if (pf%save_timings > 1) call pf_start_timer(pf, T_ITERATION)

       pf%state%iter = iter

       !  Do a v_cycle
       call pf_MGRIT_v_cycle(pf, mg_ld, q0, iter)

       !  Check for convergence
       call pf_check_convergence_block(pf, pf%state%finest_level, send_tag=1111*k+iter)

       if (pf%save_timings > 1) call pf_stop_timer(pf, T_ITERATION)
       call call_hooks(pf, -1, PF_POST_ITERATION)

       !  If we are converged, exit block
       if (pf%state%status == PF_STATUS_CONVERGED)  then
          call call_hooks(pf, -1, PF_POST_CONVERGENCE)
          call pf_set_iter(pf,iter)
          exit
       end if
    end do  !  Loop over j, the iterations in this block

    qc = mg_ld(pf%nlevels)%qc
    call qend%copy(qc(size(qc)))

    !> End timer    
    if (pf%save_timings > 0) call pf_stop_timer(pf, T_TOTAL)

    !>  Output stats
    call pf_dump_stats(pf)
  end subroutine pf_MGRIT_run

  !> Execute a MGRIT V-cycle (iteration)
  !!  It is assumed that we have two levels and two nodes here
  !!  When this is called the previous coarse integrator result should be stored in Q(1)
  !!  and the MGRIT iteration in qend (both on coarse level).  If this is called
  !!  directly after the predictor, these will be the same thing
  subroutine pf_MGRIT_v_cycle(pf, mg_ld, Q0, iteration)
    type(pf_pfasst_t), pointer, intent(inout) :: pf
    type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
    class(pf_encap_t), intent(in) :: Q0
    integer, intent(in) :: iteration
    type(pf_level_t) :: pf_lev, pf_f_lev, pf_c_lev
    type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_c_lev
    logical :: zero_rhs_flag, zero_c_pts_flag, interp_flag
    integer :: qc_len
    integer :: i, j, k, n
    integer :: nlevels, level_index, level_index_f

    nlevels = pf%nlevels;

    !> FCF-relaxation on finest grid
    level_index = nlevels
    mg_lev => mg_ld(level_index)

    zero_rhs_flag = .true.
    interp_flag = .false.
    zero_c_pts_flag = .false.
    qc_len = size(mg_lev%qc)
    do i = 1,qc_len
       call mg_lev%qc_prev(i)%copy(mg_lev%qc(i))
    end do
    call F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
    call C_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag)
    do i = 1,qc_len
       call mg_lev%qc_prev(i)%copy(mg_lev%qc(i))
    end do
    call F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
    zero_rhs_flag = .false.
    do level_index = nlevels-1,2,-1
        level_index_f = level_index+1
        mg_lev => mg_ld(level_index)
        mg_f_lev => mg_ld(level_index_f)

        !> Restrict residual
        call Restrict(pf, mg_ld, level_index, level_index_f)
        qc_len = size(mg_f_lev%qc)
        do i = 1,qc_len
           call mg_f_lev%qc(i)%copy(mg_f_lev%qc_prev(i))
        end do

        !> FCF-relaxation on intermediate grids
        interp_flag = .false.
        zero_c_pts_flag = .true.
        qc_len = size(mg_lev%qc)
        do i = 1,qc_len
           call mg_lev%qc_prev(i)%copy(mg_lev%qc(i))
        end do
        call F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
        interp_flag = .true.
        zero_c_pts_flag = .false.
        call C_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag)
        do i = 1,qc_len
           call mg_lev%qc_prev(i)%copy(mg_lev%qc(i))
        end do
        call F_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
    end do
    !> Coarsest grid solve
    level_index = nlevels
    call ExactSolve(pf, mg_ld, level_index)

  end subroutine pf_MGRIT_v_cycle

  subroutine F_Relax(pf, mg_ld, level_index, interp_flag, zero_rhs_flag, zero_c_pts_flag, Q0)
     type(pf_pfasst_t), pointer, intent(inout) :: pf
     type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     logical, intent(in) :: zero_rhs_flag, zero_c_pts_flag, interp_flag
     class(pf_encap_t), intent(in) :: Q0
     type(mgrit_level_data), pointer :: mg_finest_lev, mg_lev
     type(pf_level_t), pointer :: pf_lev
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: q_zero_flag
     real(pfdp) :: t0, dt

     nlevels = pf%nlevels

     mg_lev => mg_ld(level_index)
     mg_finest_lev => mg_ld(nlevels)
     
     pf_lev => pf%levels(level_index)

     q_zero_flag = .false.
     do i = 1, size(mg_lev%f_blocks)
        !> Check if inital values of C-points are zero
        if (zero_c_pts_flag .eqv. .true.) then
           call pf_lev%q0%setval(0.0_pfdp)
           q_zero_flag = .true.
        else
           if (i .eq. 1) then
               !> If we are at the first F-point on the fine level, use ODE initial value.
               !> Otherwise, set initial value to zero.
               if (k .eq. nlevels) then
                   call pf_lev%q0%copy(Q0);
               else
                   call pf_lev%q0%setval(0.0_pfdp)
                   q_zero_flag = .true.
               end if
           else
               call pf_lev%q0%copy(mg_lev%qc_prev(i-1))
           end if
        end if
        do n = 1,size(mg_lev%f_blocks(i)%val)
           j = mg_lev%f_blocks(i)%val(n)
           if (q_zero_flag .eqv. .true.) then
              !> Ifi intial C-points are zero, just copy RHS into q
              if (zero_rhs_flag .eqv. .false.) then
                 call pf_lev%qend%axpy(1.0_pfdp, mg_lev%g(j))
              end if
              q_zero_flag = .false.;
           else 
              !> Do a single step
              call Point_Relax(pf, mg_ld, level_index, j, pf_lev%q0, pf_lev%qend)
              !> Add g to the result
              if (zero_rhs_flag .eqv. .false.) then
                 call pf_lev%qend%axpy(1.0_pfdp, mg_lev%g(j))
              end if
           end if

           !> Interpolate to finest level
           if (interp_flag .eqv. .true.) then
              ii = mg_lev%interp_map(j)
              call mg_finest_lev%qc(ii)%axpy(1.0_pfdp, pf_lev%qend)
           end if
           call pf_lev%q0%copy(pf_lev%qend)
        end do
        !> Store result (used for C-relaxation)
        call mg_lev%qc(i)%copy(pf_lev%qend)
     end do
  end subroutine F_Relax

  subroutine C_Relax(pf, mg_ld, level_index, zero_rhs_flag, interp_flag)
     type(pf_pfasst_t), pointer, intent(inout) :: pf
     type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     logical, intent(in) :: zero_rhs_flag, interp_flag
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: q_zero_flag
     real(pfdp) :: t0, dt
     type(mgrit_level_data), pointer :: mg_lev, mg_finest_lev
     type(pf_level_t), pointer :: pf_lev

     nlevels = pf%nlevels

     mg_lev => mg_ld(level_index)
     mg_finest_lev => mg_ld(nlevels)

     pf_lev => pf%levels(level_index)

     do i = 1,size(mg_lev%c_pts)
        j = mg_lev%c_pts(i)
        call pf_lev%q0%copy(mg_lev%qc(i))
        call Point_Relax(pf, mg_ld, level_index, j, pf_lev%q0, pf_lev%qend)
        call mg_lev%qc(i)%copy(pf_lev%qend)
        if (zero_rhs_flag .eqv. .false.) then
           call mg_lev%qc(i)%axpy(1.0_pfdp, mg_lev%g(j))
        end if
        if (interp_flag .eqv. .true.) then
           ii = mg_lev%interp_map(j)
           call mg_finest_lev%qc(ii)%axpy(1.0_pfdp, mg_lev%qc(i))
        end if
    end do
  end subroutine C_Relax

  subroutine ExactSolve(pf, mg_ld, level_index)
     type(pf_pfasst_t), pointer, intent(inout) :: pf
     type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index
     integer :: level_index_f
     type(pf_level_t), pointer :: lev
     class(pf_encap_t), allocatable :: gi
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: zero_rhs_flag
     type(mgrit_level_data), pointer :: mg_lev, mg_f_lev, mg_finest_lev
     type(pf_level_t), pointer :: pf_lev

     nlevels = pf%nlevels
     level_index_f = level_index + 1

     mg_lev => mg_ld(level_index)
     mg_finest_lev => mg_ld(nlevels)
     mg_f_lev => mg_ld(level_index_f)

     pf_lev => pf%levels(level_index)

    ! if (pf%rank .gt. 0) then
    !    call pf_recv(pf, lev, 10000+iteration, .true.)        
    ! end if
     call pf_lev%ulevel%factory%create_single(gi, level_index, pf_lev%lev_shape)
     do i = 1,mg_lev%Nt
         ii = mg_f_lev%c_pts(i)
         if (level_index_f .eq. nlevels) then
             zero_rhs_flag = .true.;
         else
             zero_rhs_flag = .false.;
         end if
         call InjectRestrictPoint(pf, mg_ld, gi, level_index, level_index_f, i, ii, zero_rhs_flag);
         call mg_f_lev%qc(i)%copy(mg_f_lev%qc_prev(i));
         if (i .eq. 1) then
            call pf_lev%qend%copy(gi)   
         else
            call Point_Relax(pf, mg_ld, level_index, i, pf_lev%q0, pf_lev%qend)
            call pf_lev%qend%axpy(1.0_pfdp, gi)
         end if
         ii = mg_lev%interp_map(i)
         call mg_finest_lev%qc(ii)%axpy(1.0_pfdp, pf_lev%qend)
     end do
     call pf_lev%ulevel%factory%destroy_single(gi)
    ! if (pf%rank .lt. pf%nprocs-1) then
    !    call pf_send(pf, lev, 10000+iteration, .false.)
    ! end if
  end subroutine ExactSolve

  subroutine Restrict(pf, mg_ld, level_index_c, level_index_f)
     type(pf_pfasst_t), pointer, intent(inout) :: pf
     type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index_c, level_index_f
     class(pf_encap_t), pointer :: gi
     integer :: i, j, k, n, ii
     integer :: nlevels
     logical :: zero_rhs_flag
     type(mgrit_level_data), pointer :: mg_c_lev, mg_f_lev

     mg_c_lev => mg_ld(level_index_c)
     mg_f_lev => mg_ld(level_index_f)

     do i = 1,mg_c_lev%Nt
         ii = mg_f_lev%c_pts(i)
         if (level_index_f .eq. nlevels) then
            zero_rhs_flag = .true.;
         else
            zero_rhs_flag = .false.;
         end if
         call InjectRestrictPoint(pf, mg_ld, mg_c_lev%g(i), level_index_c, level_index_f, i, ii, zero_rhs_flag);
     end do
  end subroutine Restrict

  subroutine InjectRestrictPoint(pf, mg_ld, g_c, level_index_c, level_index_f, i_c, i_f, zero_rhs_flag)
     type(pf_pfasst_t), pointer, intent(inout) :: pf
     type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
     class(pf_encap_t), intent(inout) :: g_c
     integer, intent(in) :: level_index_c, level_index_f, i_c, i_f
     logical, intent(in) :: zero_rhs_flag

     type(mgrit_level_data), pointer :: mg_c_lev, mg_f_lev
     type(pf_level_t), pointer :: pf_f_lev

     mg_c_lev => mg_ld(level_index_c)
     mg_f_lev => mg_ld(level_index_f)

     pf_f_lev => pf%levels(level_index_f)

     call pf_f_lev%q0%copy(mg_f_lev%qc(i_c))
     call Point_Relax(pf, mg_ld, level_index_f, i_f, pf_f_lev%q0, pf_f_lev%qend)
     call g_c%copy(pf_f_lev%qend)
     call g_c%axpy(-1.0_pfdp, mg_f_lev%qc_prev(i_c))
     if (zero_rhs_flag .eqv. .false.) then
         call g_c%axpy(1.0_pfdp, mg_f_lev%g(i_f));
     end if
  end subroutine InjectRestrictPoint

  subroutine Point_Relax(pf, mg_ld, level_index, n, q0, qend)
     type(pf_pfasst_t), pointer, intent(inout) :: pf
     type(mgrit_level_data), pointer, intent(inout) :: mg_ld(:)
     integer, intent(in) :: level_index, n
     class(pf_encap_t), intent(inout) :: q0, qend
     type(pf_level_t), pointer :: pf_lev
     type(mgrit_level_data), pointer :: mg_lev
     real(pfdp) :: t0n

     mg_lev => mg_ld(level_index)
     pf_lev => pf%levels(level_index)

     t0n = mg_lev%t0 + mg_lev%dt * real(n-1,pfdp)
     call pf_lev%ulevel%stepper%do_n_steps(pf, level_index, t0n, q0, qend, mg_lev%dt, 1)
  end subroutine Point_Relax


  !> Subroutine to check if the current processor has converged and
  !> to update the next processor on the status
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor can't be converged yet either
  subroutine pf_check_convergence_block(pf, level_index, send_tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    integer,           intent(in)    :: send_tag  !! identifier for status send and receive
    
    logical           :: residual_converged, converged


    ! Shortcut for fixed iteration mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)

    !> Check to see if tolerances are met
    call pf_check_residual(pf, level_index, residual_converged)

    !>  Until I hear the previous processor is done, recieve it's status
    if (pf%state%pstatus /= PF_STATUS_CONVERGED) call pf_recv_status(pf, send_tag)

    !>  Check to see if I am converged
    converged = .false.
    if (residual_converged) then
       if (pf%rank == 0) then
          converged = .true.
       else  !  I am not the first processor, so I need to check the previous one
          if (pf%state%pstatus == PF_STATUS_CONVERGED) converged = .true.
       end if
    end if ! (residual_converged)

    !>  For MGRIT, Proc N is converged after iteration N
    if (pf%rank .lt. pf%state%iter) then
       converged = .true.
    end if
    
    !> Assign status and send it forward
    if (converged) then
       if (pf%state%status == PF_STATUS_ITERATING) then
          !  If I am converged for the first time
          !  then flip my flag and send the last status update
          pf%state%status = PF_STATUS_CONVERGED
          call pf_send_status(pf, send_tag)
       end if
    else
       !  I am not converged, send the news
       pf%state%status = PF_STATUS_ITERATING
       call pf_send_status(pf, send_tag)
    end if

  end subroutine pf_check_convergence_block

  !> Subroutine to test residuals to determine if the current processor has converged.
  subroutine pf_check_residual(pf, level_index, residual_converged)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    logical,           intent(out)   :: residual_converged  !! Return true if residual is below tolerances

    residual_converged = .false.

    ! Check to see if absolute tolerance is met
    if   (pf%levels(level_index)%residual     < pf%abs_res_tol)  then
       if (pf%debug) print*, 'DEBUG --',pf%rank, 'residual tol met',pf%levels(level_index)%residual
       residual_converged = .true.
    end if

  end subroutine pf_check_residual

end module pf_mod_MGRIT
