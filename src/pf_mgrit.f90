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
  
contains
  !>  Do the MGRIT algorithm
  subroutine pf_MGRIT_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout), target   :: pf   !!  The complete PFASST structure
    class(pf_encap_t), intent(in   )           :: q0   !!  The initial condition
    real(pfdp),        intent(in   )           :: dt   !!  The time step for each processor
    real(pfdp),        intent(in   )           :: tend !!  The final time of run
    integer,           intent(in   ), optional :: nsteps  !!  The number of time steps
    class(pf_encap_t), intent(inout), optional :: qend    !!  The computed solution at tend

    !  Local variables
    integer :: nproc  !!  Total number of processors
    integer :: nsteps_loc  !!  local number of time steps
    real(pfdp) :: tend_loc !!  The final time of run
    integer :: ierr 


    ! make a local copy of nproc
    nproc = pf%comm%nproc

    !>  Set the number of time steps to do
    !!  The user can either pass in the number of time steps or
    !!  pass in the time step size and length of run
    if (present(nsteps)) then
      nsteps_loc = nsteps
      tend_loc=real(nsteps_loc*dt,pfdp)
    else
      nsteps_loc = ceiling(tend/dt)
      !  Do  sanity check on steps
      if (abs(real(nsteps_loc,pfdp)-tend/dt) > dt/1d-7) then
        print *,'dt=',dt
        print *,'nsteps=',nsteps_loc
        print *,'tend=',tend
       call pf_stop(__FILE__,__LINE__,'Invalid nsteps ,nsteps=',nsteps)
      end if
    end if
    pf%state%nsteps = nsteps_loc


    !  do sanity checks on Nproc
    if (mod(nsteps,nproc) > 0)  call pf_stop(__FILE__,__LINE__,'nsteps must be multiple of nproc ,nsteps=',nsteps)

    !>  Allocate stuff for holding results 
    call initialize_results(pf)
    
    !>  Try to sync everyone
    call mpi_barrier(pf%comm%comm, ierr)

    !> Start timer
    if (pf%save_timings > 0) call pf_start_timer(pf, T_TOTAL)
    if (present(qend)) then
       call pf_MGRIT_block_run(pf, q0, dt, nsteps_loc,qend=qend)
    else
       call pf_MGRIT_block_run(pf, q0, dt,  nsteps_loc)
    end if

    !> End timer    
    if (pf%save_timings > 0) call pf_stop_timer(pf, T_TOTAL)

    !>  Output stats
    call pf_dump_stats(pf)


  end subroutine pf_MGRIT_run

  !>  MGRIT controller for block mode
  subroutine pf_MGRIT_block_run(pf, q0, dt, nsteps, qend,flags)
    use pf_mod_mpi, only: MPI_REQUEST_NULL
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend
    integer,           intent(in   ), optional :: flags(:)

    class(pf_level_t), pointer :: lev  !!  pointer to the one level we are operating on
    integer                   :: j, k, ierr
    integer                   :: nblocks !!  The number of blocks of steps to do
    integer                   :: nproc   !!  The number of processors being used
    integer                   :: level_index_c !!  Coarsest level in V (Lambda)-cycle
    integer                   :: level_max_depth !!  Finest level in V-cycle
!    integer::  nsteps_c,nsteps_f  

    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt

    ! set finest level to visit in the following run
    pf%state%finest_level = pf%nlevels

    !  pointer to finest  level to start
    lev => pf%levels(pf%state%finest_level)

    !  Stick the initial condition into q0 (will happen on all processors)
    call lev%q0%copy(q0, flags=0)


    nproc = pf%comm%nproc
    nblocks = nsteps/nproc

    !  Decide what the coarsest level in the V-cycle is
    level_index_c=1
    if (.not. pf%Vcycle)     level_index_c=pf%state%finest_level

    do k = 1, nblocks   !  Loop over blocks of time steps
       if (pf%save_timings > 1) call pf_start_timer(pf, T_BLOCK)
       call call_hooks(pf, -1, PF_PRE_BLOCK)
       ! print *,'Starting  step=',pf%state%step,'  block k=',k
       ! Each block will consist of
       !  1.  predictor
       !  2.  Vcycle until max iterations, or tolerances met
       !  3.  Move solution to next block

       !  Reset some flags
       !>  For initial block, this is done when initial conditions are set
       pf%state%iter    = 0
       pf%state%itcnt   = 0
       pf%state%mysteps = 0
       pf%state%status  = PF_STATUS_PREDICTOR
       pf%state%pstatus = PF_STATUS_PREDICTOR
       pf%comm%statreq  = MPI_REQUEST_NULL
       pf%state%pfblock = k
       pf%state%sweep = 1   !  Needed for compatibility of residual storage       


       if (k > 1) then
          !>  When starting a new block, broadcast new initial conditions to all procs
          if (pf%debug) print *,'DEBUG-rank=',pf%rank, ' at barrier at k=',k
          call mpi_barrier(pf%comm%comm, ierr)
          if (pf%debug) print *,'DEBUG-rank=',pf%rank, ' past barrier at k=',k
          if (nproc > 1)  then
             call lev%qend%pack(lev%send)    !!  Pack away your last solution
             call pf_broadcast(pf, lev%send, lev%mpibuflen, pf%comm%nproc-1)
             call lev%q0%unpack(lev%send)    !!  Everyone resets their q0
          else
             call lev%q0%copy(lev%qend, flags=0)    !!  Just stick qend in q0
          end if

          !>  Update the step and t0 variables for new block
          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt
       end if

       !> Call the predictor to get an initial guess on all levels and all processors
       call pf_MGRIT_predictor(pf, pf%state%t0, dt, flags)
       ! After the predictor, the residual and delta_q0 are just zero
       call pf_set_delta_q0(pf,1,0.0_pfdp)       
       call pf_set_resid(pf,pf%nlevels,0.0_pfdp)       
       call call_hooks(pf, -1, PF_POST_ITERATION)       !  This is the zero iteration
       
       if (pf%nlevels > 1) then
          !>  Start the MGRIT iterations
          do j = 1, pf%niters
             call call_hooks(pf, -1, PF_PRE_ITERATION)
             if (pf%save_timings > 1) call pf_start_timer(pf, T_ITERATION)
             
             pf%state%iter = j
             
             !  Do a v_cycle
             call pf_MGRIT_v_cycle(pf, k, pf%state%t0, dt, 1,2)
             
             !  Check for convergence
             call pf_check_convergence_block(pf, pf%state%finest_level, send_tag=1111*k+j)
             
             if (pf%save_timings > 1) call pf_stop_timer(pf, T_ITERATION)
             call call_hooks(pf, -1, PF_POST_ITERATION)
             
             !  If we are converged, exit block
             if (pf%state%status == PF_STATUS_CONVERGED)  then
                call call_hooks(pf, -1, PF_POST_CONVERGENCE)
                call pf_set_iter(pf,j)                 
                exit
             end if
          end do  !  Loop over j, the iterations in this block
       if (pf%save_timings > 1) call pf_stop_timer(pf, T_BLOCK)
       call call_hooks(pf, -1, PF_POST_BLOCK)
    end if
    
    
    end do !  Loop over the blocks
    call call_hooks(pf, -1, PF_POST_ALL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       call qend%copy(lev%qend, flags=0)
    end if
  end subroutine pf_MGRIT_block_run 

  !> Execute a MGRIT V-cycle (iteration)
  !!  It is assumed that we have two levels and two nodes here
  !!  When this is called the previous coarse integrator result should be stored in Q(1)
  !!  and the MGRIT iteration in qend (both on coarse level).  If this is called
  !!  directly after the predictor, these will be the same thing
  subroutine pf_MGRIT_v_cycle(pf, iteration, t0, dt, level_index_c, level_index_f, flags)

    !> FCF-relaxation on finest grid
    zero_rhs_flag = 1
    interp_flag = 0
    zero_c_pts_flag = 0
    qc_len = size(lev%qc)
    do i,qc_len
       call lev%qc_prev(i)%copy(lev%qc(i))
    end do
    call F_Relax(pf, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
    call C_Relax(pf, level_index, zero_rhs_flag, interp_flag)
    do i,qc_len
       call lev%qc_prev(i)%copy(lev%qc(i))
    end do
    call F_Relax(pf, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
    zero_rhs_flag = 0
    do level_index = nlevels-1,2,-1
        f_lev => pf%levels(level_index+1)
        lev => pf%levels(level_index) 

        !> Restrict residual
        call Restrict(pf, level_index, level_index_f)
        qc_len = size(f_lev%qc)
        do i,qc_len
           call f_lev%qc(i)%copy(f_lev%qc_prev(i))
        end do

        !> FCF-relaxation on intermediate grids
        interp_flag = 0
        zero_c_pts_flag = 1
        qc_len = size(lev%qc)
        do i,qc_len
           call lev%qc_prev(i)%copy(lev%qc(i))
        end do
        call F_Relax(pf, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
        interp_flag = 1
        zero_c_pts_flag = 0
        call C_Relax(pf, level_index, zero_rhs_flag, interp_flag)
        do i,qc_len
           call lev%qc_prev(i)%copy(lev%qc(i))
        end do
        call F_Relax(pf, level_index, zero_rhs_flag, interp_flag, zero_c_pts_flag, Q0)
    end do
    !> Coarsest grid solve
    level_index = nlevels
    call ExactSolve(pf, level_index)

  end subroutine pf_MGRIT_v_cycle

  subroutine F_Relax(pf, level_index, interp_flag, zero_rhs_flag, zero_init_c_pts_flag, Q0)

     nlevels = pf%nlevels
     dt = lev%ulevel%dt
     t0 = pf%state%t0 
     finest_lev => pf%levels(nlevels)

     q_zero_flag = 0;
     do i = 1, size(f_blocks)
        !> Check if inital values of C-points are zero
        if (zero_c_pts_flag .eq. .true.) then
           call lev%q0%setval(0.0_pfdp)
           q_zero_flag = 1;
        else
           if (i .eq. 1) then
               !> If we are at the first F-point on the fine level, use ODE initial value.
               !> Otherwise, set initial value to zero.
               if (k .eq. nlevels)
                   call lev%q0%copy(Q0);
               else
                   call lev%q0%setval(0.0_pfdp)
                   q_zero_flag = 1;
               end if
           else
               call lev%q0%copy(lev%qc_prev(i-1))
           end if
        end if
        do n = 1, size(f_blocks(i))
           j = f_blocks(i,n)
           if (q_zero_flag .eq. .true.) then
              !> Ifi intial C-points are zero, just copy RHS into q
              if (zero_rhs_flag .eq. .false.)
                 call lev%qend%axpy(1.0_pfdp, lev%g(j))
              end if
              q_zero_flag = 0;
           else 
              !> Do a single step
              call Point_Relax(pf, lev, level_index, j, lev%q0, lev%qend)
              !> Add g to the result
              if (zero_rhs_flag .eq. .false.) then
                 call lev%qend%axpy(1.0_pfdp, lev%g(j))
              end if
           end if

           !> Interpolate to finest level
           if (interp_flag .eq. .true.) then
              ii = lev%interp_map(j)
              call finest_lev%qc(ii)%axpy(1.0_pfdp, lev%qend)
           end if
           call lev%q0%copy(lev%qend)
        end do
        !> Store result (used for C-relaxation)
        call lev%qc(i)%copy(lev%qend)
     end do
  end subroutine F_Relax

  subroutine C_Relax(pf, level_index)
     do i = 1,size(c_pts)
        j = c_pts(i)
        call lev%q0%copy(lev%qc(i))
        call Point_Relax(pf, lev, level_index, j, lev%q0, lev%qend)
        call lev%qc(i)%copy(lev%qend)
        if (zero_rhs_flag .eq. .false.) then
           call lev%qc(i)%axpy(1.0_pfdp, lev%g(j))
        end if
        if (interp_flag .eq. .true.) then
           ii = lev%interp_map(j)
           call finest_lev%qc(ii)%axpy(1.0_pfdp, lev%qc(i))
        end if
    end do
  end subroutine C_Relax

  subroutine ExactSolve(pf, level_index)
     do i = 1,Nt(k)
         ii = restrict_map(i)
         if (k+1 .eq. nlevels)
             zero_rhs_flag = 1;
             call InjectRestrictPoint(c_lev%g(i), f_lev%qc(i), f_lev%qc_prev(i), A, null, ht(k+1), zero_rhs_flag);
         else
             zero_rhs_flag = 0;
             call InjectRestrictPoint(c_lev%g(i), f_lev%qc(i), f_lev%qc_prev(i), A, f_lev%g(ii), ht(k+1), zero_rhs_flag);
         end if
         call c_lev%qc(i)%copy(c_lev%qc_prev(i));
         if (i == 1)
            call lev%qend%copy(lev%g(i))   
         else
            call Point_Relax(pf, lev, level_index, i, lev%q0, lev%qend)
            call lev%qend%axpy(1.0_pfdp, lev%g(i))
         end if
         ii = lev%interp_map(i)
         call finest_lev%qc(ii)%axpy(1.0_pfdp, lev%qend)
     end do
  end subroutine ExactSolve


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> MATLAB TRANSLATION STOPPED HERE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Restrict(pf, level_index_c, level_index_f)
     do i = 1:,Nt(level_index_c)
         i_f = restrict_map(i)
         if (level_index_f == nlevels)
            zero_rhs_flag = 1;
         else
            zero_rhs_flag = 0;
         end if
         call InjectRestrictPoint(pf, level_index_c, level_index_f, i_c, i_f, zero_rhs_flag);
     end
  end subroutine Restrict

  subroutine InjectRestrictPoint(pf, level_index_c, level_index_f, i_c, i_f, zero_rhs_flag)
     call c_lev%g(i_c)%copy(f_lev%qc_prev(i_c))
     call c_lev%g(i_c)%scale(-1.0_pfdp)
     call f_lev%q0%copy(q_f)
     call Point_Relax(pf, lev, level_index, i, f_lev%q0, f_lev%qend) 
     call g_c%axpy(1.0_pfdp, f_lev%qend)
     if (zero_rhs_flag == 0)
         call g_c%axpy(1.0_pfdp, g_f);
     end if
  end function InjectRestrictPoint

  subroutine Point_Relax(pf, lev, level_index, n, q0, qend)
     dt = lev%ulevel%dt
     t0 = pf%state%t0
     t0n = t0 + dt * real(n-1,pfdp)
     lev%ulevel%stepper%do_n_steps(pf, level_index, t0n, q0, qend, dt, 1)
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