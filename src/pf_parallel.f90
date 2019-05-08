!!  Routines that run the PFASST algorithm
!
! This file is part of LIBPFASST.
!

!> Module of routines to run  PFASST
module pf_mod_parallel
  use pf_mod_pfasst
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_comm
  use pf_mod_results
  implicit none
contains

  !>  This is the main interface to pfasst.
  !!  It examines the parameters and decides which subroutine to call
  !!  to execute the code correctly
  subroutine pf_pfasst_run(pf, q0, dt, tend, nsteps, qend, flags)
    type(pf_pfasst_t), intent(inout), target   :: pf   !!  The complete PFASST structure
    class(pf_encap_t), intent(inout   )           :: q0   !!  The initial condition
    real(pfdp),        intent(in   )           :: dt   !!  The time step for each processor
    real(pfdp),        intent(in   )           :: tend !!  The final time of run
    integer,           intent(in   ), optional :: nsteps  !!  The number of time steps
    class(pf_encap_t), intent(inout), optional :: qend    !!  The computed solution at tend
    integer,           intent(in   ), optional :: flags(:)!!  User defnined flags


    !  Local variables
    integer :: nproc  !!  Total number of processors
    integer :: nsteps_loc  !!  local number of time steps
    real(pfdp) :: tend_loc !!  The final time of run


    ! make a local copy of nproc
    nproc = pf%comm%nproc

    !>  Set the number of time steps to do
    !!  The user can either pass in the number of time steps or
    !!  pass in the time step size and length of run
    if (present(nsteps)) then
      nsteps_loc = nsteps
      tend_loc=dble(nsteps_loc*dt)
    else
      nsteps_loc = ceiling(tend/dt)
      !  Do  sanity check on steps
      if (abs(real(nsteps_loc,pfdp)-tend/dt) > dt/100.0) then
        print *,'dt=',dt
        print *,'nsteps=',nsteps_loc
        print *,'tend=',tend
        stop "Invalid nsteps"
      end if
    end if
    pf%state%nsteps = nsteps_loc

    !>  Allocate stuff for holding results
    call pf_initialize_results(pf)

    !  do sanity checks on Nproc
    if (mod(nsteps,nproc) > 0) stop "ERROR: nsteps must be multiple of nproc (pf_parallel.f90)."

    if (present(qend)) then
       call pf_block_run(pf, q0, dt, nsteps_loc,qend=qend,flags=flags)
    else
       call pf_block_run(pf, q0, dt,  nsteps_loc,flags=flags)
    end if

    call pf_dump_results(pf)

    !>   deallocate results data
    call pf_destroy_results(pf)

    !  What we would like to do is check for
    !  1.  nlevels==1  and nprocs ==1 -> Serial SDC
    !      Predictor is either spreadQ or nothing
    !      Then we just call a loop on sweeps
    !      Communication is copy
    !  2.  nlevels > 1  and nprocs ==1 -> Serial MLSDC
    !      Predictor is needed to populate levels (or nothing)
    !      Then we just call a loop on MLSDC sweeps
    !      Communication is copy
    !  3.  nlevels == 1  and nprocs > 1 -> Pipelined SDC
    !      Predictor is just like PFASST, but on finest (only) level (or nothing)
    !  4.  nlevels > 1  and nprocs > 1 -> PFASST
  end subroutine pf_pfasst_run
  !
  !> PFASST Predictor.
  !>  Subroutine  to initialize the solution on each processor
  !!  The goal is to have a solution at each level and each node set to a consistent value
  !!  When this is called, the value of q0 at the fine level on each processor has been set somehow (see q0_style below)
  !!
  !! This can be broken down into four substeps
  !! 1. Get the  initial condition on the finest level at each node
  !! 2. Coarsen the initial condition to each coarser level with tau corrections
  !! 3. Do the "Burn in" step on the coarse level to make the coarse values consistent
  !!    (this is skipped if the fine initial conditions are already consistent)
  !! 4. Do some coarse grid sweeps to improve initial solutions on coarsest nodes
  !! 5. Interpolating coarse correction back to finer levels sweeping along the way.
  !!
  !! There are several parameters or flags that determine how things are done:
  !!  integer  q0_style:    can take 3 values
  !!           0:  Only the q0 at t=0 is valid  (default)
  !!           1:  The q0 at each processor is valid
  !!           2:  q0 and all nodes at each processor is valid
  !! logical  PFASST_pred:  If true, the burn-in step uses the "PFASST predictor" trick
  !! integer  nsweeps_burn: Determines how many sweeps are done on the coarse level during burn in
  !! integer  nsweeps_pred: Determines how many sweeps are done at the coarse level after burn in
  !! logical Pipeline_burn: True if coarse sweeps during burn in are pipelined  (meaningless if nsweeps_burn>1 on coarse level)
  !! logical Pipeline_pred: True if coarse sweeps after burn in are pipelined  (meaningless if nsweeps_pred>1 on coarse level)
  !!    Pipeline variables do nothing if there is only one processor
  !! logical  RK_pred:      If true, the coarse level is initialized with Runge-Kutta instead of the  PFASST burn in.
  !!                        We  will still do coarse sweeps after and correct finer levels
  !!
  !! The user defined flags(:) parameter is used to determine whether we are in a (standard) forward-in-time run (flags(1) == 1)
  !! or backward-in-time (for the adjoint) with a given terminal condition qend instead of initial condition q0  (flags(1) == 2).
  !! In the latter case, e.g., sweeper%spreadq0 has to do the correct thing (i.e., spread qend instead of q0).
  !!
  !! No time communication is performed during the predictor since all
  !! procesors can do the work themselves
  !!
  !!  The iteration count is reset to 0, and the status is reset to
  !!  ITERATING.
  subroutine pf_predictor(pf, t0, dt, flags)
    type(pf_pfasst_t), intent(inout), target :: pf     !! PFASST main data structure
    real(pfdp),        intent(in   )         :: t0     !! Initial time of this processor
    real(pfdp),        intent(in   )         :: dt     !! time step
    integer,           intent(in   ), optional :: flags(:)  !!  User defined flags

    class(pf_level_t), pointer :: c_lev_p
    class(pf_level_t), pointer :: f_lev_p     !!
    integer                   :: k               !!  Loop indices
    integer                   :: level_index     !!  Local variable for looping over levels
    real(pfdp)                :: t0k             !!  Initial time at time step k
    pf%state%iter = -1          

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    if (pf%debug) print*, 'DEBUG --', pf%rank, 'beginning predictor'
    !!
    !! Step 1. Getting the  initial condition on the finest level at each processor
    !!         If we are doing multiple levels, then we need to coarsen to fine level
    f_lev_p => pf%levels(pf%state%finest_level)

    if (pf%q0_style < 2) then  !  Spread q0 to all the nodes
       call f_lev_p%ulevel%sweeper%spreadq0(pf,pf%state%finest_level, t0)
    endif
    !!
    !!  Step 2:   Proceed fine to coarse levels coarsening the fine solution and computing tau correction
    if (pf%debug) print*,  'DEBUG --', pf%rank, 'do coarsen  in predictor'
    if (pf%state%finest_level > 1) then
       do level_index = pf%state%finest_level, 2, -1
          f_lev_p => pf%levels(level_index);
          c_lev_p => pf%levels(level_index-1)
          call pf_residual(pf, f_lev_p%index, dt)
          call f_lev_p%ulevel%restrict(f_lev_p, c_lev_p, f_lev_p%q0, c_lev_p%q0, t0)
          call restrict_time_space_fas(pf, t0, dt, level_index)  !  Restrict
          call save(pf,c_lev_p)
       end do  !  level_index = pf%state%finest_level, 2, -1
    else
      level_index = 1
      c_lev_p => pf%levels(1)
    end if

    !!
    !! Step 3. Do the "Burn in" step on the coarse level to make the coarse values consistent
    !!         (this is skipped if the fine initial conditions are already consistent)
    !! The first processor does nothing, the second does one set of sweeps, the third two, etc
    !! Hence, this is skipped completely if nprocs=1
    if (pf%debug) print*,  'DEBUG --', pf%rank, 'do burnin  in predictor'
    if (pf%q0_style .eq. 0) then  !  The coarse level needs burn in
       !! If RK_pred is true, just do some RK_steps
       if (pf%RK_pred) then  !  Use Runge-Kutta to get the coarse initial data
          !  Get new initial conditions
          call pf_recv(pf, c_lev_p, 100000+pf%rank, .true.)

          !  Do a RK_step
          call c_lev_p%ulevel%stepper%do_n_steps(pf, level_index,t0, c_lev_p%q0,c_lev_p%qend, dt, 1)       
          !  Send forward
          call pf_send(pf, c_lev_p,  100000+pf%rank+1, .false.)
          print *,'woo hoo'
       else  !  Normal PFASST burn in
          level_index=1
          c_lev_p => pf%levels(level_index)
          do k = 1, pf%rank + 1
             pf%state%iter = -k
             t0k = t0-(pf%rank)*dt + (k-1)*dt   ! Remember t0=pf%rank*dt is the beginning of this time slice so t0-(pf%rank)*dt is 0
                                                ! and we iterate up to the correct time step.
                                                ! for optimal control problem t, t0k has no influence on f_eval, so there this does something else

             ! Get new initial value (skip on first iteration)
             if (k > 1) then
                call c_lev_p%q0%copy(c_lev_p%qend,flags=0)
                ! If we are doing PFASST_pred, we use the old values at nodes, otherwise spread q0
                if (.not. pf%PFASST_pred) then
                   call c_lev_p%ulevel%sweeper%spreadq0(pf,level_index, t0k)
                end if
             end if
             !  Do some sweeps
             call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0k, dt,pf%nsweeps_burn)
          end do
       endif  !  RK_pred
    end if  ! (q0_style .eq. 0)

    !!
    !! Step 4: Now we have everyone burned in, so do some coarse sweeps
    if (pf%state%finest_level > 1) then
       if (pf%debug) print*,  'DEBUG --', pf%rank, 'do sweeps  in predictor', 'Pipeline_pred',pf%Pipeline_pred
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status = PF_STATUS_ITERATING
       if (pf%Pipeline_pred) then
          do k = 1, c_lev_p%nsweeps_pred
             pf%state%iter =-(pf%rank + 1) -k

             !  Get new initial conditions
             call pf_recv(pf, c_lev_p, c_lev_p%index*110000+pf%rank+k, .true.)

             !  Do a sweep
             call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1)
             !  Send forward
             call pf_send(pf, c_lev_p,  c_lev_p%index*110000+pf%rank+1+k, .false.)
          end do ! k = 1, c_lev_p%nsweeps_pred-1
       else  !  Don't pipeline
          !  Get new initial conditions
          call pf_recv(pf, c_lev_p, c_lev_p%index*110000+pf%rank, .true.)

          !  Do a sweeps
          call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, c_lev_p%nsweeps_pred)
          !  Send forward
          call pf_send(pf, c_lev_p,  c_lev_p%index*110000+pf%rank+1, .false.)
       endif  ! (Pipeline_pred .eq. .true) then
    end if

    if (pf%debug) print*,  'DEBUG --', pf%rank, 'returning to fine level in predictor'
    !!
    !!  Step 5:  Return to fine level sweeping on any level in between coarsest and finest
    do level_index = 2, pf%state%finest_level  !  Will do nothing with one level
       f_lev_p => pf%levels(level_index);
       c_lev_p => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, c_lev_p%Finterp)
       call f_lev_p%qend%copy(f_lev_p%Q(f_lev_p%nnodes), flags=0)
       if (pf%rank /= 0) call interpolate_q0(pf, f_lev_p, c_lev_p,flags=0)

       !  Do a sweep on level unless we are at the finest level
       if (level_index < pf%state%finest_level) then
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps_pred)
       end if
    end do

    call end_timer(pf, TPREDICTOR)
    call call_hooks(pf, -1, PF_POST_PREDICTOR)

    pf%state%iter   = 0
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING
    if (pf%debug) print*,  'DEBUG --', pf%rank, 'ending predictor'
  end subroutine pf_predictor


  !> Subroutine to test residuals to determine if the current processor has converged.
  subroutine pf_check_residual(pf, level_index, residual_converged)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    logical,           intent(out)   :: residual_converged  !! Return true if residual is below tolerances

    residual_converged = .false.

    ! Check to see if relative tolerance is met
    if (pf%levels(level_index)%residual_rel < pf%rel_res_tol) then
       if (pf%debug) print*, 'DEBUG --', pf%rank, ' residual relative tol met',pf%levels(level_index)%residual_rel
       residual_converged = .true.
    end if
    ! Check to see if relative tolerance is met
    if   (pf%levels(level_index)%residual     < pf%abs_res_tol)  then
       if (pf%debug) print*, 'DEBUG --',pf%rank, 'residual tol met',pf%levels(level_index)%residual
       residual_converged = .true.
    end if

  end subroutine pf_check_residual

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

  !

  !>  PFASST controller for block mode
  subroutine pf_block_run(pf, q0, dt, nsteps, qend,flags)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend
    integer,           intent(in   ), optional :: flags(:)

    class(pf_level_t), pointer :: lev_p  !!  pointer to the one level we are operating on
    integer                   :: j, k
    integer                   :: nblocks !!  The number of blocks of steps to do
    integer                   :: nproc   !!  The number of processors being used
    integer                   :: level_index_c !!  Coarsest level in V (Lambda)-cycle
    integer                   :: level_max_depth !!  Finest level in V-cycle


    call start_timer(pf, TTOTAL)

    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt

    ! set finest level to visit in the following run
    pf%state%finest_level = pf%nlevels

    !  pointer to finest  level to start
    lev_p => pf%levels(pf%state%finest_level)

    !  Stick the initial condition into q0 (will happen on all processors)
    call lev_p%q0%copy(q0, flags=0)


    nproc = pf%comm%nproc
    nblocks = nsteps/nproc

    !  Decide what the coarsest level in the V-cycle is
    level_index_c=1
    if (.not. pf%Vcycle)     level_index_c=pf%state%finest_level

    do k = 1, nblocks   !  Loop over blocks of time steps
       ! print *,'Starting  step=',pf%state%step,'  block k=',k
       ! Each block will consist of
       !  1.  predictor
       !  2.  Vcycle until max iterations, or tolerances met
       !  3.  Move solution to next block

       !  Reset some flags
       !>  When starting a new block, broadcast new initial conditions to all procs
       !>  For initial block, this is done when initial conditions are set

       !> Reset some flags
       pf%state%iter    = -1
       pf%state%itcnt   = 0
       pf%state%mysteps = 0
       pf%state%status  = PF_STATUS_PREDICTOR
       pf%state%pstatus = PF_STATUS_PREDICTOR
       pf%comm%statreq  = -66
       pf%state%pfblock = k


       if (k > 1) then
          if (nproc > 1)  then
             call lev_p%qend%pack(lev_p%send)    !!  Pack away your last solution
             call pf_broadcast(pf, lev_p%send, lev_p%mpibuflen, pf%comm%nproc-1)
             call lev_p%q0%unpack(lev_p%send)    !!  Everyone resets their q0
          else
             call lev_p%q0%copy(lev_p%qend, flags=0)    !!  Just stick qend in q0
          end if

          !>  Update the step and t0 variables for new block
          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt
       end if

       !> Call the predictor to get an initial guess on all levels and all processors
       call pf_predictor(pf, pf%state%t0, dt, flags)

       !>  Start the loops over SDC sweeps
       pf%state%iter = 0
       call call_hooks(pf, -1, PF_POST_ITERATION)

       call start_timer(pf, TITERATION)
       do j = 1, pf%niters

          call call_hooks(pf, -1, PF_PRE_ITERATION)

          pf%state%iter = j

          !  Do a v_cycle
          call pf_v_cycle(pf, k, pf%state%t0, dt, level_index_c, pf%state%finest_level)

          !  Check for convergence
          call pf_check_convergence_block(pf, pf%state%finest_level, send_tag=1111*k+j)

!          print *,pf%rank, ' post res'
          call call_hooks(pf, -1, PF_POST_ITERATION)

          !  If we are converged, exit block
          if (pf%state%status == PF_STATUS_CONVERGED)  then
             call pf%levels(pf%nlevels)%ulevel%sweeper%sweep(pf, pf%nlevels, pf%state%t0, dt, 1)
             exit             
          end if

       end do  !  Loop over the iteration in this bloc
       call call_hooks(pf, -1, PF_POST_CONVERGENCE)
       call end_timer(pf, TITERATION)
       call call_hooks(pf, -1, PF_POST_STEP)
    end do !  Loop over the blocks

    call end_timer(pf, TTOTAL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       call qend%copy(lev_p%qend, flags=0)
    end if
  end subroutine pf_block_run



  !> Execute a V-cycle between levels nfine and ncoarse
  subroutine pf_v_cycle(pf, iteration, t0, dt,level_index_c,level_index_f, flags)


    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration
    integer,           intent(in)    :: level_index_c  !! Coarsest level of V-cycle
    integer,           intent(in)    :: level_index_f  !! Finest level of V-cycle
    integer, optional, intent(in)    :: flags

    type(pf_level_t), pointer :: f_lev_p, c_lev_p
    integer :: level_index, j

    !>  Post the nonblocking receives on the all the levels that will be recieving later
    !>    (for single level this will be skipped)
    do level_index = level_index_c+1, level_index_f
       f_lev_p => pf%levels(level_index)
       call pf_post(pf, f_lev_p, f_lev_p%index*10000+iteration)
    end do


    !> move from fine to coarse doing sweeps
    do level_index = level_index_f, level_index_c+1, -1
       f_lev_p => pf%levels(level_index);
       c_lev_p => pf%levels(level_index-1)
       call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps)
       call pf_send(pf, f_lev_p, level_index*10000+iteration, .false.)
       call restrict_time_space_fas(pf, t0, dt, level_index)
       call save(pf,c_lev_p)
    end do


    ! Do the coarsest level
    level_index=level_index_c
    f_lev_p => pf%levels(level_index)
    if (pf%pipeline_pred) then
       do j = 1, f_lev_p%nsweeps
          call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration+j, .true.)
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1)
          call pf_send(pf, f_lev_p, f_lev_p%index*10000+iteration+j, .false.)
       end do
    else
       call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration, .true.)
       call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps)
       call pf_send(pf, f_lev_p, f_lev_p%index*10000+iteration, .false.)
    endif

    ! Now move coarse to fine interpolating and sweeping
    do level_index = level_index_c+1,level_index_f
       f_lev_p => pf%levels(level_index);
       c_lev_p => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, c_lev_p%Finterp)
       call f_lev_p%qend%copy(f_lev_p%Q(f_lev_p%nnodes), flags=0)
       call pf_recv(pf, f_lev_p, level_index*10000+iteration, .false.)

       if (pf%rank /= 0) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          call interpolate_q0(pf, f_lev_p, c_lev_p,flags=0)
       end if
       ! don't sweep on the finest level since that is only done at beginning
       if (level_index < level_index_f) then
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps)
       else  !  compute residual for diagnostics since we didn't sweep
          pf%state%sweep=1
          call pf_residual(pf, f_lev_p%index, dt)
       end if
    end do

  end subroutine pf_v_cycle

end module pf_mod_parallel
