!!  Main controllers for optimal control problems
!
! This file is part of LIBPFASST.
!
!> Module of parallel PFASST routines for optimal control problems.
module pf_mod_parallel_oc
  use pf_mod_pfasst
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_comm
  implicit none
contains

  subroutine pf_predictor_oc(pf, t0, dt, flags)
    type(pf_pfasst_t), intent(inout), target :: pf     !! PFASST main data structure
    real(pfdp),        intent(in   )         :: t0     !! Initial time of this processor
    real(pfdp),        intent(in   )         :: dt     !! time step
    integer,           intent(in   ), optional :: flags(:)  !!  User defined flags

    class(pf_level_t), pointer :: c_lev
    class(pf_level_t), pointer :: f_lev
    integer                   :: k               !!  Loop indices
    integer                   :: level_index     !!  Local variable for looping over levels
    real(pfdp)                :: t0k             !!  Initial time at time step k
    integer :: which, dir, send_tag, burnin_sweeps, my_coarse_sweeps

    which = 1                   ! standard: predict and sweep forward-in-time
    dir = 1                     ! for MPI communication, standard is forward-in-time
    if(present(flags)) then
       if(flags(1)==2) then
          which = 2                ! if we are computing an adjoint, predict and sweep backward-in-time
          dir = -1                 ! communication has to be backwards as well
       end if
       if(flags(1)==0) which = 0   ! sweep forward and backward simultaneously on two components, communication only forwards
    end if
    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call pf_start_timer(pf, T_PREDICTOR)

    if (pf%debug) print*, 'DEBUG --', pf%rank, 'beginning predictor'

    !! Step 1. Getting the  initial condition on the finest level at each processor
    !!         If we are doing multiple levels, then we need to coarsen to fine level
    f_lev => pf%levels(pf%state%finest_level)
    if (pf%q0_style < 2) then  !  Spread q0 to all the nodes
       if( (which == 0) .or. (which == 1)) call f_lev%ulevel%sweeper%spreadq0(pf,pf%state%finest_level, t0, 1, pf%state%step+1)
       if( (which == 0) .or. (which == 2)) call f_lev%ulevel%sweeper%spreadq0(pf,pf%state%finest_level, t0+dt, 2, pf%state%step+1)
    endif


    !!  Step 2:   Proceed fine to coarse levels coarsening the fine solution and computing tau correction
    if (pf%debug) print*,  'DEBUG --', pf%rank, 'do coarsen  in predictor'
    if (pf%state%finest_level > 1) then
       do level_index = pf%state%finest_level, 2, -1
          f_lev => pf%levels(level_index);
          c_lev => pf%levels(level_index-1)
          call pf_residual(pf, f_lev%index, dt, which)
          if( (which == 0) .or. (which == 1)) &
               call f_lev%ulevel%restrict(f_lev, c_lev, f_lev%q0, c_lev%q0, t0, flags=1)
          if( (which == 0) .or. (which == 2)) &
               call f_lev%ulevel%restrict(f_lev, c_lev, f_lev%qend, c_lev%qend, t0+dt, flags=2)
          call restrict_time_space_fas(pf, t0, dt, level_index, flags=which)  !  Restrict
          call save(pf, c_lev, which)
       end do  !  level_index = pf%state%finest_level, 2, -1
    end if
    level_index = 1
    c_lev => pf%levels(1)

  if (pf%q0_style < 3) then

    ! Step 3. Do the "Burn in" step on the coarse level to make the coarse values consistent
    !         (this is skipped if the fine initial conditions are already consistent)
    ! The first processor does nothing, the second does one set of sweeps, the 2nd two, etc
    ! Hence, this is skipped completely if nprocs=1
    if (pf%q0_style .eq. 0) then  !  The coarse level needs burn in
       if (pf%debug) print*,  'DEBUG --', pf%rank, 'do burnin in pred', ' RK_pred', pf%RK_pred, ' PFASST_pred', pf%PFASST_pred
       !! If RK_pred is true, just do some RK_steps
       if (pf%RK_pred .or. which==2) then  !  Use Runge-Kutta to get the coarse initial data
          !  Get new initial conditions
          call pf_recv(pf, c_lev, 100000+pf%rank, .true., dir=dir,which=which)

          !  Do a RK_step
          call c_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1, which )
          !  Send forward/backward
!          if (dir == 1) send_tag = 100000+pf%rank+1
          !          if (dir == 2) send_tag = 100000+pf%rank-1
          send_tag = 100000+pf%rank+dir          
          call pf_send(pf, c_lev, send_tag, .false., dir=dir,which=which)
       else  !  Normal PFASST burn in
          burnin_sweeps = pf%rank+1
          if(which == 2) then
            if(pf%rank == 0) &
              print *, 'WARNING --- normal PFASST burn in is not suitable for adjoint as rhs cannot be evaluated for [t0k, t0k+dt]'
            burnin_sweeps = pf%comm%nproc-pf%rank
          end if
          if (pf%debug) print *, 'DEBUG ---', pf%rank, 'which = ', which, 'burnin_sweeps = ', burnin_sweeps
          do k = 1, burnin_sweeps !pf%rank + 1
             pf%state%iter = -k
             t0k = t0-(pf%rank)*dt + (k-1)*dt   ! Remember t0=pf%rank*dt is the beginning of this time slice so t0-(pf%rank)*dt is 0
                                                ! and we iterate up to the correct time step.
                                                ! for optimal control problem t, t0k has no influence on f_eval, so there this does something else
             if(which == 2) t0k = t0 + (burnin_sweeps-1)*dt - (k-1)*dt
             if(pf%debug) print *, 'DEBUG ----', pf%rank, 't0k = ', t0k

             ! Get new initial value (skip on first iteration)
             if (k > 1) then
                if ((which == 0) .or. (which == 1)) call c_lev%q0%copy(c_lev%qend, 1)
!                 if ((which == 0) .or. (which == 2)) call c_lev%qend%copy(c_lev%q0, 2) ! for which==0, we solve with zero terminal conditions,
                                                                                            ! but q0,2 is not zero (source term due to state sweeps)
                if (which == 2) call c_lev%qend%copy(c_lev%q0, 2)
                ! If we are doing PFASST_pred, we use the old values at nodes, otherwise spread q0
                if (.not. pf%PFASST_pred) then
                   if( (which == 0) .or. (which == 1)) call c_lev%ulevel%sweeper%spreadq0(pf,level_index, t0k, 1, pf%state%step+1)
!                    if( (which == 0) .or. (which == 2)) call c_lev%ulevel%sweeper%spreadq0(pf,level_index, t0k+dt, 2, pf%state%step+1)
                   if( which == 2) call c_lev%ulevel%sweeper%spreadq0(pf,level_index, t0k+dt, 2, pf%state%step+1)
                end if
             end if
             !  Do some sweeps
             if( which == 0 .or. which == 1 ) call c_lev%ulevel%sweeper%sweep(pf, level_index, t0k, dt, pf%nsweeps_burn, 1) ! was: 1 not which
             if( which == 2 )                 call c_lev%ulevel%sweeper%sweep(pf, level_index, t0k, dt, pf%nsweeps_burn, 2)
          end do
       endif  !  RK_pred
    end if  ! (q0_style .eq. 0)

    if (pf%q0_style > 0) then
      my_coarse_sweeps = pf%rank+1 ! for warm start do pipelining
      if(which == 2) my_coarse_sweeps = pf%comm%nproc-pf%rank
    else
      my_coarse_sweeps = c_lev%nsweeps_pred
    end if

    ! Step 4: Now we have everyone burned in, so do some coarse sweeps
    ! Modification: each processor does sweeps according to its rank
    if (pf%state%finest_level > 1) then
      pf%state%pstatus = PF_STATUS_ITERATING
      pf%state%status  = PF_STATUS_ITERATING
      if (pf%debug) print*,  'DEBUG --', pf%rank, 'do sweeps  in predictor', ' Pipeline_pred', pf%Pipeline_pred
      level_index=1
      c_lev => pf%levels(level_index)

      if (pf%Pipeline_pred) then
        do k = 1, my_coarse_sweeps !c_lev%nsweeps_pred
          pf%state%iter =-(pf%rank + 1) -k

          !  Get new initial conditions
          call pf_recv(pf, c_lev, c_lev%index*110000+pf%rank+k, .true., dir=dir,which=which)

          !  Do a sweep
          call c_lev%ulevel%sweeper%sweep(pf, c_lev%index, t0, dt, 1, which)
          !  Send forward/backward
!          if (dir == 1) send_tag = c_lev%index*1110000+pf%rank+1+k
          !          if (dir == 2) send_tag = c_lev%index*1110000+pf%rank-1+k
          send_tag = c_lev%index*1110000+pf%rank+dir+k          
          call pf_send(pf, c_lev, send_tag, .false., dir=dir,which=which)
       end do ! k = 1, c_lev%nsweeps_pred-1
      else  !  Don't pipeline
        !  Get new initial conditions
        call pf_recv(pf, c_lev, c_lev%index*100000+pf%rank, .true., dir=dir,which=which)

        !  Do sweeps
!         if(which == 0 .or. which == 1) call c_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, c_lev%nsweeps_pred, 1) !1 ! why only state?
        if(which == 0 .or. which == 1) call c_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, my_coarse_sweeps, which) !1 ! why only state?
!         if(which == 2)                 call c_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, c_lev%nsweeps_pred, 2) !which
        if(which == 2)                 call c_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, my_coarse_sweeps, 2) !which
        !  Send forward/backward
!        if (dir == 1) send_tag = c_lev%index*100000+pf%rank+1
!        if (dir == 2) send_tag = c_lev%index*100000+pf%rank-1
        send_tag = c_lev%index*100000+pf%rank+dir
        call pf_send(pf, c_lev, send_tag, .false., dir=dir,which=which)
      endif  ! (Pipeline_pred .eq. .true) then
    end if ! pf%state%finest_level > 1

    !  Step 5:  Return to fine level sweeping on any level in between coarsest and finest
    if (pf%debug) print*,  'DEBUG --', pf%rank, 'returning to fine level in predictor'
    do level_index = 2, pf%state%finest_level  !  Will do nothing with one level
       f_lev => pf%levels(level_index);
       c_lev => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, c_lev%Finterp, flags=which)
       if ((which == 0) .or. (which == 1)) then
          call f_lev%qend%copy(f_lev%Q(f_lev%nnodes), flags=1)
           if (pf%rank /= 0) call interpolate_q0(pf, f_lev, c_lev, flags=1)
       end if
       if (which == 2) then ! for which==0, qend never changes, so don't need to interpolate
          call f_lev%q0%copy(f_lev%Q(1), flags=2)
          if (pf%rank /= pf%comm%nproc-1) call interpolate_qend(pf, f_lev, c_lev)
       end if
       !  Do sweeps on level unless we are at the finest level
       if (level_index < pf%state%finest_level) then
          if ((which == 0) .or. (which == 1)) &
            call f_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev%nsweeps_pred, which) !which was 1
          if (which == 2)                     &
            call f_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev%nsweeps_pred, 2)
       end if
    end do

  end if

    call pf_stop_timer(pf, T_PREDICTOR)
    call call_hooks(pf, -1, PF_POST_PREDICTOR)

    pf%state%iter   = 0
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING

    if (pf%debug) print*,  'DEBUG --', pf%rank, 'ending predictor'
  end subroutine pf_predictor_oc


  !> Subroutine to test residuals to determine if the current processor has converged.
  subroutine pf_check_residual_oc(pf, level_index, residual_converged)
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

  end subroutine pf_check_residual_oc

  !>
  !> Test residuals to determine if the current processor has converged,
  !> adapted to optimal control. Can probably be removed, when pf_pfasst_block_oc
  !> is changed to use pf_check_convergence of pf_check_convergence_old.
  !>
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor hasn't converged yet either,
  !> regardless of the residual.
  !>
!   subroutine pf_check_convergence_oc(pf, k, residual,converged, flags)
  subroutine pf_check_convergence_oc(pf, level_index, send_tag, flags)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    integer,           intent(in)    :: send_tag
!     real(pfdp),        intent(inout) :: residual
!     integer,           intent(in)    :: k
!     logical,           intent(out)   :: converged   !!  True if this processor is done
    integer, optional, intent(in)    :: flags
!     real(pfdp)     :: residual1
    integer :: dir, which
    logical :: residual_converged, converged

    converged = .false.


    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    ! in first sweep: always continue
    if (pf%state%iter == 1) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    which = 1
    if (present(flags)) which = flags
    ! send forward by default, even if sweeping on both components; send backwards if sweeping on p only
    dir = 1
    if(which == 2) dir = -1

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)

    !> Check to see if tolerances are met
    call pf_check_residual_oc(pf, level_index, residual_converged)

    !>  Until I hear the previous processor is done, recieve it's status
    if (pf%state%pstatus /= PF_STATUS_CONVERGED) call pf_recv_status(pf, send_tag, dir=dir)

    !>  Check to see if I am converged
    converged = .false.
    if (residual_converged) then
       if (pf%rank == 0 .and. dir==1) then
          converged = .true.
       elseif (pf%rank == pf%comm%nproc-1 .and. dir==-1) then
          converged = .true.
       else  !  I am not the first/last processor, so I need to check the previous one
          if (pf%state%pstatus == PF_STATUS_CONVERGED) converged = .true.
       end if
    end if ! (residual_converged)


    !> Assign status and send it forward
    if (converged) then
       if (pf%state%status == PF_STATUS_ITERATING) then
          !  If I am converged for the first time
          !  then flip my flag and send the last status update
          pf%state%status = PF_STATUS_CONVERGED
          call pf_send_status(pf, send_tag, dir=dir)
       end if
    else
       !  I am not converged, send the news
       pf%state%status = PF_STATUS_ITERATING
       call pf_send_status(pf, send_tag, dir=dir)
    end if

    call call_hooks(pf, 1, PF_POST_CONVERGENCE)
  end subroutine pf_check_convergence_oc


  !>  Routine to do the pfasst iterations for optimal control problems on one block of processors until completion.
  !>  Each processor will do either a fixed number of iterations, or iterate until a tolerance is met
  !>  On calling, it is assumed that the levels are already loaded with the initial guesses
  !>
  subroutine pf_pfasst_block_oc(pf, dt, nsteps, predict, flags, step)
    use pf_mod_mpi, only: MPI_REQUEST_NULL
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: nsteps
    logical,           intent(in)    :: predict
    integer, optional, intent(in)    :: flags    !0 (default): sweep on y and p, 1: just y, 2: just p
    integer, optional, intent(in)    :: step
    ! not yet clear how to handle send and receive for forward and backward combined

    type(pf_level_t), pointer :: fine_lev, coarse_lev
    integer                   :: k, j, l, which, pred_flags(1), dir, ierror !dir to choose forward or backward send
    real(pfdp)                :: residual

    logical :: converged, qbroadcast
    logical :: did_post_step_hook
    call pf_start_timer(pf, T_TOTAL)

    which = 1
    if (present(flags)) which = flags
    ! send forward by default, even if sweeping on both components; send backwards if sweeping on p only
    dir = 1
    if(which == 2) dir = -1
    pred_flags(1) = which

    if( present(step) ) then
      pf%state%step    = step
    else
      pf%state%step    = pf%rank
    end if

!     print *, pf%state%step

    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
!     pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = MPI_REQUEST_NULL
    pf%state%nsteps  = nsteps
!     pf%state%component = which


    residual = -1
    did_post_step_hook = .false.


!       if (pf%state%status == PF_STATUS_PREDICTOR) then
!         !print *, 'pf%state%status == PF_STATUS_PREDICTOR', pf%state%t0, dt, which
        if (predict) then
          !print *, 'calling predictor'
           call pf_predictor_oc(pf, pf%state%t0, dt, pred_flags)
        end if
!       end if
        call call_hooks(pf, -1, PF_POST_ITERATION)


    pf%state%iter = 0
    pf%state%status  = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING


    k = pf%state%pfblock  !
    !pf%state%pfblock = k ! has to be set in pf_optimization_flex to current step
                          ! this is relevant for save_residuals
    do j = 1, pf%niters
      call pf_start_timer(pf, T_ITERATION)
      call call_hooks(pf, -1, PF_PRE_ITERATION)

      pf%state%iter = j

      !  Do a v_cycle
!       call pf_v_cycle(pf, k, pf%state%t0, dt, 1 ,pf%nlevels)
      call pf_v_cycle_oc(pf, j, pf%state%t0, dt, 1, pf%state%finest_level,  which)

      !  Check for convergence
      call pf_check_convergence_oc(pf, pf%state%finest_level, send_tag=1111*k+j, flags=dir)

      call pf_stop_timer(pf, T_ITERATION)
      call call_hooks(pf, -1, PF_POST_ITERATION)

      !  If we are converged, exit block
      if (pf%state%status == PF_STATUS_CONVERGED) then
         call call_hooks(pf, -1, PF_POST_CONVERGENCE)
         exit
      end if
      
    end do  !  Loop over the iteration in this block
    pf%state%itcnt = pf%state%itcnt + pf%state%iter
    call call_hooks(pf, -1, PF_POST_BLOCK)

    
    call pf_stop_timer(pf, T_TOTAL)
    ! call dump_results(pf%results)
    call dump_timingsl(pf%results,pf)
  end subroutine pf_pfasst_block_oc

  subroutine pf_v_cycle_oc(pf, iteration, t0, dt, level_index_c,level_index_f, flags)
  ! Execute a V-cycle between levels nfine and ncoarse

    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration
    integer,           intent(in)    :: level_index_c  !! Coarsest level of V-cycle
    integer,           intent(in)    :: level_index_f  !! Finest level of V-cycle
    integer, optional, intent(in)    :: flags

    type(pf_level_t), pointer :: f_lev, c_lev
    integer :: level_index, j, which, dir

    which = 1
    if(present(flags)) which = flags
    ! send forward by default, even if sweeping on both components; send backwards if sweeping on p only
    dir = 1
    if(which == 2) dir = -1 !


    !>  Post the nonblocking receives on the all the levels that will be recieving later
    !>    (for single level this will be skipped)
    do level_index = level_index_c+1, level_index_f
       f_lev => pf%levels(level_index)
       call pf_post(pf, f_lev, f_lev%index*10000+iteration, dir=dir)
       !call pf_post(pf, f_lev, f_lev%index*10000+iteration, dir)       
    end do

    !> move from fine to coarse doing sweeps
    do level_index = level_index_f, level_index_c+1, -1
       f_lev => pf%levels(level_index);
       c_lev => pf%levels(level_index-1)
       call f_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev%nsweeps, which)
       call pf_send(pf, f_lev, level_index*10000+iteration, .false., dir=dir,which=which)
       call restrict_time_space_fas(pf, t0, dt, level_index, flags=which)
       call save(pf, c_lev, which)
    end do

    !> Do the coarsest level
    level_index=level_index_c
    f_lev => pf%levels(level_index)
    if (pf%pipeline_pred) then
       do j = 1, f_lev%nsweeps
          call pf_recv(pf, f_lev, f_lev%index*10000+iteration+j, .true., dir=dir,which=which)
          call f_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1, which)
          call pf_send(pf, f_lev, f_lev%index*10000+iteration+j, .false., dir=dir,which=which)
       end do
    else
       call pf_recv(pf, f_lev, f_lev%index*10000+iteration, .true., dir=dir,which=which)
       call f_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev%nsweeps, which)
       call pf_send(pf, f_lev, f_lev%index*10000+iteration, .false., dir=dir,which=which)
    endif


    !> Now move coarse to fine interpolating and sweeping
    do level_index = level_index_c+1,level_index_f
       f_lev => pf%levels(level_index);
       c_lev => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, c_lev%Finterp, flags=which)

       if ((flags .eq. 0) .or. (flags .eq. 1))  call f_lev%qend%copy(f_lev%Q(f_lev%nnodes), flags=1)
       if (flags .eq. 2)                        call f_lev%q0%copy(f_lev%Q(1), flags=2)

       call pf_recv(pf, f_lev, level_index*10000+iteration, .false., dir=dir,which=which)

       if (pf%rank /= 0) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          if ((which .eq. 0) .or. (which .eq. 1)) call interpolate_q0(pf, f_lev, c_lev, flags=1)
       end if
       if (pf%rank /= pf%comm%nproc-1) then
          if (which .eq. 2)                       call interpolate_qend(pf, f_lev, c_lev)
       end if

       ! don't sweep on the finest level since that is only done at beginning
       if (level_index < level_index_f) then
          call f_lev%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev%nsweeps, which)
       else  !  compute residual for diagnostics since we didn't sweep
          call pf_residual(pf, f_lev%index, dt, which)
       end if
    end do


  end subroutine pf_v_cycle_oc

end module pf_mod_parallel_oc
