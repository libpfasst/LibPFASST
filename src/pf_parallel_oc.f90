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
    type(pf_pfasst_t), intent(inout), target :: pf     !< PFASST main data structure
    real(pfdp),        intent(in   )         :: t0     !< Initial time of this processor
    real(pfdp),        intent(in   )         :: dt     !< time step
    integer,           intent(in   ), optional :: flags(:)  !<  User defined flags

    class(pf_level_t), pointer :: c_lev_p
    class(pf_level_t), pointer :: f_lev_p     !<
    integer                   :: j, k            !<  Loop indices
    integer                   :: level_index     !<  Local variable for looping over levels
    real(pfdp)                :: t0k             !<  Initial time at time step k
    integer :: which, dir

    which = 1                   ! standard: predict and sweep forward-in-time
    dir = 1                     ! for MPI communication, standard is forward-in-time
    if(present(flags)) then
       if(flags(1)==2) then
          which = 2                ! if we are computing an adjoint, predict and sweep backward-in-time
          dir = 2                  ! communication has to be backwards as well
       end if
       if(flags(1)==0) which = 0   ! sweep forward and backward simultaneously on two components, communication only forwards
    end if
    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    !! Step 1. Getting the  initial condition on the finest level at each processor
    !!         If we are doing multiple levels, then we need to coarsen to fine level
    f_lev_p => pf%levels(pf%nlevels)
    if (pf%q0_style < 2) then  !  Spread q0 to all the nodes
       if( (which == 0) .or. (which == 1)) call f_lev_p%ulevel%sweeper%spreadq0(f_lev_p, t0, 1, pf%state%step+1)
       if( (which == 0) .or. (which == 2)) call f_lev_p%ulevel%sweeper%spreadq0(f_lev_p, t0+dt, 2, pf%state%step+1)
    endif
    
    !!  Step 2:   Proceed fine to coarse levels coarsening the fine solution and computing tau correction
    if (pf%nlevels > 1) then  
       do level_index = pf%nlevels, 2, -1
          f_lev_p => pf%levels(level_index);
          c_lev_p => pf%levels(level_index-1)
          call pf_residual(pf, f_lev_p, dt, which)  
          call restrict_time_space_fas(pf, t0, dt, level_index, which)  !  Restrict
          call save(c_lev_p, which)
          if( (which == 0) .or. (which == 1)) call f_lev_p%ulevel%restrict(f_lev_p, c_lev_p, f_lev_p%q0, c_lev_p%q0, t0, 1)
          if( (which == 0) .or. (which == 2)) call f_lev_p%ulevel%restrict(f_lev_p, c_lev_p, f_lev_p%qend, c_lev_p%qend, t0+dt, 2)
       end do  !  level_index = pf%nlevels, 2, -1
    else
      level_index = 1
      c_lev_p => pf%levels(1)
    end if
    
    ! Step 3. Do the "Burn in" step on the coarse level to make the coarse values consistent
    !         (this is skipped if the fine initial conditions are already consistent)
    ! The first processor does nothing, the second does one set of sweeps, the 2nd two, etc
    ! Hence, this is skipped completely if nprocs=1
    if (pf%q0_style .eq. 0) then  !  The coarse level needs burn in
       !! If RK_pred is true, just do some RK_steps
       if (pf%RK_pred) then  !  Use Runge-Kutta to get the coarse initial data
          !  Get new initial conditions
          call pf_recv(pf, c_lev_p, 30000+pf%rank+k, .true., dir)

          !  Do a RK_step
          call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1, which )
          !  Send forward
          call pf_send(pf, c_lev_p,  30000+pf%rank+1+k, .false., dir)
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
                if ((which == 0) .or. (which == 1)) call c_lev_p%q0%copy(c_lev_p%qend, 1)
!                 if ((which == 0) .or. (which == 2)) call c_lev_p%qend%copy(c_lev_p%q0, 2) ! for which==0, we solve with zero terminal conditions,
                                                                                            ! but q0,2 is not zero (source term due to state sweeps)
                if (which == 2) call c_lev_p%qend%copy(c_lev_p%q0, 2)
                ! If we are doing PFASST_pred, we use the old values at nodes, otherwise spread q0
                if (.not. pf%PFASST_pred) then
                   if( (which == 0) .or. (which == 1)) call c_lev_p%ulevel%sweeper%spreadq0(c_lev_p, t0k, 1, pf%state%step+1)
!                    if( (which == 0) .or. (which == 2)) call c_lev_p%ulevel%sweeper%spreadq0(c_lev_p, t0k+dt, 2, pf%state%step+1)
                   if( which == 2) call c_lev_p%ulevel%sweeper%spreadq0(c_lev_p, t0k+dt, 2, pf%state%step+1)
                end if
             end if
             !  Do some sweeps
             if( which == 0 .or. which == 1 ) call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0k, dt,pf%nsweeps_burn, 1)
             if( which == 2 ) call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0k, dt,pf%nsweeps_burn, 2)
          end do
       endif  !  RK_pred
    end if  ! (q0_style .eq. 0)

    ! Step 4: Now we have everyone burned in, so do some coarse sweeps
    if (pf%Pipeline_pred) then
       do k = 1, c_lev_p%nsweeps_pred
          pf%state%pstatus = PF_STATUS_ITERATING
          pf%state%status = PF_STATUS_ITERATING
          pf%state%iter =-(pf%rank + 1) -k
       
          !  Get new initial conditions
          call pf_recv(pf, c_lev_p, c_lev_p%index*20000+pf%rank+k, .true., dir)
       
          !  Do a sweep
          call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1, which)
          !  Send forward
          call pf_send(pf, c_lev_p,  c_lev_p%index*20000+pf%rank+1+k, .false., dir)
       end do ! k = 1, c_lev_p%nsweeps_pred-1
    else  !  Don't pipeline
       !  Get new initial conditions
       call pf_recv(pf, c_lev_p, c_lev_p%index*20000+pf%rank, .true., dir)
       
       !  Do a sweeps
       if(which == 0 .or. which == 1) call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, c_lev_p%nsweeps_pred, 1) !which
       if(which == 2)                 call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, c_lev_p%nsweeps_pred, 2) !which
       !  Send forward
       call pf_send(pf, c_lev_p,  c_lev_p%index*20000+pf%rank+1, .false., dir)
    endif  ! (Pipeline_pred .eq. .true) then

    !  Step 6:  Return to fine level sweeping on any level in between coarsest and finest
    do level_index = 2, pf%nlevels  !  Will do nothing with one level
       f_lev_p => pf%levels(level_index);
       c_lev_p => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, c_lev_p%Finterp, which)
       if ((which == 0) .or. (which == 1)) call interpolate_q0(pf, f_lev_p, c_lev_p, 1)
!        if ((which == 0) .or. (which == 2)) call interpolate_qend(pf, f_lev_p, c_lev_p)
       if (which == 2) call interpolate_qend(pf, f_lev_p, c_lev_p) ! for which==0, qend never changes, so don't need to interpolate
       !  Do a sweep on unless we are at the finest level
       if (level_index < pf%nlevels) then
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps_pred, which)
       end if
    end do

    call end_timer(pf, TPREDICTOR)
    call call_hooks(pf, -1, PF_POST_PREDICTOR)

    pf%state%iter   = 0
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING
    
  end subroutine pf_predictor_oc

     !>
  !> Test residuals to determine if the current processor has converged,
  !> adapted to optimal control. Can probably be removed, when pf_pfasst_block_oc
  !> is changed to use pf_check_convergence of pf_check_convergence_old.
  !>
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor hasn't converged yet either,
  !> regardless of the residual.
  !>
  subroutine pf_check_convergence_oc(pf, k, dt, residual, energy, converged, flags)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual, energy
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: k
    logical,           intent(out)   :: converged   !<  True if this processor is done
    integer, optional, intent(in)    :: flags
    real(pfdp)     :: residual1
    integer :: dir, which
    
    converged = .false.

    
    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if
    
    which = 1
    if (present(flags)) which = flags
    ! send forward by default, even if sweeping on both components; send backwards if sweeping on p only
    dir = 1
    if(which == 2) dir = 2


    ! Check to see if tolerances are met
    residual1 = pf%levels(pf%nlevels)%residual
    if (pf%state%status == PF_STATUS_ITERATING .and. residual > 0.0d0) then
       if ( (abs(1.0_pfdp - abs(residual1/residual)) < pf%rel_res_tol) .or. &
            (abs(residual1)                          < pf%abs_res_tol) ) then
          pf%state%status = PF_STATUS_CONVERGED
       end if
    end if
        
!     !->why? how to do that more cleanly?
!     if (pf%state%status == PF_STATUS_ITERATING .and. residual >= 0.0d0) then    
!                 ! if do_mixed, adjoint on last time step will be constant zero, so residual will be zero
!                 ! need to stop in that case as well, but not in the very first iteration
!       if( abs(residual1) < pf%abs_res_tol ) then
!           pf%state%status = PF_STATUS_CONVERGED
!       end if
!     end if
!     !<-
    
    residual = residual1

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)
    if (pf%state%pstatus /= PF_STATUS_CONVERGED) call pf_recv_status(pf, 1+k, dir)

    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING .and. dir == 1) &
         pf%state%status = PF_STATUS_ITERATING
    if (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING .and. dir == 2) &
         pf%state%status = PF_STATUS_ITERATING
         
!     if (pf%state%status .ne. PF_STATUS_CONVERGED) 
    call pf_send_status(pf, 1+k, dir)
    call call_hooks(pf, 1, PF_POST_CONVERGENCE)

    ! XXX: this ain't so pretty, perhaps we should use the
    ! 'nmoved' thinger to break this cycle if everyone is
    ! done...

    if (pf%state%status == PF_STATUS_CONVERGED) then
       converged = .true.
       return
    end if

    if (0 == pf%comm%nproc) then
       pf%state%status = PF_STATUS_PREDICTOR
       converged = .true.
       return
    end if

  end subroutine pf_check_convergence_oc

    
  !>  Routine to do the pfasst iterations for optimal control problems on one block of processors until completion.
  !>  Each processor will do either a fixed number of iterations, or iterate until a tolerance is met
  !>  On calling, it is assumed that the levels are already loaded with the initial guesses
  !> 
  subroutine pf_pfasst_block_oc(pf, dt, nsteps, predict, flags, step)
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: nsteps 
    logical,           intent(in)    :: predict
    integer, optional, intent(in)    :: flags    !0 (default): sweep on y and p, 1: just y, 2: just p
    integer, optional, intent(in)    :: step
    ! not yet clear how to handle send and receive for forward and backward combined 

    type(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer                   :: j, k, l, which, pred_flags(1), dir !dir to choose forward or backward send 
    real(pfdp)                :: residual, energy
    real(pfdp), pointer       :: z(:)

    logical :: converged, qbroadcast
    logical :: did_post_step_hook

    call start_timer(pf, TTOTAL)

    which = 1
    if (present(flags)) which = flags
    ! send forward by default, even if sweeping on both components; send backwards if sweeping on p only
    dir = 1
    if(which == 2) dir = 2
    pred_flags(1) = which
    
    if( present(step) ) then
      pf%state%step    = step
    else
      pf%state%step    = pf%rank
    end if
      
    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
    pf%state%cycle   = -1
!     pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66
    pf%state%nsteps  = nsteps
    

    residual = -1
    energy   = -1
    did_post_step_hook = .false.
    
!    call pf%results%initialize(nsteps, pf%niters, pf%comm%nproc, pf%nlevels)

    
    do k = 1, 666666666 

       qbroadcast = .false.

       if (pf%state%status == PF_STATUS_CONVERGED .and. .not. did_post_step_hook) then
         call call_hooks(pf, -1, PF_POST_STEP)
         did_post_step_hook = .true.
         pf%state%itcnt = pf%state%itcnt + pf%state%iter
         pf%state%mysteps = pf%state%mysteps + 1
         exit
       end if

       ! jump to next block if we've reached the max iteration count
       if (pf%state%iter >= pf%niters) then
!           print *, pf%rank, 'pf%state%iter >= pf%niters'
          if (.not. did_post_step_hook) then
            call call_hooks(pf, -1, PF_POST_STEP)
            pf%state%itcnt = pf%state%itcnt + pf%state%iter
            pf%state%mysteps = pf%state%mysteps + 1
          end if
          did_post_step_hook = .false.

          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt
          
          if (pf%state%step >= pf%state%nsteps) exit  ! for optimal control this exit should always happen
          
          pf%state%status = PF_STATUS_PREDICTOR
          !pf%state%block  = pf%state%block + 1
          residual = -1
          qbroadcast = .true.
       end if

       if (k > 1 .and. qbroadcast) then
          if (pf%comm%nproc > 1) then
             stop "broadcast not supported" 
             !fine_lev_p => pf%levels(pf%nlevels)
             !call pf%comm%wait(pf, pf%nlevels)
             !call fine_lev_p%encap%pack(fine_lev_p%send, fine_lev_p%qend)
             !call pf_broadcast(pf, fine_lev_p%send, fine_lev_p%nvars, pf%comm%nproc-1)
             !call fine_lev_p%encap%unpack(fine_lev_p%q0,fine_lev_p%send)
          else
             stop "we should not be here I guess"
             ! for sequential optimal control, we need to save the Q(m) values for state solution
             ! and load them when solving the adjoint
             ! additionally, state solution is needed for objective, adjoint for gradient
             
             !print *, 'copying initial/terminal value'
             fine_lev_p => pf%levels(pf%nlevels)
             if ((which .eq. 0) .or. (which .eq. 1)) call fine_lev_p%q0%copy(fine_lev_p%qend, 1)
             if (which .eq. 2) call fine_lev_p%qend%copy(fine_lev_p%q0, 2)
          end if
       end if

      if (pf%state%status == PF_STATUS_PREDICTOR) then
        !print *, 'pf%state%status == PF_STATUS_PREDICTOR', pf%state%t0, dt, which
        if (predict) then
          !print *, 'calling predictor'
           call pf_predictor_oc(pf, pf%state%t0, dt, pred_flags)
        else
           pf%state%iter = 0
           pf%state%status  = PF_STATUS_ITERATING
           pf%state%pstatus = PF_STATUS_ITERATING
        end if
      end if    

      pf%state%iter  = pf%state%iter + 1

!       exit! just do predictor
      
      call start_timer(pf, TITERATION)
      call call_hooks(pf, -1, PF_PRE_ITERATION)
      
      if (pf%state%status /= PF_STATUS_CONVERGED) then
          fine_lev_p => pf%levels(pf%nlevels)
          call fine_lev_p%ulevel%sweeper%sweep(pf, pf%nlevels, pf%state%t0, dt, fine_lev_p%nsweeps, which)
       end if
      
      ! check convergence  (should always be not converged)
      call pf_check_convergence_oc(pf, k, dt, residual, energy, converged, dir)
      
      if (pf%state%step >= pf%state%nsteps) exit
      
      if (.not. converged) then
        !   non-blocking receive at all but the coarsest level
        do l = 2, pf%nlevels
          fine_lev_p => pf%levels(l)
          call pf_post(pf, fine_lev_p, fine_lev_p%index*10000+k, dir)
        end do

        if (pf%state%status /= PF_STATUS_CONVERGED) then
          fine_lev_p => pf%levels(pf%nlevels)
          call pf_send(pf, fine_lev_p, fine_lev_p%index*10000+k, .false., dir)
          if (pf%nlevels > 1) then
            coarse_lev_p => pf%levels(pf%nlevels-1)
            call restrict_time_space_fas(pf, pf%state%t0, dt, pf%nlevels, which)
            call save(coarse_lev_p, which)
          end if             
        end if
        
        call pf_v_cycle_oc(pf, k, pf%state%t0, dt, which)
        call call_hooks(pf, -1, PF_POST_ITERATION)
        call end_timer(pf, TITERATION)
      end if
    end do  !  Niter loop

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)
    
  end subroutine pf_pfasst_block_oc

  subroutine pf_v_cycle_oc(pf, iteration, t0, dt, flags)
  ! Execute a V-cycle between levels nfine and ncoarse

    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration
    integer, optional, intent(in)    :: flags

    type(pf_level_t), pointer :: f_lev_p, c_lev_p
    integer :: level_index, j, which, dir

    which = 1
    if(present(flags)) which = flags
    ! send forward by default, even if sweeping on both components; send backwards if sweeping on p only
    dir = 1
    if(which == 2) dir = 2 !

    !  For a single level, just get new initial conditions and return
    if (pf%nlevels == 1) then
       f_lev_p => pf%levels(1)
       call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration, .true., dir)
       return
    end if
    
    !
    ! down (fine to coarse)
    !
    do level_index = pf%nlevels-1, 2, -1
      f_lev_p => pf%levels(level_index);
      c_lev_p => pf%levels(level_index-1)
      call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps, which)
      call pf_send(pf, f_lev_p, level_index*10000+iteration, .false., dir)
      call restrict_time_space_fas(pf, t0, dt, level_index, which)
      call save(c_lev_p, which)
    end do

    !
    ! bottom  (coarsest level)
    !
    level_index=1
    f_lev_p => pf%levels(level_index)
    if (pf%pipeline_pred) then
       do j = 1, f_lev_p%nsweeps
          call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration+j, .true., dir)
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1, which)
          call pf_send(pf, f_lev_p, f_lev_p%index*10000+iteration+j, .false., dir)
       end do
    else
!       if (which == 0) then
!         call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration, .true., dir)
!         call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps, 1)
!         call pf_send(pf, f_lev_p, level_index*10000+iteration, .false., dir)
!         call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps, 2) ! this interferes with skipping y sweeps: have to check 
!                                                                                        ! state residual in case of which==1 in sweeper as well
!       else
        call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration, .true., dir)
        call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps, which)
        call pf_send(pf, f_lev_p, level_index*10000+iteration, .false., dir)
!       endif
    endif
    
    !
    ! up  (coarse to fine)
    !
    do level_index = 2, pf%nlevels
      f_lev_p => pf%levels(level_index);
      c_lev_p => pf%levels(level_index-1)
      call interpolate_time_space(pf, t0, dt, level_index, c_lev_p%Finterp, which)
      call pf_recv(pf, f_lev_p, level_index*10000+iteration, .false., dir)

       if (pf%rank /= 0) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          if ((which .eq. 0) .or. (which .eq. 1)) call interpolate_q0(pf, f_lev_p, c_lev_p, 1)
       end if
       if (pf%rank /= pf%comm%nproc-1) then
          if (which .eq. 2) call interpolate_qend(pf, f_lev_p, c_lev_p)
       end if

       if (level_index < pf%nlevels) then
          call call_hooks(pf, level_index, PF_PRE_SWEEP)
          ! compute residual
          ! do while residual > tol and j < nswps
          ! assuming residual computed at end of sweep 
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps, which)
       end if
    end do

  end subroutine pf_v_cycle_oc

end module pf_mod_parallel_oc
