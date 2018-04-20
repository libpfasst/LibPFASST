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

!> Module of parallel PFASST routines.

module pf_mod_parallel
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_pfasst

  implicit none
contains

  !>  This is the main interface to pfasst.
  !!  It examines the parameters and decides which subroutine to call
  !!  to execute the code correctly
  subroutine pf_pfasst_run(pf, q0, dt, tend, nsteps, qend, flags)
    type(pf_pfasst_t), intent(inout), target   :: pf   !<  The complete PFASST structure
    class(pf_encap_t), intent(in   )           :: q0   !<  The initial condition
    real(pfdp),        intent(in   )           :: dt   !<  The time step for each processor
    real(pfdp),        intent(in   )           :: tend !<  The final time of run
    integer,           intent(in   ), optional :: nsteps  !<  The number of time steps
    class(pf_encap_t), intent(inout), optional :: qend    !<  The computed solution at tend
    integer,           intent(in   ), optional :: flags(:)!<  User defnined flags


    !  Local variables
    integer :: nproc  !<  Total number of processors
    integer :: nsteps_loc  !<  local number of time steps    
    real(pfdp) :: tend_loc !<  The final time of run

    
    
    ! make a local copy of nproc
    nproc = pf%comm%nproc

    !>  Set the number of time steps to do
    !!  The user can either pass in the number of time steps or
    !!  pass in the time step size and length of run
    if (present(nsteps)) then
      nsteps_loc = nsteps
      tend_loc=dble(nsteps_loc*dt)
    else
      nsteps_loc = ceiling(1.0*tend/dt)
      !  Do  sanity check on steps
      if (abs(dble(nsteps_loc)-tend/dt) > dt/100.0) then
        print *,'dt=',dt
        print *,'nsteps=',nsteps_loc
        print *,'tend=',tend
        stop "Invalid nsteps"
      end if
    end if
    pf%state%nsteps = nsteps_loc

    !  do sanity checks on Nproc
    if (mod(nsteps,nproc) > 0) stop "ERROR: nsteps must be multiple of nproc (pf_parallel.f90)."

    call pf%results%initialize(nsteps_loc, pf%niters, pf%comm%nproc, pf%nlevels)

    ! figure out what routine to call
    if (pf%nlevels .eq. 1) then
       print *,'Calling pipelined SDC with 1 level'
       if (present(qend)) then
          call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc, qend=qend,flags=flags)
       else
          call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc,flags=flags)
       end if
    else
       if (pf%Vcycle) then
          !  Right now, we just call the old routine
          if (present(qend)) then
             call pf_pfasst_run_old(pf, q0, dt, tend_loc, nsteps_loc, qend=qend,flags=flags)
          else
             call pf_pfasst_run_old(pf, q0, dt, tend_loc, nsteps_loc,flags=flags)
          end if
       else
          print *,'Calling pipelined SDC with multiple levels'
          if (present(qend)) then
             call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc, qend=qend,flags=flags)
          else
             call pf_pipeline_run(pf, q0, dt, tend_loc, nsteps_loc,flags=flags)
          end if
       end if
    end if

    if (pf%save_results) then
       call pf%results%dump()
    endif

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
    end if
    if (which == 2) then
      call end_timer(pf, TPREDICTOR)
      call call_hooks(pf, -1, PF_POST_PREDICTOR)

      pf%state%iter   = 0
      pf%state%status = PF_STATUS_ITERATING
      pf%state%pstatus = PF_STATUS_ITERATING
      return
    end if
        
    !! Step 3. Do the "Burn in" step on the coarse level to make the coarse values consistent
    !!         (this is skipped if the fine initial conditions are already consistent)
    !! The first processor does nothing, the second does one set of sweeps, the 2nd two, etc
    !! Hence, this is skipped completely if nprocs=1
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
             t0k = t0-(pf%rank)*dt + (k-1)*dt   !  Remember t0 is the beginning of this time slice so t0-(pf%rank)*dt is t0 of problem

             ! Get new initial value (skip on first iteration)
             if (k > 1) then
                if ((which == 0) .or. (which == 1)) call c_lev_p%q0%copy(c_lev_p%qend, 1)
                if ((which == 0) .or. (which == 2)) call c_lev_p%qend%copy(c_lev_p%q0, 2)
                ! If we are doing PFASST_pred, we use the old values at nodes, otherwise spread q0
                if (.not. pf%PFASST_pred) then
                   if( (which == 0) .or. (which == 1)) call c_lev_p%ulevel%sweeper%spreadq0(c_lev_p, t0k, 1, pf%state%step+1)
                   if( (which == 0) .or. (which == 2)) call c_lev_p%ulevel%sweeper%spreadq0(c_lev_p, t0k+dt, 2, pf%state%step+1)
                end if
             end if
             !  Do some sweeps
             call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0k, dt,pf%nsweeps_burn, which)
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
       call c_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, c_lev_p%nsweeps_pred, which)
       !  Send forward
       call pf_send(pf, c_lev_p,  c_lev_p%index*20000+pf%rank+1, .false., dir)
    endif  ! (Pipeline_pred .eq. .true) then

    !  Step 6:  Return to fine level sweeping on any level in between coarsest and finest
    do level_index = 2, pf%nlevels  !  Will do nothing with one level
       f_lev_p => pf%levels(level_index);
       c_lev_p => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, c_lev_p%Finterp, which)
       if ((which == 0) .or. (which == 1)) call interpolate_q0(pf, f_lev_p, c_lev_p)
       if ((which == 0) .or. (which == 2)) call interpolate_qend(pf, f_lev_p, c_lev_p)
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

  end subroutine pf_predictor

  !>
  !> Test residuals to determine if the current processor has converged.
  !>
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor hasn't converged yet either,
  !> regardless of the residual.
  !>
  subroutine pf_check_convergence_old(pf, k, dt, residual, energy, qcycle)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual, energy
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: k
    logical,           intent(out)   :: qcycle
    real(pfdp)     :: residual1
    qcycle = .false.

    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    ! Check to see if tolerances are met
    residual1 = pf%levels(pf%nlevels)%residual
    if (pf%state%status == PF_STATUS_ITERATING .and. residual > 0.0d0) then
       if ( (abs(1.0_pfdp - abs(residual1/residual)) < pf%rel_res_tol) .or. &
            (abs(residual1)                          < pf%abs_res_tol) ) then
          pf%state%status = PF_STATUS_CONVERGED
       end if
    end if
    residual = residual1

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)
    call pf_recv_status(pf, 8000+k)

    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING) &
         pf%state%status = PF_STATUS_ITERATING

    call pf_send_status(pf, 8000+k)
    call call_hooks(pf, 1, PF_POST_CONVERGENCE)

    ! XXX: this ain't so pretty, perhaps we should use the
    ! 'nmoved' thinger to break this cycle if everyone is
    ! done...

    if (pf%state%status == PF_STATUS_CONVERGED) then
       qcycle = .true.
       return
    end if

    if (0 == pf%comm%nproc) then
       pf%state%status = PF_STATUS_PREDICTOR
       qcycle = .true.
       return
    end if

  end subroutine pf_check_convergence_old

    subroutine pf_check_convergence(pf, k, dt, residual, energy, converged)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual, energy
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: k
    logical,           intent(out)   :: converged   !<  True if this processor is done
    real(pfdp)     :: residual1
    converged = .false.

    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    ! Check to see if tolerances are met
    residual1 = pf%levels(pf%nlevels)%residual
    if (pf%state%status == PF_STATUS_ITERATING .and. residual > 0.0d0) then
       if ( (abs(1.0_pfdp - abs(residual1/residual)) < pf%rel_res_tol) .or. &
            (abs(residual1)                          < pf%abs_res_tol) ) then
          pf%state%status = PF_STATUS_CONVERGED
       end if
    end if
    residual = residual1

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)
    call pf_recv_status(pf, 8000+k)

    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING) &
         pf%state%status = PF_STATUS_ITERATING

    call pf_send_status(pf, 8000+k)
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

  end subroutine pf_check_convergence

  !
  ! Run in parallel using PFASST.
  !
  subroutine pf_pfasst_run_old(pf, q0, dt, tend, nsteps, qend,flags)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt, tend
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend
    integer,           intent(in   ), optional :: flags(:)

    class(pf_level_t), pointer :: f_lev_p, c_lev_p
    integer                   :: j, k
    integer                   :: level_index
    real(pfdp)                :: residual, energy
    integer                   ::  ierror  !<  Warning flag for communication routines

    logical :: qcycle, qbroadcast
    logical :: did_post_step_hook

    call start_timer(pf, TTOTAL)


    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
    pf%state%cycle   = -1
    pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66

    residual = -1
    energy   = -1
    did_post_step_hook = .false.

    f_lev_p => pf%levels(pf%nlevels)
    call f_lev_p%q0%copy(q0, flags=1)

    do k = 1, 666666666   !  Loop over blocks of time steps

       qbroadcast = .false.
       
       !  Check to see if we should do one more hook 
       if (pf%state%status == PF_STATUS_CONVERGED .and. .not. did_post_step_hook) then
         call call_hooks(pf, -1, PF_POST_STEP)
         did_post_step_hook = .true.
         pf%state%itcnt = pf%state%itcnt + pf%state%iter-1
         pf%state%mysteps = pf%state%mysteps + 1
       end if

       ! in block mode, jump to next block if we've reached the max iteration count
       if (pf%state%iter >= pf%niters) then

          if (.not. did_post_step_hook) then
            call call_hooks(pf, -1, PF_POST_STEP)
            pf%state%itcnt = pf%state%itcnt + pf%state%iter-1
            pf%state%mysteps = pf%state%mysteps + 1
          end if
          did_post_step_hook = .false.

          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt

          if (pf%state%step >= pf%state%nsteps) exit

          pf%state%status = PF_STATUS_PREDICTOR
          qbroadcast = .true.
       end if

       !  Do this when starting a new block, broadcast new initial conditions to all procs
       if (k > 1 .and. qbroadcast) then
          f_lev_p => pf%levels(pf%nlevels)
          call pf%comm%wait(pf, pf%nlevels,ierror)             !<  make sure everyone is done
          call f_lev_p%qend%pack(f_lev_p%send)    !<  Pack away your last solution
          call pf_broadcast(pf, f_lev_p%send, f_lev_p%mpibuflen, pf%comm%nproc-1)
          call f_lev_p%q0%unpack(f_lev_p%send)    !<  Everyone resets their q0
       end if

       ! predictor, if requested or we are starting new bloc
       if (pf%state%status == PF_STATUS_PREDICTOR) &
            call pf_predictor(pf, pf%state%t0, dt,flags)

       !
       ! perform fine sweeps
       !
       pf%state%iter  = pf%state%iter + 1
       pf%state%cycle = 1
       
       call start_timer(pf, TITERATION)

       ! XXX: this if statement is necessary for block mode cycling...
       if (pf%state%status /= PF_STATUS_CONVERGED) then

          f_lev_p => pf%levels(pf%nlevels)
          call f_lev_p%ulevel%sweeper%sweep(pf, pf%nlevels, pf%state%t0, dt, f_lev_p%nsweeps)
       end if

       !
       ! check convergence, continue with iteration
       !

       !       call pf_check_convergence(pf, k, dt, residual, energy, converged)
       call pf_check_convergence_old(pf, k, dt, residual, energy, qcycle)


       if (pf%state%step >= pf%state%nsteps) exit
       if (qcycle) cycle

!       if (.not. converged) then
       do level_index = 2, pf%nlevels
          f_lev_p => pf%levels(level_index)
          call pf_post(pf, f_lev_p, f_lev_p%index*10000+k)
       end do

       if (pf%state%status /= PF_STATUS_CONVERGED) then

          f_lev_p => pf%levels(pf%nlevels)
          call pf_send(pf, f_lev_p, f_lev_p%index*10000+k, .false.)

          if (pf%nlevels > 1) then
             c_lev_p => pf%levels(pf%nlevels-1)
             call restrict_time_space_fas(pf, pf%state%t0, dt, pf%nlevels)
             call save(c_lev_p)
          end if

       end if

       call pf_v_cycle(pf, k, pf%state%t0, dt)
       call call_hooks(pf, -1, PF_POST_ITERATION)
       call end_timer(pf, TITERATION)

    end do  !   Loop on k over blocks of time steps

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       f_lev_p => pf%levels(pf%nlevels)
       call qend%copy(f_lev_p%qend, flags=1)
    end if
  end subroutine pf_pfasst_run_old

  !
  !> Run single level SDC in pipeline fashion
  !
  subroutine pf_pipeline_run(pf, q0, dt, tend, nsteps, qend,flags)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt, tend
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend
    integer,           intent(in   ), optional :: flags(:)
    
    class(pf_level_t), pointer :: lev_p  !<  pointer to the one level we are operating on
    integer                   :: j, k
    real(pfdp)                :: residual
    integer                   :: nblocks !<  The number of blocks of steps to do
    integer                   :: nproc   !<  The number of processors being used
    integer                   :: ierror  !<  Warning flag for communication routines

    logical :: converged   !<  True when this processor is converged to residual

    call start_timer(pf, TTOTAL)

    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
    pf%state%cycle   = -1
    pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66

!    if (pf%nlevels > 1) stop "ERROR: nlevels  must be 1 to run pipeline mode (pf_parallel.f90)"

    !  pointer to fine level on which we will iterate
    lev_p => pf%levels(pf%nlevels)
    call lev_p%q0%copy(q0, flags=1)

    nproc = pf%comm%nproc
    nblocks = nsteps/nproc
    do k = 1, nblocks   !  Loop over blocks of time steps
       ! print *,'Starting  step=',pf%state%step,'  block k=',k      
       ! Each block will consist of
       !  1.  predictor
       !  2.  A loop until max iterations, or tolerances met
       !      2a.  Recieve
       !      2b.  Do SDC sweep(s)1
       !      2c.  Send
       !  3.  Move solution to next block


       !>  When starting a new block, broadcast new initial conditions to all procs
       !>  For initial block, this is done when initial conditions are set
       if (k > 1) then
          if (nproc > 1)  then
             call lev_p%qend%pack(lev_p%send)    !<  Pack away your last solution
             call pf_broadcast(pf, lev_p%send, lev_p%mpibuflen, pf%comm%nproc-1)
             call lev_p%q0%unpack(lev_p%send)    !<  Everyone resets their q0
          else
             call lev_p%q0%copy(lev_p%qend, flags=1)    !<  Just stick qend in q0
          end if

          !>  Update the step and t0 variables for new block
          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt

          pf%state%status = PF_STATUS_PREDICTOR
          pf%state%iter    = -1
          pf%state%cycle   = -1
          pf%state%itcnt   = 0
          pf%state%mysteps = 0
          pf%state%status  = PF_STATUS_PREDICTOR
          pf%state%pstatus = PF_STATUS_PREDICTOR
          pf%comm%statreq  = -66
       end if

       !> Call the predictor
       !> Currently the predictor will do nothing but spread q0 to all the nodes
       if (pf%state%status == PF_STATUS_PREDICTOR) then
          call pf_predictor(pf, pf%state%t0, dt, flags)
       end if

       !>  Start the loops over SDC sweeps
       pf%state%iter = 0
       converged = .FALSE.
       pf%state%status = PF_STATUS_ITERATING
       pf%state%pstatus = PF_STATUS_ITERATING

       pf%results%times(pf%state%step+1, lev_p%index, pf%rank+1) = pf%state%t0

       call start_timer(pf, TITERATION)
       do j = 1, pf%niters
          call call_hooks(pf, -1, PF_PRE_ITERATION)

          pf%state%iter = j

          if (j>1 .and.  pf%state%pstatus /= PF_STATUS_CONVERGED) then
             call pf_recv_status(pf, 8000+k)
             call pf_recv(pf, lev_p, lev_p%index*10000+100*k+pf%state%iter, .true.)
          endif

          call lev_p%ulevel%sweeper%sweep(pf, pf%nlevels, pf%state%t0, dt, lev_p%nsweeps)
          call pf_check_convergence_pipeline(pf, lev_p%residual, converged)

          if (pf%state%status .ne. PF_STATUS_CONVERGED) then

             call pf_send(pf, lev_p, lev_p%index*10000+100*k+pf%state%iter, .false.)

             if (converged) then
                pf%state%status = PF_STATUS_CONVERGED
             endif
             call pf_send_status(pf, 8000+k)
          endif

          call call_hooks(pf, -1, PF_POST_ITERATION)

          pf%results%residuals(pf%state%iter, pf%state%step+1, lev_p%index) = lev_p%residual

          if (pf%state%status == PF_STATUS_CONVERGED) exit
       end do  !  Loop over the iteration in this bloc

       call end_timer(pf, TITERATION)

    end do

    call end_timer(pf, TTOTAL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       call qend%copy(lev_p%qend, flags=1)
    end if
  end subroutine pf_pipeline_run

  !>
  !> Test residuals to determine if the current processor has converged.
  !>
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor hasn't converged yet either,
  !> regardless of the residual.
  !>
  subroutine pf_check_convergence_pipeline(pf, residual, converged)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual
    logical,           intent(out)   :: converged   !<  True if this processor is done
    real(pfdp)     :: residual1

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)

    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    ! Check to see if tolerances are met
    if (pf%state%status == PF_STATUS_ITERATING .and. residual > 0.0d0) then
       if (pf%rank > 0) then
          if (abs(residual) < pf%abs_res_tol .and. &
               pf%state%pstatus == PF_STATUS_CONVERGED) then
             converged = .true.
          endif
       else
          if (abs(residual) < pf%abs_res_tol) converged = .true.
       endif
    end if

    call call_hooks(pf, 1, PF_POST_CONVERGENCE)

  end subroutine pf_check_convergence_pipeline
  
  
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
    call pf_recv_status(pf, 8000+k, dir)

    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING .and. dir == 1) &
         pf%state%status = PF_STATUS_ITERATING
    if (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING .and. dir == 2) &
         pf%state%status = PF_STATUS_ITERATING
         
    call pf_send_status(pf, 8000+k, dir)
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
    pf%state%itcnt   = 0
    pf%state%mysteps = 0
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66
    pf%state%nsteps  = nsteps
    

    residual = -1
    energy   = -1
    did_post_step_hook = .false.
    
    call pf%results%initialize(nsteps, pf%niters, pf%comm%nproc, pf%nlevels)

    
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
             stop "broadcast not supported yet" 
             !fine_lev_p => pf%levels(pf%nlevels)
             !call pf%comm%wait(pf, pf%nlevels)
             !call fine_lev_p%encap%pack(fine_lev_p%send, fine_lev_p%qend)
             !call pf_broadcast(pf, fine_lev_p%send, fine_lev_p%nvars, pf%comm%nproc-1)
             !call fine_lev_p%encap%unpack(fine_lev_p%q0,fine_lev_p%send)
          else
             ! for sequential optimal control, we need to save the Q(m) values for state solution
             ! and load them when solving the adjoint
             ! additionally, state solution is needed for objective, adjoint for gradient
             
             !print *, 'copying initial/terminal value'
             fine_lev_p => pf%levels(pf%nlevels)
             if ((which .eq. 0) .or. (which .eq. 1)) call fine_lev_p%q0%copy(fine_lev_p%qend, 1)
             if ((which .eq. 0) .or. (which .eq. 2)) call fine_lev_p%qend%copy(fine_lev_p%q0, 2)
          end if
       end if

      if (pf%state%status == PF_STATUS_PREDICTOR) then
        !print *, 'pf%state%status == PF_STATUS_PREDICTOR', pf%state%t0, dt, which
        if (predict) then
          !print *, 'calling predictor'
           call pf_predictor(pf, pf%state%t0, dt, pred_flags)
        else
           pf%state%iter = 0
           pf%state%status  = PF_STATUS_ITERATING
           pf%state%pstatus = PF_STATUS_ITERATING
        end if
      end if    

      pf%state%iter  = pf%state%iter + 1

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
        
        call pf_v_cycle(pf, k, pf%state%t0, dt, which)
        call call_hooks(pf, -1, PF_POST_ITERATION)
        call end_timer(pf, TITERATION)
      end if
    end do  !  Niter loop

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)
    
  end subroutine pf_pfasst_block_oc

  

  subroutine pf_v_cycle(pf, iteration, t0, dt, flags)
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
    if (pf%Pipeline_G) then
       do j = 1, f_lev_p%nsweeps
          call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration+j, .true., dir)
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, 1, which)
          call pf_send(pf, f_lev_p, f_lev_p%index*10000+iteration+j, .false., dir)
       end do
    else
       call pf_recv(pf, f_lev_p, f_lev_p%index*10000+iteration, .true., dir)
       call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps, which)
       call pf_send(pf, f_lev_p, level_index*10000+iteration, .false., dir)
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
          if ((which .eq. 0) .or. (which .eq. 1)) call interpolate_q0(pf, f_lev_p, c_lev_p)
       end if
       if (pf%rank /= pf%comm%nproc-1) then
          if ((which .eq. 0) .or. (which .eq. 2)) call interpolate_qend(pf, f_lev_p, c_lev_p)
       end if

       if (level_index < pf%nlevels) then
          call call_hooks(pf, level_index, PF_PRE_SWEEP)
          ! compute residual
          ! do while residual > tol and j < nswps
          ! assuming residual computed at end of sweep 
          call f_lev_p%ulevel%sweeper%sweep(pf, level_index, t0, dt, f_lev_p%nsweeps, which)
       end if
    end do

  end subroutine pf_v_cycle


    !
  !> Communication helpers
  !
  !>  Subroutine to post a receive request for a new initial condition to be received after doing some work
  subroutine pf_post(pf, level, tag, direction)
    type(pf_pfasst_t), intent(in)    :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: direction
    integer                          :: dir
    integer ::  ierror 
    
    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction
    
    ierror = 0
    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                  .and. dir == 1) then
       call pf%comm%post(pf, level, tag, ierror, dir)
    elseif (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                  .and. dir == 2) then
       call pf%comm%post(pf, level, tag, ierror, dir)
    end if
    
    if (ierror /= 0) then
      print *, pf%rank, 'warning: error during post', ierror
      stop "pf_parallel:pf_post"
    endif
  end subroutine pf_post
  
  !>  Subroutine to send this processor's convergence status to the next processor
  subroutine pf_send_status(pf, tag, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: direction
    integer ::  dir
    integer ::  istatus
    integer ::  ierror
    
    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction

    ierror = 0
    istatus = pf%state%status
    if (pf%rank /= pf%comm%nproc-1 .and. dir == 1) then
       if (pf%debug) print*, pf%rank, 'is sending status', pf%state%status, 'with tag =', tag 
       call pf%comm%send_status(pf, tag, istatus, ierror, dir)
       if (pf%debug) print*, pf%rank, 'status sent' 
    elseif (pf%rank /= 0 .and. dir == 2) then
       if (pf%debug) print*, pf%rank, 'is sending status', pf%state%status, 'backwards with tag =', tag 
       call pf%comm%send_status(pf, tag, istatus, ierror, dir)
       if (pf%debug) print*, pf%rank, 'status sent' 
    end if
    
    if (ierror /= 0) then
      print *, pf%rank, 'warning: error during send_status', ierror
      stop "pf_parallel:pf_send_status"
    endif
    
  end subroutine pf_send_status
  
  !>  Subroutine to receive the convergence status from the previous processor
  subroutine pf_recv_status(pf, tag, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    integer, optional, intent(in)    :: direction
    integer ::  dir
    integer ::  ierror, istatus
    
    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction
    
    ierror = 0
!     if (pf%rank /= 0 .and. dir == 1) then
   if (pf%rank /= 0 .and. pf%state%pstatus .ne. PF_STATUS_CONVERGED .and. dir == 1) then
       if (pf%debug) print*, pf%rank, 'my status = ', pf%state%status
       if (pf%debug) print*, pf%rank,  'is receiving status with tag ', tag  
       call pf%comm%recv_status(pf, tag, istatus, ierror, dir)
       if (pf%debug) print *, pf%rank, 'status recvd = ', istatus 
       if (ierror .eq. 0) pf%state%pstatus = istatus
!     elseif (pf%rank /= pf%comm%nproc-1 .and. dir == 2) then
   elseif (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus .ne. PF_STATUS_CONVERGED .and. dir == 2) then
       if (pf%debug) print*, pf%rank, 'my status = ', pf%state%status
       if (pf%debug) print*, pf%rank,  'is receiving status backwards with tag ', tag  
       call pf%comm%recv_status(pf, tag, istatus, ierror, dir)
       if (pf%debug) print *, pf%rank, 'status recvd = ', istatus 
       if (ierror .eq. 0) pf%state%pstatus = istatus
    end if
       
    if (ierror .ne. 0) then
      print *, pf%rank, 'warning: error during recv_status', ierror
      stop "pf_parallel_oc:pf_recv_status"
    endif    
  end subroutine pf_recv_status

  !>  Subroutine to send the solution to the next processor
  subroutine pf_send(pf, level, tag, blocking, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    integer, optional, intent(in)    :: direction
    integer                          :: dir, ierror
    
    dir = 1 ! default: send forward
    if(present(direction)) dir = direction

    ierror = 0
    call start_timer(pf, TSEND + level%index - 1)
    if (pf%rank /= pf%comm%nproc-1 &
         .and. pf%state%status == PF_STATUS_ITERATING & 
         .and. dir == 1 ) then
       ! print*, pf%rank,  'is sending soln'
       call pf%comm%send(pf, level, tag, blocking, ierror, dir)
    elseif (pf%rank /= 0 &
         .and. pf%state%status == PF_STATUS_ITERATING & 
         .and. dir == 2 ) then
       !print *, pf%rank, 'sending backward',tag,blocking,pf%rank
       call pf%comm%send(pf, level, tag, blocking, ierror, dir)
    end if
    if (ierror /= 0) then
      print *, pf%rank, 'warning: error during send', ierror
      stop "pf_parallel:pf_send"
    endif
    call end_timer(pf, TSEND + level%index - 1)
  end subroutine pf_send
  
  !>  Subroutine to recieve the solution from the previous processor
  subroutine pf_recv(pf, level, tag, blocking, direction)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    integer, optional, intent(in)    :: direction
    integer                          :: dir, ierror
    
    dir = 1 ! default: send forward
    if(present(direction)) dir = direction

    ierror = 0
    call start_timer(pf, TRECEIVE + level%index - 1)
    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                  .and. dir == 1) then
       !print *,'recv forward',tag,blocking,pf%rank
       call pf%comm%recv(pf, level,tag, blocking, ierror, dir)
       !print *,'done recv',tag,blocking,pf%rank
       if (ierror .eq. 0) call level%q0%unpack(level%recv, 1)
    elseif (pf%rank /= pf%comm%nproc-1 .and. pf%state%pstatus == PF_STATUS_ITERATING &
                                     .and. dir == 2) then
       !print *, pf%rank, 'recv backward',tag,blocking,pf%rank
       call pf%comm%recv(pf, level, tag, blocking, ierror, dir)
       !print *, pf%rank, 'done recv',tag,blocking,pf%rank
       if (ierror .eq. 0) call level%qend%unpack(level%recv, 2)
    end if
    
    if(ierror .ne. 0) then
      print *, pf%rank, 'warning: mpi error during receive', ierror
      stop "pf_parallel:pf_recv"     
    end if
    
    call end_timer(pf, TRECEIVE + level%index - 1)
  end subroutine pf_recv


  !>  Subroutine to broadcast the initial condition to all processors
  subroutine pf_broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root
    integer :: ierror
    call start_timer(pf, TBROADCAST)
    call pf%comm%broadcast(pf, y, nvar, root, ierror)
       if (ierror /= 0) then
          print *, pf%rank, 'warning:  error during broadcast', ierror
          stop "pf_parallel:pf_broadcast"
       endif
    call end_timer(pf, TBROADCAST)
  end subroutine pf_broadcast

    !> Save current solution and function value so that future corrections can be computed
  subroutine save(lev, flags)
    class(pf_level_t), intent(inout) :: lev  !<  Level to save on
    integer, optional, intent(in)   :: flags !<  which component to save (state/adjoint)
    integer :: m, p, which
    
    which = 1
    if(present(flags)) which = flags
    
    if (lev%Finterp) then
       if (allocated(lev%pFflt)) then
          do m = 1, lev%nnodes
             do p = 1,size(lev%F(1,:))
                call lev%pF(m,p)%copy(lev%F(m,p), which)
             end do
             call lev%pQ(m)%copy(lev%Q(m), which)
          end do
       end if
    else
       if (allocated(lev%pQ)) then
          do m = 1, lev%nnodes
             call lev%pQ(m)%copy(lev%Q(m), which)
          end do
       end if
    end if
  end subroutine save
  
end module pf_mod_parallel
