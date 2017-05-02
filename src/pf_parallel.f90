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

! Parallel PFASST routines.

module pf_mod_parallel
  use pf_mod_dtype
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_hooks
  implicit none
contains

  !
  ! Predictor.
  !
  !
  !
  !>  This is the routine to initialize the solution on each processor
  !>
  !> Spreads the fine initial condition (F%q0) to all levels and all
  !> nodes.  If we're running with more than one processor, performs
  !> sweeps on the coarsest level.
  !> This routine is in need of reorganizing, since we now have cases
  !> where we already have a solution that we just want to load on levels.
  !>
  !> No time communication is performed during the predictor since all
  !> procesors can do the work themselves
  !>
  !>  The iteration count is reset to 0, and the status is reset to
  !>  ITERATING.
  !>
  subroutine pf_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf     !< PFASST main data structure
    real(pfdp),        intent(in   )         :: t0     !< Initial time of this processor
    real(pfdp),        intent(in   )         :: dt     !< time step

    class(pf_level_t), pointer :: coarse_lev_p
    class(pf_level_t), pointer :: fine_lev_p
    integer                   :: j, k
    integer                   :: level_index
    real(pfdp)                :: t0k

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    fine_lev_p => pf%levels(pf%nlevels)
    call spreadq0(fine_lev_p, t0)

    if (pf%nlevels > 1) then

       do level_index = pf%nlevels, 2, -1
         fine_lev_p => pf%levels(level_index);
         coarse_lev_p => pf%levels(level_index-1)
         call pf_residual(pf, fine_lev_p, dt)
         call restrict_time_space_fas(pf, t0, dt, level_index)
         call save(coarse_lev_p)
         call coarse_lev_p%q0%copy(coarse_lev_p%Q(1))
       end do

       if (pf%comm%nproc > 1) then
          coarse_lev_p => pf%levels(1)
          if (pf%Pipeline_G .and. (coarse_lev_p%nsweeps_pred > 1)) then
             !  This is the weird choice.  We burn in without communication, then do extra sweeps
             coarse_lev_p => pf%levels(1)
             do k = 1, pf%rank + 1
                pf%state%iter = -k

                ! Get new initial value (skip on first iteration)
                if (k > 1) then
                   call coarse_lev_p%q0%copy(coarse_lev_p%qend)
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(coarse_lev_p, t0)
                   end if
                end if

                call call_hooks(pf, coarse_lev_p%level, PF_PRE_SWEEP)
                call coarse_lev_p%ulevel%sweeper%sweep(pf, coarse_lev_p, t0, dt)
                call pf_residual(pf, coarse_lev_p, dt)  !  why is this here?
                call call_hooks(pf, coarse_lev_p%level, PF_POST_SWEEP)
             end do
             ! Now we have mimicked the burn in and we must do pipe-lined sweeps
             do k = 1, coarse_lev_p%nsweeps_pred-1
                pf%state%pstatus = PF_STATUS_ITERATING
                pf%state%status = PF_STATUS_ITERATING
                pf%state%iter =-(pf%rank + 1) -k

                !  Get new initial conditions
                call pf_recv(pf, coarse_lev_p, coarse_lev_p%level*20000+pf%rank+k, .true.)

                !  Do a sweep
                call call_hooks(pf, coarse_lev_p%level, PF_PRE_SWEEP)
                call coarse_lev_p%ulevel%sweeper%sweep(pf, coarse_lev_p, t0, dt )
                call pf_residual(pf, coarse_lev_p, dt)  !  why is this here?
                call call_hooks(pf, coarse_lev_p%level, PF_POST_SWEEP)
                !  Send forward
                call pf_send(pf, coarse_lev_p,  coarse_lev_p%level*20000+pf%rank+1+k, .false.)

             end do
             call pf_residual(pf, coarse_lev_p, dt)
          else
             ! Normal predictor burn in
             coarse_lev_p => pf%levels(1)
             do k = 1, pf%rank + 1
                pf%state%iter = -k
                t0k = t0-(pf%rank)*dt + (k-1)*dt

                ! Get new initial value (skip on first iteration)
                if (k > 1) then
                   call coarse_lev_p%q0%copy(coarse_lev_p%qend)
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(coarse_lev_p, t0k)
                   end if
                end if

                call call_hooks(pf, coarse_lev_p%level, PF_PRE_SWEEP)
                do j = 1, coarse_lev_p%nsweeps_pred
                   call coarse_lev_p%ulevel%sweeper%sweep(pf, coarse_lev_p, t0k, dt)
                end do
                call pf_residual(pf, coarse_lev_p, dt)
                call call_hooks(pf, coarse_lev_p%level, PF_POST_SWEEP)
             end do
          end if

          ! Return to fine level...
          call pf_v_cycle_post_predictor(pf, t0, dt)

       else

          ! Single processor... sweep on coarse and return to fine level.

          coarse_lev_p => pf%levels(1)
          do k = 1, pf%rank + 1
             pf%state%iter = -k
             t0k = t0-(pf%rank)*dt + (k-1)*dt

             call call_hooks(pf, coarse_lev_p%level, PF_PRE_SWEEP)
             do j = 1, coarse_lev_p%nsweeps_pred
                call coarse_lev_p%ulevel%sweeper%sweep(pf, coarse_lev_p, t0k, dt)
                call call_hooks(pf, coarse_lev_p%level, PF_POST_SWEEP)
             end do
             call pf_residual(pf, coarse_lev_p, dt)

          end do

          ! Return to fine level...
          call pf_v_cycle_post_predictor(pf, t0, dt)

       end if

    end if

    call end_timer(pf, TPREDICTOR)
    call call_hooks(pf, -1, PF_POST_PREDICTOR)

    pf%state%iter   = 0
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING

  end subroutine pf_predictor

  subroutine pf_check_tolerances(pf, residual, energy)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual, energy

    real(pfdp) :: residual1

    residual1 = pf%levels(pf%nlevels)%residual
    if (pf%state%status == PF_STATUS_ITERATING .and. residual > 0.0d0) then
       if ( (abs(1.0_pfdp - abs(residual1/residual)) < pf%rel_res_tol) .or. &
            (abs(residual1)                          < pf%abs_res_tol) ) then
          pf%state%status = PF_STATUS_CONVERGED
       end if
    end if
    residual = residual1

    pf%state%res = residual
  end subroutine pf_check_tolerances

  !
  ! Test residuals to determine if the current processor has converged.
  !
  ! Note that if the previous processor hasn't converged yet
  ! (pstatus), the current processor hasn't converged yet either,
  ! regardless of the residual.
  !
  subroutine pf_check_convergence(pf, k, dt, residual, energy, qcycle)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: residual, energy
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: k
    logical,           intent(out)   :: qcycle

    qcycle = .false.

    ! shortcut for fixed block mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    call pf_check_tolerances(pf, residual, energy)

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

  end subroutine pf_check_convergence

  !
  ! Run in parallel using PFASST.
  !
  subroutine pf_pfasst_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt, tend
    class(pf_encap_t), intent(inout), optional :: qend
    integer,           intent(in   ), optional :: nsteps

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer                   :: j, k
    integer                   :: level_index
    real(pfdp)                :: residual, energy

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

    fine_lev_p => pf%levels(pf%nlevels)
    call fine_lev_p%q0%copy(q0)

    if (present(nsteps)) then
       pf%state%nsteps = nsteps
    else
       pf%state%nsteps = ceiling(1.0*tend/dt)
    end if
    do k = 1, 666666666

       qbroadcast = .false.

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

       if (k > 1 .and. qbroadcast) then
          fine_lev_p => pf%levels(pf%nlevels)
          call pf%comm%wait(pf, pf%nlevels)
          call fine_lev_p%qend%pack(fine_lev_p%send)
          call pf_broadcast(pf, fine_lev_p%send, fine_lev_p%nvars, pf%comm%nproc-1)
          call fine_lev_p%q0%unpack(fine_lev_p%send)
       end if

       ! predictor, if requested
       if (pf%state%status == PF_STATUS_PREDICTOR) &
            call pf_predictor(pf, pf%state%t0, dt)

       !
       ! perform fine sweeps
       !

       pf%state%iter  = pf%state%iter + 1
       pf%state%cycle = 1

       call start_timer(pf, TITERATION)
       call call_hooks(pf, -1, PF_PRE_ITERATION)

       ! XXX: this if statement is necessary for block mode cycling...
       if (pf%state%status /= PF_STATUS_CONVERGED) then

          fine_lev_p => pf%levels(pf%nlevels)
          call call_hooks(pf, fine_lev_p%level, PF_PRE_SWEEP)
          do j = 1, fine_lev_p%nsweeps_pred
             call fine_lev_p%ulevel%sweeper%sweep(pf, fine_lev_p, pf%state%t0, dt)

             call pf_residual(pf, fine_lev_p, dt)
             call call_hooks(pf, fine_lev_p%level, PF_POST_SWEEP)
          end do
       end if

       !
       ! check convergence, continue with iteration
       !

       call pf_check_convergence(pf, k, dt, residual, energy, qcycle)

       if (pf%state%step >= pf%state%nsteps) exit

       if (qcycle) cycle
       do level_index = 2, pf%nlevels
          fine_lev_p => pf%levels(level_index)
          call pf_post(pf, fine_lev_p, fine_lev_p%level*10000+k)
       end do

       if (pf%state%status /= PF_STATUS_CONVERGED) then

          fine_lev_p => pf%levels(pf%nlevels)
          call pf_send(pf, fine_lev_p, fine_lev_p%level*10000+k, .false.)

          if (pf%nlevels > 1) then
             coarse_lev_p => pf%levels(pf%nlevels-1)
             call restrict_time_space_fas(pf, pf%state%t0, dt, pf%nlevels)
             call save(coarse_lev_p)
          end if

       end if

       call pf_v_cycle(pf, k, pf%state%t0, dt)
       call call_hooks(pf, -1, PF_POST_ITERATION)
       call end_timer(pf, TITERATION)

    end do

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)

    if (present(qend)) then
       fine_lev_p => pf%levels(pf%nlevels)
       call qend%copy(fine_lev_p%qend)
    end if
  end subroutine pf_pfasst_run

  !
  ! After predictor, return to fine level.
  !
  subroutine pf_v_cycle_post_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer :: j
    integer :: level_index

    if (pf%nlevels <= 1) return

    do level_index = 2, pf%nlevels-1
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
       call interpolate_time_space(pf, t0, dt, level_index, coarse_lev_p%Finterp)
       !       call fine_lev_p%Q(1)%pack(fine_lev_p%q0)
       call fine_lev_p%q0%copy(fine_lev_p%Q(1))
       call call_hooks(pf, fine_lev_p%level, PF_PRE_SWEEP)
       do j = 1, fine_lev_p%nsweeps_pred
          call fine_lev_p%ulevel%sweeper%sweep(pf, fine_lev_p, t0, dt)
          call pf_residual(pf, fine_lev_p, dt)
          call call_hooks(pf, fine_lev_p%level, PF_POST_SWEEP)
       end do

    end do
    coarse_lev_p => pf%levels(pf%nlevels-1)
    fine_lev_p => pf%levels(pf%nlevels)
    call interpolate_time_space(pf, t0, dt, pf%nlevels, coarse_lev_p%Finterp)
    call fine_lev_p%q0%copy(fine_lev_p%Q(1))

  end subroutine pf_v_cycle_post_predictor

  !
  ! Execute a V-cycle, starting and ending from the middle level.
  !
  subroutine pf_v_cycle(pf, iteration, t0, dt)
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration

    class(pf_level_t), pointer :: fine_lev_p, coarse_lev_p
    integer :: j
    integer :: level_index

    if (pf%nlevels == 1) then
       fine_lev_p => pf%levels(1)
       call pf_recv(pf, fine_lev_p, fine_lev_p%level*10000+iteration, .true.)
       return
    end if

    !
    ! down
    !
    do level_index = pf%nlevels-1, 2, -1
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
       call call_hooks(pf, level_index, PF_PRE_SWEEP)
       do j = 1, fine_lev_p%nsweeps
          call fine_lev_p%ulevel%sweeper%sweep(pf, fine_lev_p, t0, dt)
          call pf_residual(pf, fine_lev_p, dt)
          call call_hooks(pf, level_index, PF_POST_SWEEP)
       end do
       call pf_send(pf, fine_lev_p, level_index*10000+iteration, .false.)
       call restrict_time_space_fas(pf, t0, dt, level_index)
       call save(coarse_lev_p)
    end do

    !
    ! bottom
    !
    level_index=1
    fine_lev_p => pf%levels(level_index)
    if (pf%Pipeline_G) then
       do j = 1, fine_lev_p%nsweeps
          call pf_recv(pf, fine_lev_p, fine_lev_p%level*10000+iteration+j, .true.)
          call call_hooks(pf, level_index, PF_PRE_SWEEP)
          call fine_lev_p%ulevel%sweeper%sweep(pf, fine_lev_p, t0, dt)
          call pf_residual(pf, fine_lev_p, dt)
          call call_hooks(pf, level_index, PF_POST_SWEEP)
          call pf_send(pf, fine_lev_p, fine_lev_p%level*10000+iteration+j, .false.)
       end do
    else
       call pf_recv(pf, fine_lev_p, fine_lev_p%level*10000+iteration, .true.)
       call call_hooks(pf, level_index, PF_PRE_SWEEP)
       do j = 1, fine_lev_p%nsweeps
          call fine_lev_p%ulevel%sweeper%sweep(pf, fine_lev_p, t0, dt)
       end do
       call pf_residual(pf, fine_lev_p, dt)
       call call_hooks(pf, 1, PF_POST_SWEEP)
       call pf_send(pf, fine_lev_p, level_index*10000+iteration, .false.)
    endif
    !
    ! up
    !
    do level_index = 2, pf%nlevels
      fine_lev_p => pf%levels(level_index);
      coarse_lev_p => pf%levels(level_index-1)
      call interpolate_time_space(pf, t0, dt, level_index,coarse_lev_p%Finterp)
      call pf_recv(pf, fine_lev_p, level_index*10000+iteration, .false.)

       if (pf%rank /= 0) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          call interpolate_q0(pf,fine_lev_p, coarse_lev_p)
       end if

       if (level_index < pf%nlevels) then
          call call_hooks(pf, level_index, PF_PRE_SWEEP)
          do j = 1, fine_lev_p%nsweeps
             call fine_lev_p%ulevel%sweeper%sweep(pf, fine_lev_p, t0, dt)
          end do
          call pf_residual(pf, fine_lev_p, dt)
          call call_hooks(pf, level_index, PF_POST_SWEEP)
       end if
    end do

  end subroutine pf_v_cycle

  !
  ! Communication helpers
  !
  subroutine pf_post(pf, level, tag)
    type(pf_pfasst_t), intent(in)    :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    if (pf%rank /= 0 .and. pf%state%pstatus == PF_STATUS_ITERATING) then
       call pf%comm%post(pf, level, tag)
    end if
  end subroutine pf_post

  subroutine pf_send_status(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    if (pf%rank /= pf%comm%nproc-1) then
       call pf%comm%send_status(pf, tag)
    end if
  end subroutine pf_send_status

  subroutine pf_recv_status(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    if (pf%rank /= 0) then
       call pf%comm%recv_status(pf, tag)
    end if
  end subroutine pf_recv_status

  subroutine pf_send(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    call start_timer(pf, TSEND + level%level - 1)
    if (pf%rank /= pf%comm%nproc-1 &
         .and. pf%state%status == PF_STATUS_ITERATING) then
       call pf%comm%send(pf, level, tag, blocking)
    end if
    call end_timer(pf, TSEND + level%level - 1)
  end subroutine pf_send

  subroutine pf_recv(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    call start_timer(pf, TRECEIVE + level%level - 1)
    if (pf%rank /= 0 .and.  pf%state%pstatus == PF_STATUS_ITERATING) then

       call pf%comm%recv(pf, level,tag, blocking)

    end if
    call end_timer(pf, TRECEIVE + level%level - 1)
  end subroutine pf_recv

  subroutine pf_broadcast(pf, y, nvar, root)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root
    call start_timer(pf, TBROADCAST)
    call pf%comm%broadcast(pf, y, nvar, root)
    call end_timer(pf, TBROADCAST)
  end subroutine pf_broadcast

end module pf_mod_parallel
