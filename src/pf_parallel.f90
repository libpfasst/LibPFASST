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
  ! Spreads the fine initial condition (F%q0) to all levels and all
  ! nodes.  If we're running with more than one processor, performs
  ! sweeps on the coarse accordingly and execute "start cycles".
  !
  ! No time communication is performed during the predictor.
  !
  ! The iteration count is reset to 0, and the status is reset to
  ! ITERATING.
  !
  subroutine pf_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: t0, dt

    type(pf_level_t), pointer :: F, G
    integer                   :: j, k, l
    real(pfdp)                :: t0k

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    F => pf%levels(pf%nlevels)
    call spreadq0(F, t0)

    do l = pf%nlevels, 2, -1
       F => pf%levels(l); G => pf%levels(l-1)
       call pf_residual(pf, F, dt)
       call restrict_time_space_fas(pf, t0, dt, F, G)
       call save(G)
       call G%encap%pack(G%q0, G%Q(1))
    end do


    if (pf%comm%nproc > 1) then
       if (pf%Pipeline_G .and. (G%nsweeps > 1)) then
          !  This is the weird choice.  We burn in without communication, then do extra sweeps
          G => pf%levels(1)
          do k = 1, pf%rank + 1
             pf%state%iter = -k

             ! get new initial value (skip on first iteration)
             if (k > 1) then
!                call G%encap%pack(G%q0, G%Q(G%nnodes))
                call G%encap%pack(G%q0, G%qend)
                if (.not. pf%PFASST_pred) then
                   call spreadq0(G, t0)
                end if
             end if

             call call_hooks(pf, G%level, PF_PRE_SWEEP)
             call G%sweeper%sweep(pf, G, t0, dt)
             call call_hooks(pf, G%level, PF_POST_SWEEP)
             call pf_residual(pf, G, dt)  !  why is this here?
          end do
          ! Now we have mimicked the burn in and we must do pipe-lined sweeps
          do k = 1, G%nsweeps-1
             pf%state%iter =-(pf%rank + 1) -k

             !  Get new initial conditions
             call pf_recv(pf, G, G%level*20000+pf%rank, .true.)
             !  Do a sweep
             call call_hooks(pf, G%level, PF_PRE_SWEEP)
             call G%sweeper%sweep(pf, G, t0, dt)
             call call_hooks(pf, G%level, PF_POST_SWEEP)
             !  Send forward
             call pf_send(pf, G,  G%level*20000+pf%rank+1, .true.)
          end do
          call pf_residual(pf, G, dt)  !  why is this here?
       else
          ! Normal predictor burn in
          G => pf%levels(1)
          do k = 1, pf%rank + 1
             pf%state%iter = -k
             t0k = t0-(pf%rank)*dt + (k-1)*dt

             ! get new initial value (skip on first iteration)
             if (k > 1) then
!                call G%encap%pack(G%q0, G%Q(G%nnodes))
                call G%encap%pack(G%q0, G%qend)
                if (.not. pf%PFASST_pred) then
                   call spreadq0(G, t0k)
                end if
             end if

             call call_hooks(pf, G%level, PF_PRE_SWEEP)
             do j = 1, G%nsweeps
                call G%sweeper%sweep(pf, G, t0k, dt)
             end do
             call call_hooks(pf, G%level, PF_POST_SWEEP)
             call pf_residual(pf, G, dt)
          end do
       endif    ! (Pipeline_G)

       ! return to fine level...
       call pf_v_cycle_post_predictor(pf, t0, dt)

    end if

    call end_timer(pf, TPREDICTOR)
    call call_hooks(pf, 1, PF_POST_PREDICTOR)

    pf%state%iter   = 0
    pf%state%status = PF_STATUS_ITERATING

  end subroutine pf_predictor

  subroutine pf_check_convergence(pf, k, dt, res, qexit, qcycle)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(inout) :: res
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: k
    logical,           intent(out)   :: qexit, qcycle

    integer    :: steps_to_last
    real(pfdp) :: res1

    pf%state%nmoved = 0

    qexit  = .false.
    qcycle = .false.

    if (pf%window == PF_WINDOW_BLOCK .and. pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    res1 = pf%levels(pf%nlevels)%residual
    if (pf%state%status == PF_STATUS_ITERATING .and. res > 0.0d0) then
       if ( (abs(1.0_pfdp - abs(res1/res)) < pf%rel_res_tol) .or. &
            (abs(res1)                     < pf%abs_res_tol) ) then
          pf%state%status = PF_STATUS_CONVERGED
       end if
    end if
    res = res1

    pf%state%res = res

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)
    call pf_recv_status(pf, 8000+k)

    if (pf%rank /= pf%state%first .and. pf%state%pstatus == PF_STATUS_ITERATING) &
         pf%state%status = PF_STATUS_ITERATING

    call pf_send_status(pf, 8000+k)
    call call_hooks(pf, 1, PF_POST_CONVERGENCE)

    if (pf%window == PF_WINDOW_BLOCK) then

       ! XXX: this ain't so pretty, perhaps we should use the
       ! 'nmoved' thinger to break this cycle if everyone is
       ! done...

       if (pf%state%status == PF_STATUS_CONVERGED) then
          qcycle = .true.
          return
       end if

    else

       if (pf%state%status == PF_STATUS_ITERATING .and. pf%state%iter > pf%niters) then
          stop "failed to converge before max iteration count"
       end if

       if (pf%state%status == PF_STATUS_CONVERGED) then

          if (pf%rank == pf%state%last .and. pf%rank == pf%state%first) then
             qexit = .true.
             return
          end if

          if (pf%rank == pf%state%last) then
             pf%state%nmoved = pf%comm%nproc
             call pf_send_nmoved(pf, PF_TAG_NMOVED)
          else
             call pf_recv_nmoved(pf, PF_TAG_NMOVED)
          end if

          pf%state%pstatus = PF_STATUS_ITERATING
          pf%state%status  = PF_STATUS_ITERATING
          pf%state%step    = pf%state%step + pf%comm%nproc
          pf%state%t0      = pf%state%step * dt
          pf%state%iter    = 1
          res = -1

       else if (pf%state%pstatus == PF_STATUS_CONVERGED) then

          call pf_send_nmoved(pf, PF_TAG_NMOVED)

       end if

    end if

    pf%state%first  = modulo(pf%state%first + pf%state%nmoved, pf%comm%nproc)
    pf%state%last   = modulo(pf%state%last  + pf%state%nmoved, pf%comm%nproc)

    if (pf%state%step >= pf%state%nsteps) then
       qexit = .true.
       return
    end if

    ! roll back "last" processor
    steps_to_last = modulo(pf%state%last - pf%rank, pf%comm%nproc)
    do while (pf%state%step + steps_to_last >= pf%state%nsteps)
       pf%state%last = modulo(pf%state%last - 1, pf%comm%nproc)
       steps_to_last = modulo(pf%state%last - pf%rank, pf%comm%nproc)
    end do

    if (pf%state%nmoved == pf%comm%nproc) then
       pf%state%status = PF_STATUS_PREDICTOR
       qcycle = .true.
       return
    end if

  end subroutine pf_check_convergence

  !
  ! Run in parallel using PFASST.
  !
  subroutine pf_pfasst_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout) :: pf
    type(c_ptr),       intent(in)    :: q0
    real(pfdp),        intent(in)    :: dt, tend
    type(c_ptr),       intent(in), optional :: qend
    integer,           intent(in), optional :: nsteps

    type(pf_level_t), pointer :: F, G
    integer                   :: j, k, l
    real(pfdp)                :: res1

    logical :: qexit, qcycle, qbroadcast

    call start_timer(pf, TTOTAL)


    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%block   = 1
    pf%state%t0      = pf%state%step * dt
    pf%state%iter    = -1
    pf%state%cycle   = -1
    pf%state%first   = 0
    pf%state%last    = pf%comm%nproc - 1
    pf%state%status  = PF_STATUS_PREDICTOR
    pf%state%pstatus = PF_STATUS_PREDICTOR
    pf%comm%statreq  = -66

    res1 = -1

    F => pf%levels(pf%nlevels)
    call F%encap%pack(F%q0, q0)

    if (present(nsteps)) then
       pf%state%nsteps = nsteps
    else
       pf%state%nsteps = ceiling(1.0*tend/dt)
    end if
    do k = 1, 666666666

       qbroadcast = .false.

       ! in block mode, jump to next block if we've reached the max iteration count
       if (pf%window == PF_WINDOW_BLOCK .and. pf%state%iter >= pf%niters) then
          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt

          if (pf%state%step >= pf%state%nsteps) exit

          pf%state%status = PF_STATUS_PREDICTOR
          pf%state%block  = pf%state%block + 1
          qbroadcast = .true.
       end if

       ! in ring mode, if all procs moved at once, broadcast
       if (pf%window == PF_WINDOW_RING .and. pf%state%status == PF_STATUS_PREDICTOR) then
          qbroadcast = .true.
       end if

       if (k > 1 .and. qbroadcast) then
          F => pf%levels(pf%nlevels)
          call pf%comm%wait(pf, pf%nlevels)
          call F%encap%pack(F%send, F%qend)
          call pf_broadcast(pf, F%send, F%nvars, pf%comm%nproc-1)
          F%q0 = F%send
       end if
       ! predictor, if requested
       if (pf%state%status == PF_STATUS_PREDICTOR) &
            call pf_predictor(pf, pf%state%t0, dt)

       !
       ! perform fine sweeps
       !

       pf%state%iter   = pf%state%iter + 1
       pf%state%cycle  = 1

       call start_timer(pf, TITERATION)
       call call_hooks(pf, -1, PF_PRE_ITERATION)

       ! this if statement is necessary for block mode cycling...
       if (pf%state%status /= PF_STATUS_CONVERGED) then

          F => pf%levels(pf%nlevels)
          call call_hooks(pf, F%level, PF_PRE_SWEEP)
          do j = 1, F%nsweeps
             call F%sweeper%sweep(pf, F, pf%state%t0, dt)
          end do
          call call_hooks(pf, F%level, PF_POST_SWEEP)
          call pf_residual(pf, F, dt)

       end if
       !
       ! check convergence, continue with iteration
       !

       call pf_check_convergence(pf, k, dt, res1, qexit, qcycle)

       if (qexit)  exit
       if (qcycle) cycle
       do l = 2, pf%nlevels
          F => pf%levels(l)
          call pf_post(pf, F, F%level*10000+k)
       end do

       if (pf%state%status /= PF_STATUS_CONVERGED) then

          F => pf%levels(pf%nlevels)
          call pf_send(pf, F, F%level*10000+k, .false.)

          if (pf%nlevels > 1) then
             G => pf%levels(pf%nlevels-1)
             call restrict_time_space_fas(pf, pf%state%t0, dt, F, G)
             call save(G)
          end if

       end if

       call pf_v_cycle(pf, k, pf%state%t0, dt)
       call call_hooks(pf, -1, PF_POST_ITERATION)
       call end_timer(pf, TITERATION)

    end do

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)

    if (present(qend)) then
       F => pf%levels(pf%nlevels)
       call F%encap%copy(qend, F%qend)
    end if
  end subroutine pf_pfasst_run

  !
  ! After predictor, return to fine level.
  !
  subroutine pf_v_cycle_post_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt

    type(pf_level_t), pointer :: F, G
    integer :: l, j

    if (pf%nlevels <= 1) return

    do l = 2, pf%nlevels-1
       F => pf%levels(l); G => pf%levels(l-1)
       call interpolate_time_space(pf, t0, dt, F, G, G%Finterp)
       call G%encap%pack(F%q0, F%Q(1))
       call call_hooks(pf, F%level, PF_PRE_SWEEP)
       do j = 1, F%nsweeps
          call F%sweeper%sweep(pf, F, t0, dt)
       end do
       call call_hooks(pf, F%level, PF_POST_SWEEP)
       call pf_residual(pf, F, dt)
    end do

    F => pf%levels(pf%nlevels); G => pf%levels(pf%nlevels-1)
    call interpolate_time_space(pf, t0, dt, F, G, G%Finterp)
    call G%encap%pack(F%q0, F%Q(1))

  end subroutine pf_v_cycle_post_predictor

  !
  ! Execute a V-cycle, starting and ending from the middle level.
  !
  subroutine pf_v_cycle(pf, iteration, t0, dt)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration

    type(pf_level_t), pointer :: F, G
    integer :: l, j

    if (pf%nlevels == 1) then
       F => pf%levels(1)
       call pf_recv(pf, F, F%level*10000+iteration, .true.)
       return
    end if

    !
    ! down
    !
    do l = pf%nlevels-1, 2, -1
       F => pf%levels(l); G => pf%levels(l-1)
       call call_hooks(pf, F%level, PF_PRE_SWEEP)
       do j = 1, F%nsweeps
          call F%sweeper%sweep(pf, F, t0, dt)
       end do
       call call_hooks(pf, F%level, PF_POST_SWEEP)
       call pf_residual(pf, F, dt)
       call pf_send(pf, F, F%level*10000+iteration, .false.)
       call restrict_time_space_fas(pf, t0, dt, F, G)
       call save(G)
    end do

    !
    ! bottom
    !
    F => pf%levels(1)
    call pf_recv(pf, F, F%level*10000+iteration, .true.)
    call call_hooks(pf, F%level, PF_PRE_SWEEP)
    do j = 1, F%nsweeps
       call F%sweeper%sweep(pf, F, t0, dt)
    end do
    call call_hooks(pf, F%level, PF_POST_SWEEP)
    call pf_residual(pf, F, dt)
    call pf_send(pf, F, F%level*10000+iteration, .true.)

    !
    ! up
    !
    do l = 2, pf%nlevels
       F => pf%levels(l); G => pf%levels(l-1)

       call interpolate_time_space(pf, t0, dt, F, G,G%Finterp)
       call pf_recv(pf, F, F%level*10000+iteration, .false.)

       if (pf%rank /= pf%state%first) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          call interpolate_q0(pf,F, G)
       end if

       if (F%level < pf%nlevels) then
          call call_hooks(pf, F%level, PF_PRE_SWEEP)
          do j = 1, F%nsweeps
             call F%sweeper%sweep(pf, F, t0, dt)
          end do
          call call_hooks(pf, F%level, PF_POST_SWEEP)
          call pf_residual(pf, F, dt)
       end if
    end do

  end subroutine pf_v_cycle

  !
  ! Communication helpers
  !
  subroutine pf_post(pf, level, tag)
    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    if (pf%rank /= pf%state%first .and. pf%state%pstatus == PF_STATUS_ITERATING) then
       call pf%comm%post(pf, level, tag)
    end if
  end subroutine pf_post

  subroutine pf_send_status(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    if (pf%rank /= pf%state%last) then
       call pf%comm%send_status(pf, tag)
    end if
  end subroutine pf_send_status

  subroutine pf_recv_status(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    if (pf%rank /= pf%state%first) then
       call pf%comm%recv_status(pf, tag)
    end if
  end subroutine pf_recv_status

  subroutine pf_send_nmoved(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    call pf%comm%send_nmoved(pf, tag)
  end subroutine pf_send_nmoved

  subroutine pf_recv_nmoved(pf, tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: tag
    call pf%comm%recv_nmoved(pf, tag)
  end subroutine pf_recv_nmoved

  subroutine pf_send(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    call start_timer(pf, TSEND + level%level - 1)
    if (pf%rank /= pf%state%last &
         .and. pf%state%status == PF_STATUS_ITERATING) then
       call pf%comm%send(pf, level, tag, blocking)
    end if
    call end_timer(pf, TSEND + level%level - 1)
  end subroutine pf_send

  subroutine pf_recv(pf, level, tag, blocking)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking
    call start_timer(pf, TRECEIVE + level%level - 1)
    if (pf%rank /= pf%state%first .and.  pf%state%pstatus == PF_STATUS_ITERATING) then
       call pf%comm%recv(pf, level, tag, blocking)
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
