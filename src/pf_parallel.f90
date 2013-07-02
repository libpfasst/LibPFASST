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
  use pf_mod_cycle
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
  subroutine pf_predictor(pf, t0, dt)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt

    type(pf_level_t), pointer :: F, G
    integer :: j, k, l

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    F => pf%levels(pf%nlevels)
    call spreadq0(F, t0)

    do l = pf%nlevels, 2, -1
       F => pf%levels(l); G => pf%levels(l-1)
       call pf_residual(F, dt)
       call restrict_time_space_fas(pf, t0, dt, F, G)
       call save(G)
       call G%encap%pack(G%q0, G%Q(1))
    end do

    if (pf%comm%nproc > 1) then

       G => pf%levels(1)
       do k = 1, pf%rank + 1
          pf%state%iter = -k

          ! get new initial value (skip on first iteration)
          if (k > 1) &
               call G%encap%pack(G%q0, G%Q(G%nnodes))

          call call_hooks(pf, G%level, PF_PRE_SWEEP)
          do j = 1, G%nsweeps
             call G%sweeper%sweep(pf, G, t0, dt)
          end do
          call call_hooks(pf, G%level, PF_POST_SWEEP)
          call pf_residual(G, dt)
       end do

    end if

    pf%state%iter  = -1
    pf%state%cycle = -1
    call end_timer(pf, TPREDICTOR)
    call call_hooks(pf, 1, PF_POST_PREDICTOR)

  end subroutine pf_predictor


  !
  ! Execute a cycle "stage".
  !
  subroutine pf_do_stage(pf, stage, iteration, t0, dt)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_stage_t),  intent(in)    :: stage
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration

    type(pf_level_t), pointer :: F, G
    integer :: j

    select case(stage%type)
    case (SDC_CYCLE_UP)

       F => pf%levels(stage%F)
       G => pf%levels(stage%G)

       call interpolate_time_space(pf, t0, dt, F, G,G%Finterp)

       call pf%comm%recv(pf, F, F%level*100+iteration, .false.)

       if (pf%rank > 0) then
          ! interpolate increment to q0 -- the fine initial condition
          ! needs the same increment that Q(1) got, but applied to the
          ! new fine initial condition
          call interpolate_q0(pf,F, G)
       end if

       if (F%level < pf%nlevels) then
          ! don't sweep on the finest level on the way up -- if you
          ! need to do this, use a PF_CYCLE_SWEEP

          call call_hooks(pf, F%level, PF_PRE_SWEEP)
          do j = 1, F%nsweeps
             call F%sweeper%sweep(pf, F, t0, dt)
          end do
          call call_hooks(pf, F%level, PF_POST_SWEEP)
          call pf_residual(F, dt)
       end if


    case (SDC_CYCLE_DOWN)

       F => pf%levels(stage%F)
       G => pf%levels(stage%G)

       call call_hooks(pf, F%level, PF_PRE_SWEEP)
       do j = 1, F%nsweeps
          call F%sweeper%sweep(pf, F, t0, dt)
       end do
       call call_hooks(pf, F%level, PF_POST_SWEEP)
       call pf_residual(F, dt)
       call pf%comm%send(pf, F, F%level*100+iteration, .false.)

       call restrict_time_space_fas(pf, t0, dt, F, G)
       call save(G)


    case (SDC_CYCLE_BOTTOM)

       F => pf%levels(stage%F)

       call pf%comm%recv(pf, F, F%level*100+iteration, .true.)
       call call_hooks(pf, F%level, PF_PRE_SWEEP)
       do j = 1, F%nsweeps
          call F%sweeper%sweep(pf, F, t0, dt)
       end do
       call call_hooks(pf, F%level, PF_POST_SWEEP)
       call pf_residual(F, dt)
       call pf%comm%send(pf, F, F%level*100+iteration, .true.)


    case (SDC_CYCLE_SWEEP)

       F => pf%levels(stage%F)

       call call_hooks(pf, F%level, PF_PRE_SWEEP)
       do j = 1, F%nsweeps
          call F%sweeper%sweep(pf, F, t0, dt)
       end do
       call call_hooks(pf, F%level, PF_POST_SWEEP)
       call pf_residual(F, dt)

    case (SDC_CYCLE_INTERP)

       F => pf%levels(stage%F)
       G => pf%levels(stage%G)

       call interpolate_time_space(pf, t0, dt, F, G, G%Finterp)
       call G%encap%pack(F%q0, F%Q(1))

    end select

  end subroutine pf_do_stage


  subroutine pf_do_window_block(pf, dt, tend, nsteps)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, tend
    integer,           intent(in), optional :: nsteps

    type(pf_level_t), pointer :: F
    real(pfdp) :: t0, res1, res0
    integer    :: nblock, b, c, k, l


    nblock = pf%state%nsteps/pf%comm%nproc
    res0   = 1.d0

    !
    ! time "block" loop
    !

    pf%comm%statreq = -66

    do b = 1, nblock
       pf%state%block  = b
       pf%state%step   = pf%rank + (b-1)*pf%comm%nproc
       pf%state%t0     = pf%state%step * dt
       pf%state%iter   = -1
       pf%state%cycle  = -1

       t0 = pf%state%t0

       call call_hooks(pf, -1, PF_PRE_BLOCK)
       call start_timer(pf, TSTEP)

       call pf_predictor(pf, t0, dt)

       ! do start cycle stages
       if (associated(pf%cycles%start)) then
          do c = 1, size(pf%cycles%start)
             pf%state%cycle = c
             call pf_do_stage(pf, pf%cycles%start(c), -1, t0, dt)
          end do
       end if

       ! pfasst iterations
       do k = 1, pf%niters
          if (.not. (pf%state%status == PF_STATUS_CONVERGED)) then
             pf%state%iter  = k
          end if
          pf%state%cycle = 1

          call start_timer(pf, TITERATION)

          do l = 1, pf%nlevels
             F => pf%levels(l)
             call call_hooks(pf, F%level, PF_PRE_ITERATION)
          end do

          ! send/receive status
          if (pf%abs_res_tol > 0.d0 .or. pf%rel_res_tol > 0.d0) then
             F => pf%levels(pf%nlevels)

             res1 = F%residual

             if ((k > 1) .and. ( (abs(1.0_pfdp - res1/res0) < pf%rel_res_tol) &
                            .or. (res1 < pf%abs_res_tol) &
                            .or. (abs(res0-res1) < pf%abs_res_tol) )) then
                pf%state%status = PF_STATUS_CONVERGED
             else
                pf%state%status = PF_STATUS_ITERATING
             end if

             res0 = res1

             call pf%comm%recv_status(pf, 8000+k)

             ! if the previous processor hasn't converged yet, keep
             ! iterating
             if (pf%rank /= pf%state%first) then
                if (pf%state%pstatus /= PF_STATUS_CONVERGED) &
                   pf%state%status = PF_STATUS_ITERATING
             end if

             call pf%comm%send_status(pf, 8000+k)
          end if

          if (pf%state%status == PF_STATUS_CONVERGED) then
             cycle
          end if

          ! post receive requests
          do l = 2, pf%nlevels
             F => pf%levels(l)
             call pf%comm%post(pf, F, F%level*100+k)
          end do

          ! do pfasst cycle stages
          do c = 1, size(pf%cycles%pfasst)
             pf%state%cycle = pf%state%cycle + 1
             call pf_do_stage(pf, pf%cycles%pfasst(c), k, t0, dt)
          end do

          call call_hooks(pf, pf%nlevels, PF_POST_ITERATION)
          call end_timer(pf, TITERATION)

       end do ! end pfasst iteration loop

       ! do end cycle stages
       if (associated(pf%cycles%end)) then
          do c = 1, size(pf%cycles%end)
             pf%state%cycle = c
             call pf_do_stage(pf, pf%cycles%end(c), -1, t0, dt)
          end do
       end if

       call call_hooks(pf, -1, PF_POST_STEP)
       call end_timer(pf, TSTEP)

       ! broadcast fine qend (non-pipelined time loop)
       if (nblock > 1) then
          F => pf%levels(pf%nlevels)

          call pf%comm%wait(pf, pf%nlevels)
          call F%encap%pack(F%send, F%qend)
          call pf%comm%broadcast(pf, F%send, F%nvars, pf%comm%nproc-1)
          F%q0 = F%send
       end if

       pf%comm%statreq = -66

    end do ! end block loop
  end subroutine pf_do_window_block

  subroutine pf_do_window_ring(pf, dt, tend, nsteps)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, tend
    integer,           intent(in), optional :: nsteps

    type(pf_level_t), pointer :: F
    real(pfdp) :: t0, res1, res0
    integer    :: c, k, l
    integer    :: steps_to_last


    pf%comm%statreq = -66
    pf%state%block  = -66

    pf%state%step   = pf%rank
    pf%state%t0     = pf%state%step * dt
    pf%state%iter   = -1
    pf%state%cycle  = -1

    pf%state%status  = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING

    t0 = pf%state%t0

    call pf_predictor(pf, t0, dt)

    ! do start cycle stages
    if (associated(pf%cycles%start)) then
       do c = 1, size(pf%cycles%start)
          pf%state%cycle = c
          call pf_do_stage(pf, pf%cycles%start(c), -1, t0, dt)
       end do
    end if

    pf%state%nmoved = 0
    pf%state%status = PF_STATUS_PREDICTOR

    ! pfasst iterations
    do k = 1, 9999
       pf%state%iter   = k
       pf%state%cycle  = 1

       call start_timer(pf, TITERATION)

       do l = 1, pf%nlevels
          F => pf%levels(l)
          call call_hooks(pf, F%level, PF_PRE_ITERATION)
       end do

       ! check convergence
       res1 = pf%levels(pf%nlevels)%residual
       if (pf%state%status == PF_STATUS_ITERATING .and. res0 > 0.0d0) then
          if ( (abs(1.0_pfdp - abs(res1/res0)) < pf%rel_res_tol) .or. &
               (abs(res1)                      < pf%abs_res_tol) ) then

             pf%state%status = PF_STATUS_CONVERGED

          end if
       end if
       res0 = res1

       if (pf%state%status /= PF_STATUS_CONVERGED) &
         pf%state%status = PF_STATUS_ITERATING

       pf%state%nmoved = 0
       call pf%comm%recv_status(pf, 8000+k)

       ! keep iterating if the previous processor hasn't converged yet
       if (pf%rank /= pf%state%first) then
          if (pf%state%pstatus == PF_STATUS_ITERATING) &
               pf%state%status = PF_STATUS_ITERATING
       end if

       call pf%comm%send_status(pf, 8000+k)

       if (pf%state%status == PF_STATUS_CONVERGED) then

          if (pf%rank == pf%state%last .and. pf%rank == pf%state%first) then
             exit
          end if

          if (pf%rank == pf%state%last) then
             call pf%comm%send_nmoved(pf, PF_TAG_NMOVED)

             pf%state%status = PF_STATUS_ITERATING
          else
             call call_hooks(pf, -1, PF_POST_STEP)
             call pf%comm%recv_nmoved(pf, PF_TAG_NMOVED)

             pf%state%status = PF_STATUS_ITERATING
             pf%state%step   = pf%state%step + pf%comm%nproc
             res0 = -1
          end if

       else if (pf%state%pstatus == PF_STATUS_CONVERGED) then

          call pf%comm%send_nmoved(pf, PF_TAG_NMOVED)

       end if

       pf%state%t0     = pf%state%step * dt
       pf%state%first  = modulo(pf%state%first + pf%state%nmoved, pf%comm%nproc)
       pf%state%last   = modulo(pf%state%last  + pf%state%nmoved, pf%comm%nproc)

       if (pf%state%step >= pf%state%nsteps) exit

       ! roll back "last" processor
       steps_to_last = modulo(pf%state%last - pf%rank, pf%comm%nproc)
       do while (pf%state%step + steps_to_last >= pf%state%nsteps)
          pf%state%last = modulo(pf%state%last - 1, pf%comm%nproc)
          steps_to_last = modulo(pf%state%last - pf%rank, pf%comm%nproc)
       end do

       ! post receive requests
       do l = 2, pf%nlevels
          F => pf%levels(l)
          call pf%comm%post(pf, F, F%level*100+k)
       end do

       ! do pfasst cycle stages
       do c = 1, size(pf%cycles%pfasst)
          pf%state%cycle = pf%state%cycle + 1
          call pf_do_stage(pf, pf%cycles%pfasst(c), k, t0, dt)
       end do

       call call_hooks(pf, pf%nlevels, PF_POST_ITERATION)
       call end_timer(pf, TITERATION)

    end do ! end pfasst iteration loop

    ! do end cycle stages
    if (associated(pf%cycles%end)) then
       do c = 1, size(pf%cycles%end)
          pf%state%cycle = c
          call pf_do_stage(pf, pf%cycles%end(c), -1, t0, dt)
       end do
    end if

  end subroutine pf_do_window_ring


  !
  ! Run in parallel using PFASST.
  !
  subroutine pf_pfasst_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout) :: pf
    type(c_ptr),       intent(in)    :: q0
    real(pfdp),        intent(in)    :: dt, tend
    type(c_ptr),       intent(in), optional :: qend
    integer,           intent(in), optional :: nsteps

    type(pf_level_t), pointer :: F
    
    call start_timer(pf, TTOTAL)

    !
    ! initialize
    !

    pf%state%t0   = 0.0d0
    pf%state%dt   = dt
    pf%state%iter = -1

    pf%state%first = 0
    pf%state%last  = pf%comm%nproc - 1

    F => pf%levels(pf%nlevels)
    call F%encap%pack(F%q0, q0)

    if(present(nsteps)) then
       pf%state%nsteps = nsteps
    else
       pf%state%nsteps = ceiling(1.0*tend/dt)
    end if

    !
    ! dispatch
    !

    select case (pf%window)
    case (PF_WINDOW_BLOCK)
       call pf_do_window_block(pf, dt, tend, nsteps)
    case (PF_WINDOW_RING)
       call pf_do_window_ring(pf, dt, tend, nsteps)
    case default
       stop "ERROR: Invalid window type."
    end select

    !
    ! finish up
    !

    pf%state%iter = -1
    call end_timer(pf, TTOTAL)

    if (present(qend)) then
       F => pf%levels(pf%nlevels)
       call F%encap%copy(qend, F%qend)
    end if
  end subroutine pf_pfasst_run

end module pf_mod_parallel
