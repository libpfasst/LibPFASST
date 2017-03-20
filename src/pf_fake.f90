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

module pf_mod_fake
  use pf_mod_dtype
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_hooks
  use pf_mod_parallel
  implicit none
contains

  !
  ! Predictor.
  !
  subroutine pf_fake_predictor(pf0, t0, dt)
    type(pf_pfasst_t), intent(inout) :: pf0
    real(pfdp),        intent(in)    :: t0, dt

    type(pf_level_t), pointer :: F, G
    integer :: j, k, l, p, nproc

    type(pf_pfasst_t), pointer :: pf
    type(c_ptr),       pointer :: pfs(:)

    nproc =  pf0%comm%nproc
    pfs   => pf0%comm%pfs

    do p = 1, nproc
       call c_f_pointer(pfs(p), pf)

       call call_hooks(pf, 1, PF_PRE_PREDICTOR)
       call start_timer(pf, TPREDICTOR)

       F => pf%levels(pf%nlevels)
       call spreadq0(F, t0)

       do l = pf%nlevels, 2, -1
          F => pf%levels(l); G => pf%levels(l-1)
          call pf_residual(pf, F, dt)
          call restrict_time_space_fas(pf, t0, dt, F, G)
          call save(G)
          call G%Q(1)%pack(G%q0)
       end do
    end do

    if (nproc > 1) then

       do k = 1, nproc
          do p = k, nproc
             call c_f_pointer(pfs(p), pf)

             pf%state%iter = -k

             G => pf%levels(1)

             ! get new initial value (skip on first iteration)
             if (k > 1) then
                call G%Q(G%nnodes)%pack(G%q0)
                call spreadq0(G, t0)
             end if

             call call_hooks(pf, G%level, PF_PRE_SWEEP)
             do j = 1, G%nsweeps
                call G%ulevel%sweeper%sweep(pf, G, t0, dt)
             end do
             call call_hooks(pf, G%level, PF_POST_SWEEP)
             call pf_residual(pf, G, dt)
          end do

       end do

    end if

    do p = 1, nproc
       call c_f_pointer(pfs(p), pf)

       pf%state%iter  = -1
       pf%state%cycle = -1
       call end_timer(pf, TPREDICTOR)
       call call_hooks(pf, 1, PF_POST_PREDICTOR)
    end do

  end subroutine pf_fake_predictor


  !
  ! Run PFASST using fake communicator.
  !
  subroutine pf_fake_run(pf0, q0, dt, nsteps)
    ! use pf_mod_parallel, only: pf_do_stage

    type(pf_pfasst_t), intent(inout) :: pf0
    type(c_ptr),       intent(in)    :: q0
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: nsteps

    ! type(pf_level_t), pointer :: F, G, Fb
    ! real(pfdp) :: t0
    ! integer    :: nblock, b, j, k, p, nproc

    ! type(pf_pfasst_t), pointer :: pf
    ! type(c_ptr),       pointer :: pfs(:)

    ! nproc =  pf0%comm%nproc
    ! pfs   => pf0%comm%pfs

    ! !
    ! ! initialize
    ! !

    ! do p = 1, nproc
    !    call c_f_pointer(pfs(p), pf)

    !    call start_timer(pf, TTOTAL)

    !    pf%state%t0   = 0.0d0
    !    pf%state%dt   = dt
    !    pf%state%iter = -1

    !    pf%state%first = 0
    !    pf%state%last  = pf%comm%nproc - 1

    !    F => pf%levels(pf%nlevels)
    !    call F%encap%pack(F%q0, q0)

    !    pf%state%nsteps = nsteps
    ! end do


    ! nblock = pf0%state%nsteps/pf0%comm%nproc

    ! !
    ! ! time "block" loop
    ! !

    ! do b = 1, nblock

    !    do p = 1, nproc
    !       call c_f_pointer(pfs(p), pf)

    !       pf%state%block  = b
    !       pf%state%step   = pf%rank + (b-1)*pf%comm%nproc
    !       pf%state%t0     = pf%state%step * dt
    !       pf%state%iter   = -1
    !       pf%state%cycle  = -1

    !       pf%state%status  = PF_STATUS_ITERATING
    !       pf%state%pstatus = PF_STATUS_ITERATING

    !       t0 = pf%state%t0

    !       call call_hooks(pf, -1, PF_PRE_BLOCK)
    !    end do

    !    call pf_fake_predictor(pf0, t0, dt)

    !    ! do start cycle stages
    !    do p = 1, nproc
    !       call c_f_pointer(pfs(p), pf)
    !       call pf_v_cycle_post_predictor(pf, t0, dt)
    !    end do

    !    ! pfasst iterations
    !    do k = 1, pf0%niters

    !       do p = 1, nproc
    !          call c_f_pointer(pfs(p), pf)

    !          pf%state%iter  = k
    !          pf%state%cycle = 1

    !          call start_timer(pf, TITERATION)
    !          call call_hooks(pf, -1, PF_PRE_ITERATION)

    !          ! XXX: skipping send/receive status

    !          F => pf%levels(pf%nlevels)
    !          call call_hooks(pf, F%level, PF_PRE_SWEEP)
    !          do j = 1, F%nsweeps
    !             call F%sweeper%sweep(pf, F, pf%state%t0, dt)
    !          end do
    !          call call_hooks(pf, F%level, PF_POST_SWEEP)
    !          call pf_residual(pf, F, dt)
    !          call pf_send(pf, F, F%level*10000+k, .false.)
    !          if (pf%nlevels > 1) then
    !             G => pf%levels(pf%nlevels-1)
    !             call restrict_time_space_fas(pf, pf%state%t0, dt, F, G)
    !             call save(G)
    !          end if

    !          call pf_v_cycle(pf, k, pf%state%t0, dt)

    !       end do

    !       do p = 1, nproc
    !          call call_hooks(pf, pf%nlevels, PF_POST_ITERATION)
    !          call end_timer(pf, TITERATION)
    !       end do

    !    end do ! end pfasst iteration loop

    !    ! do end cycle stages
    !    do p = 1, nproc
    !       call c_f_pointer(pfs(p), pf)
    !       ! ! end_stage%type = SDC_CYCLE_SWEEP
    !       ! ! end_stage%F    = pf%nlevels
    !       ! ! end_stage%G    = -1
    !       ! ! call pf_do_stage(pf, end_stage, -1, t0, dt)

    !       !    F => pf%levels(pf%nlevels)
    !       !    call call_hooks(pf, F%level, PF_PRE_SWEEP)
    !       !    do j = 1, F%nsweeps
    !       !       call F%sweeper%sweep(pf, F, pf%state%t0, dt)
    !       !    end do
    !       !    call call_hooks(pf, F%level, PF_POST_SWEEP)
    !       !    call pf_residual(pf, F, dt)
    !    end do

    !    do p = 1, nproc
    !       call c_f_pointer(pfs(p), pf)
    !       call call_hooks(pf, -1, PF_POST_STEP)
    !    end do

    !    ! broadcast fine qend (non-pipelined time loop)
    !    if (nblock > 1) then
    !       call c_f_pointer(pfs(nproc), pf)
    !       F => pf%levels(pf%nlevels)

    !       do p = 1, nproc
    !          call c_f_pointer(pfs(p), pf)
    !          Fb => pf%levels(pf%nlevels)
    !          call Fb%encap%copy(Fb%qend, F%qend)
    !          call Fb%encap%pack(Fb%q0, Fb%qend)
    !       end do
    !    end if

    ! end do ! end block loop

    ! !
    ! ! finish up
    ! !
    ! do p = 1, nproc
    !    call c_f_pointer(pfs(p), pf)
    !    pf%state%iter = -1
    !    call end_timer(pf, TTOTAL)
    ! end do

  end subroutine pf_fake_run

end module pf_mod_fake
