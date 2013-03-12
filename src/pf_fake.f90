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

! PFASST routines implmented in serial.  We fake the action of MPI by
! accessing all the "processors" through pointers stored in the PFASST
! communicator.

module pf_mod_fake
  use iso_c_binding
  implicit none
contains

  ! Run in serial using the fake communicator.
  subroutine pfasst_fake_run(pf0, q0, dt, tend, nsteps, qend)
    use encap
    use pf_mod_comm

    use pf_mod_interpolate
    use pf_mod_restrict
    use pf_mod_sweep
    use pf_mod_utils
    use pf_mod_timer
    use pf_mod_hooks
    use transfer, only: interpolate

    type(pf_pfasst_t), intent(inout) :: pf0
    type(pf_encap_t),  intent(in)    :: q0
    real(pfdp),        intent(in)    :: dt, tend
    type(pf_encap_t),  intent(inout), optional :: qend
    integer,           intent(in),    optional :: nsteps

    real(pfdp) :: t0
    integer    :: nproc, p
    integer    :: nblock, b, k, j, l

    type(pf_pfasst_t), pointer :: pf     ! used in loops
    type(pf_level_t), pointer  :: F, G, Fb

    type(c_ptr), pointer :: pfs(:)

    nproc =  pf0%comm%nproc
    pfs   => pf0%comm%pfs

    !!!! set initial conditions
    do p = 1, nproc
       call c_f_pointer(pfs(p), pf)

       pf%state%t0   = 0.0d0
       pf%state%dt   = dt
       pf%state%iter = -1

       call start_timer(pf, TTOTAL)

       F => pf%levels(pf%nlevels)
       call pack(F%q0, q0)
    end do

    if(present(nsteps)) then
       nblock = nsteps/nproc
    else
       nblock = int(ceiling(tend/dt/nproc))
    end if


    !!!! time "block" loop
    do b = 1, nblock


       !!!! spread and restrict
       do p = 1, nproc
          call c_f_pointer(pfs(p), pf)
          pf%state%step = pf%rank + (b-1)*pf%comm%nproc
          pf%state%t0   = pf%state%step * dt

          t0 = pf%state%t0

          call start_timer(pf, TPREDICTOR)

          F => pf%levels(pf%nlevels)
          call spreadq0(F, t0)

          do l = pf%nlevels, 2, -1
             F => pf%levels(l); G => pf%levels(l-1)
             call restrict_time_space_fas(pf, t0, dt, F, G)
             call save(G)
             call pack(G%q0, G%qSDC(1))
          end do
       end do ! proc loop


       !!!! predictor loop
       if (nproc > 1) then
          do k = 1, nproc
             do p = k, nproc
                call c_f_pointer(pfs(p), pf)
                pf%state%iter = -k

                G => pf%levels(1)

                ! get new initial value (skip on first iteration)
                if (k > 1) &
                     call recv(pf, G, k-1, .true.)

                do j = 1, G%nsweeps
                   call sweep(pf, pf%state%t0, dt, G)
                end do
             end do
             do p = k, nproc
                call c_f_pointer(pfs(p), pf)
                G => pf%levels(1)
                call send(pf, G, k, .true.)
             end do ! proc loop
          end do

          do p = 1, nproc
             call c_f_pointer(pfs(p), pf)

             ! interpolate
             do l = 1, pf%nlevels - 1
                F => pf%levels(l+1)
                G => pf%levels(l)
                call interpolate_time_space(pf, pf%state%t0, dt, F, G,G%Finterp)
                call pack(F%q0, F%qSDC(1))
             end do
          end do ! proc loop
       end if

       do p = 1, nproc
          call c_f_pointer(pfs(p), pf)
          call end_timer(pf, TPREDICTOR)
          call call_hooks(pf, -1, PF_POST_PREDICTOR)
       end do

       !!!! pfasst iterations (v-cycle)
       do k = 1, pf%niters
          do p = 1, nproc
             call c_f_pointer(pfs(p), pf)
             pf%state%iter  = k
             pf%state%cycle = 1

             call start_timer(pf, TITERATION)
          end do ! proc loop

          do p = 1, nproc
             call c_f_pointer(pfs(p), pf)

             !! go down the v-cycle
             do l = pf%nlevels, 2, -1
                pf%state%cycle = pf%state%cycle + 1
                F => pf%levels(l)
                G => pf%levels(l-1)

                do j = 1, F%nsweeps
                   call sweep(pf, pf%state%t0, dt, F)
                end do
                call send(pf, F, F%level*100+k, .false.)

                call call_hooks(pf, F%level, PF_POST_SWEEP)
                call restrict_time_space_fas(pf, pf%state%t0, dt, F, G)
                call save(G)
             end do ! on l:  down the vcycle

             !! bottom
             pf%state%cycle = pf%state%cycle + 1
             F => pf%levels(1)

             call recv(pf, F, F%level*100+k, .true.)
             do j = 1, F%nsweeps
                call sweep(pf, pf%state%t0, dt, F)
             end do
             call send(pf, F, F%level*100+k, .true.)

             call call_hooks(pf, F%level, PF_POST_SWEEP)

             !! go up the v-cycle
             do l = 2, pf%nlevels
                pf%state%cycle = pf%state%cycle + 1
                F => pf%levels(l)
                G => pf%levels(l-1)

                call interpolate_time_space(pf, pf%state%t0, dt, F, G,G%Finterp)

                call recv(pf, F, F%level*100+k, .false.)

                if (pf%rank > 0) then
                   ! XXX: correction...
                   ! XXX: using qSDC(1) for now but this is redundant (beginning of next sweep)
                   call unpack(G%qSDC(1), G%q0)
                   call interpolate(F%qSDC(1), G%qSDC(1), F%level, F%ctx, G%level, G%ctx)
                   call pack(F%q0, F%qSDC(1))
                end if

                if (l < pf%nlevels) then
                   do j = 1, F%nsweeps
                      call sweep(pf, pf%state%t0, dt, F)
                   end do
                   call call_hooks(pf, F%level, PF_POST_SWEEP)
                end if
             end do  ! on l:  up the v-cycle

             ! XXX
             ! if (F%residual_tol > 0.0_pfdp) then
             !    call pf_compute_residual(F%q0, dt, F%fSDC, &
             !         F%nvars, F%nnodes, F%smat, F%qend, res)
             !    if (res < F%residual_tol) then
             !       print '("serial residual condition met after: ",i3," iterations",i3)', k
             !       exit
             !    end if
             ! end if

             call call_hooks(pf, -1, PF_POST_ITERATION)
             call end_timer(pf, TITERATION)
          end do ! proc loop
       end do ! iter loop

       do p = 1, nproc
          call c_f_pointer(pfs(p), pf)
          F => pf%levels(pf%nlevels)
          call call_hooks(pf, F%level, PF_POST_STEP)
       end do     !  on p: -------   Proc loop ----------
       !   MMM not sure about theis broadcast and what it does
       ! broadcast fine qend (non-pipelined time loop)
       call pack(F%send, F%qend)
       F => pf%levels(pf%nlevels)

       ! Doing looping now.  F still has the final processor,
       ! So put it's qend into all the q0
       if (nblock > 1) then
          !  This is the broadcast of initial data
          do p = 1,nproc
             call c_f_pointer(pfs(p), pf)
             Fb => pf%levels(pf%nlevels)
             call copy(Fb%qend, F%qend)
             call pack(Fb%q0,Fb%qend)
          end do
       end if


    end do ! on b, block loops

    do p = 1,nproc
       call c_f_pointer(pfs(p), pf)
       pf%state%iter = -1
       call end_timer(pf, TTOTAL)
    end do

    if (present(qend)) then
       F => pf%levels(pf%nlevels)
       call copy(qend, F%qend)
    end if
  end subroutine pfasst_fake_run

end module pf_mod_fake
