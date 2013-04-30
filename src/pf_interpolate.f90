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

module pf_mod_interpolate
  use pf_mod_dtype
  use pf_mod_restrict
  use pf_mod_timer
  use pf_mod_hooks
  use pf_mod_utils
  implicit none
contains

  ! Interpolate (in time and space) G to F.
  !
  ! Interpolation is done by interpolating increments.  The fine
  ! function values are re-evaluated after interpolation.
  subroutine interpolate_time_space(pf, t0, dt, F, G, Finterp)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: F, G
    logical,           intent(in), optional :: Finterp !  if true, then do interp on f not q

    integer    :: m, n, p
    real(pfdp) :: tm(F%nnodes)

    type(c_ptr) :: &
         delG(G%nnodes), &
         delGF(G%nnodes)

    call call_hooks(pf, F%level, PF_PRE_INTERP_ALL)
    call start_timer(pf, TINTERPOLATE + F%level - 1)

    ! create workspaces
    do m = 1, G%nnodes
       call G%encap%create(delG(m),   G%level, SDC_KIND_CORRECTION, &
            G%nvars, G%shape, G%ctx, G%encap%ctx)
       call F%encap%create(delGF(m),  F%level, SDC_KIND_CORRECTION, &
            F%nvars, F%shape, F%ctx, F%encap%ctx)
    end do

    if(present(Finterp) .and. (Finterp)) then
       !!  Interpolating F
       do p = 1,size(G%F(1,:))
          ! needed for amr
          do m = 1, G%nnodes
             call G%encap%setval(delG(m),   0.0_pfdp)
             call F%encap%setval(delGF(m),  0.0_pfdp)
          end do

          do m = 1, G%nnodes
             call G%encap%copy(delG(m), G%F(m,p))
             call G%encap%axpy(delG(m), -1.0_pfdp, G%pF(m,p))

             call G%interpolate(delGF(m), delG(m), F%level, F%ctx, G%level, G%ctx)
          end do

          ! interpolate corrections
          do n = 1, F%nnodes
             do m = 1, G%nnodes
                call G%encap%axpy(F%F(n,p), F%tmat(n,m), delGF(m))
             end do
          end do
       end do !  Loop on npieces

       !  Do interpolation of qSDC(1) to update initial condition
       call G%encap%setval(delG(1),   0.0_pfdp)
       call G%encap%setval(delGF(1),  0.0_pfdp)
       call G%encap%copy(delG(1), G%Q(1))
       call G%encap%axpy(delG(1), -1.0_pfdp, G%pQ(1))

       call F%interpolate(delGF(1), delG(1), F%level, F%ctx, G%level, G%ctx)

       ! interpolate corrections
       do n = 1, F%nnodes
          call F%encap%axpy(F%Q(n), F%tmat(n,1), delGF(1))
       end do

       ! ! recompute fs
       ! tm = t0 + dt*F%nodes
       ! call sdceval(tm(1), 1, F)
    else
       !!  Interpolating q
       ! create workspaces
       do m = 1, G%nnodes
          ! needed for amr
          call G%encap%setval(delG(m),   0.0_pfdp)
          call G%encap%setval(delGF(m),  0.0_pfdp)
       end do

       do m = 1, G%nnodes
          call G%encap%copy(delG(m), G%Q(m))
          call G%encap%axpy(delG(m), -1.0_pfdp, G%pQ(m))

          call F%interpolate(delGF(m), delG(m), F%level, F%ctx, G%level, G%ctx)
       end do

       ! interpolate corrections
       call pf_apply_mat(F%Q, 1.0_pfdp, F%tmat, delGF, F%encap, .false.)

       ! recompute fs
       tm = t0 + dt*F%nodes
       do m = 1, F%nnodes
          call F%sweeper%evaluate(F, tm(m), m)
       end do
    end if  !  Feval

    ! destroy workspaces
    do m = 1, G%nnodes
       call G%encap%destroy(delG(m))
       call G%encap%destroy(delGF(m))
    end do

    call end_timer(pf, TINTERPOLATE + F%level - 1)
    call call_hooks(pf, F%level, PF_POST_INTERP_ALL)
  end subroutine interpolate_time_space

  subroutine interpolate_q0(pf, F, G)
    !  Use to update the fine initial condition from increment

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: F, G

    type(c_ptr) ::    delG, delF
    type(c_ptr) ::    q0F,q0G

    call call_hooks(pf, F%level, PF_PRE_INTERP_Q0)
    call start_timer(pf, TINTERPOLATE + F%level - 1)

    ! create workspaces
    call G%encap%create(q0G,  F%level, SDC_KIND_SOL_NO_FEVAL, &
         G%nvars, G%shape, G%ctx, G%encap%ctx)
    call F%encap%create(q0F,  F%level, SDC_KIND_SOL_NO_FEVAL, &
         F%nvars, F%shape, F%ctx, F%encap%ctx)
    call G%encap%create(delG, G%level, SDC_KIND_CORRECTION, &
         G%nvars, G%shape, G%ctx, G%encap%ctx)
    call F%encap%create(delF, F%level, SDC_KIND_CORRECTION, &
         F%nvars, F%shape, F%ctx, F%encap%ctx)

    ! needed for amr
    call F%encap%setval(q0F,  0.0_pfdp)
    call G%encap%setval(q0G,  0.0_pfdp)
    call G%encap%setval(delG, 0.0_pfdp)
    call F%encap%setval(delF, 0.0_pfdp)

    call G%encap%unpack(q0G, G%q0)
    call F%encap%unpack(q0F, F%q0)

    call F%restrict(q0F, delG, F%level, F%ctx, G%level, G%ctx)
    call G%encap%axpy(delG, -1.0_pfdp, q0G)

    call F%interpolate(delF, delG, F%level, F%ctx, G%level, G%ctx)
    call F%encap%axpy(q0F, -1.0_pfdp, delF)

    call F%encap%pack(F%q0, q0F)
    call G%encap%destroy(delG)
    call F%encap%destroy(delF)
    call F%encap%destroy(q0F)
    call G%encap%destroy(q0G)

    call end_timer(pf, TINTERPOLATE + F%level - 1)
    call call_hooks(pf, F%level, PF_POST_INTERP_Q0)
  end subroutine interpolate_q0

end module pf_mod_interpolate
