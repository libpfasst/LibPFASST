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
  implicit none
contains

  ! Interpolate (in time and space) G to F.
  !
  ! Interpolation is done by interpolating increments.  The fine
  ! function values are re-evaluated after interpolation.
  subroutine interpolate_time_space(pf, t0, dt, F, G, Finterp)
    use pf_mod_dtype
    use pf_mod_restrict
    use pf_mod_sweep
    use pf_mod_timer
    use transfer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: F, G
    logical,           intent(in), optional :: Finterp !  if true, then do interp on f not q

    integer    :: m, n, p, trat
    real(pfdp) :: tm(F%nnodes)

    type(pf_encap_t) :: &
         delG(G%nnodes), &
         delGF(G%nnodes)

    call start_timer(pf, TINTERPOLATE + F%level - 1)

    ! create workspaces
    do m = 1, G%nnodes
       call create(delG(m),   G%level, .true., G%nvars, G%shape, G%ctx)
       call create(delGF(m),  F%level, .true., F%nvars, F%shape, F%ctx)
    end do

    if(present(Finterp) .and. (Finterp)) then
       !!  Interpolating F
       do p = 1,size(G%fSDC(1,:))
          ! needed for amr
          do m = 1, G%nnodes
             call setval(delG(m),   0.0_pfdp)
             call setval(delGF(m),  0.0_pfdp)
          end do

          do m = 1, G%nnodes
             call copy(delG(m), G%fSDC(m,p))
             call axpy(delG(m), -1.0_pfdp, G%pfSDC(m,p))

             call interpolate(delGF(m), delG(m), F%level, F%ctx, G%level, G%ctx)
          end do

          ! interpolate corrections
          trat = (F%nnodes-1) / (G%nnodes-1)
          do n = 1, F%nnodes
             do m = 1, G%nnodes
                call axpy(F%fSDC(n,p), F%tmat(n,m), delGF(m))
             end do
          end do
       end do !  Loop on npieces

       !  Do interpolation of qSDC(1) to update initial condition
       call setval(delG(1),   0.0_pfdp)
       call setval(delGF(1),  0.0_pfdp)
       call copy(delG(1), G%qSDC(1))
       call axpy(delG(1), -1.0_pfdp, G%pSDC(1))

       call interpolate(delGF(1), delG(1), F%level, F%ctx, G%level, G%ctx)

       ! interpolate corrections
       trat = (F%nnodes-1) / (G%nnodes-1)

       do n = 1, F%nnodes
          call axpy(F%qSDC(n), F%tmat(n,1), delGF(1))
       end do

       ! ! recompute fs
       ! tm = t0 + dt*F%nodes
       ! call sdceval(tm(1), 1, F)
    else
       !!  Interpolating q
       ! create workspaces
       do m = 1, G%nnodes
          ! needed for amr
          call setval(delG(m),   0.0_pfdp)
          call setval(delGF(m),  0.0_pfdp)
       end do

       do m = 1, G%nnodes
          call copy(delG(m), G%qSDC(m))
          call axpy(delG(m), -1.0_pfdp, G%pSDC(m))

          call interpolate(delGF(m), delG(m), F%level, F%ctx, G%level, G%ctx)
       end do

       ! interpolate corrections
       trat = (F%nnodes-1) / (G%nnodes-1)

       do n = 1, F%nnodes
          do m = 1, G%nnodes
             call axpy(F%qSDC(n), F%tmat(n,m), delGF(m))
          end do
       end do

       ! recompute fs
       tm = t0 + dt*F%nodes
       do m = 1, F%nnodes
          call sdceval(tm(m), m, F)
       end do
    end if  !  Feval

    ! destroy workspaces
    do m = 1, G%nnodes
       call destroy(delG(m))
       call destroy(delGF(m))
    end do

    call end_timer(pf, TINTERPOLATE + F%level - 1)
  end subroutine interpolate_time_space

  subroutine interpolate_q0(pf, F, G)
    !  Use to update the fine initial condition from increment
    use pf_mod_dtype
    use pf_mod_restrict
    use pf_mod_sweep
    use pf_mod_timer
    use transfer

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: F, G

    type(pf_encap_t) ::    delG, delF
    type(pf_encap_t) ::    q0F,q0G

    call start_timer(pf, TINTERPOLATE + F%level - 1)

    ! create workspaces
    call create(q0G,  F%level, .true., G%nvars, G%shape, G%ctx)
    call create(q0F,  F%level, .true., F%nvars, F%shape, F%ctx)
    call create(delG,   G%level, .true., G%nvars, G%shape, G%ctx)
    call create(delF,  F%level, .true., F%nvars, F%shape, F%ctx)

    ! needed for amr
    call setval(q0F,  0.0_pfdp)
    call setval(q0G,  0.0_pfdp)
    call setval(delG,   0.0_pfdp)
    call setval(delF,  0.0_pfdp)

    call unpack(q0G,G%q0)
    call unpack(q0F,F%q0)

    call restrict(q0F, delG, F%level, F%ctx, G%level, G%ctx)

    call axpy(delG, -1.0_pfdp, q0G)

    call interpolate(delF, delG, F%level, F%ctx, G%level, G%ctx)

    call axpy(q0F, -1.0_pfdp, delF)

    call pack(F%q0, q0F)
    call destroy(delG)
    call destroy(delF)
    call destroy(q0F)
    call destroy(q0G)

    call end_timer(pf, TINTERPOLATE + F%level - 1)
  end subroutine interpolate_q0

end module pf_mod_interpolate
