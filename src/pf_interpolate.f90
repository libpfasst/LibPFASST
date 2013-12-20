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

  ! Interpolate (in time and space) LevG to LevF.
  !
  ! Interpolation is done by interpolating increments.  The fine
  ! function values are re-evaluated after interpolation.
  subroutine interpolate_time_space(pf, t0, dt, LevF, LevG, Finterp)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: LevF, LevG
    logical,           intent(in), optional :: Finterp !  if true, then do interp on f not q

    integer    :: m, n, p
    real(pfdp) :: tF(LevF%nnodes)
    real(pfdp) :: tG(LevG%nnodes)
    logical :: Finterp_loc

    type(c_ptr) ::  delG(LevG%nnodes)   !  Coarse in time and space 
    type(c_ptr) ::  delGF(LevG%nnodes)  !  Coarse in time but fine in space

    call call_hooks(pf, LevF%level, PF_PRE_INTERP_ALL)
    call start_timer(pf, TINTERPOLATE + LevF%level - 1)

    ! create workspaces
    do m = 1, LevG%nnodes
       call LevG%encap%create(delG(m),   LevG%level, SDC_KIND_CORRECTION, &
            LevG%nvars, LevG%shape, LevG%levelctx, LevG%encap%encapctx)
       call LevF%encap%create(delGF(m),  LevF%level, SDC_KIND_CORRECTION, &
            LevF%nvars, LevF%shape, LevF%levelctx, LevF%encap%encapctx)
    end do

    ! set time at coarse and fine nodes
    tG = t0 + dt*LevG%nodes
    tF = t0 + dt*LevF%nodes


    ! needed for amr
    do m = 1, LevG%nnodes
       call LevG%encap%setval(delG(m),   0.0_pfdp)
       call LevF%encap%setval(delGF(m),  0.0_pfdp)
    end do

    !!  Interpolating q
    do m = 1, LevG%nnodes
       call LevG%encap%copy(delG(m), LevG%Q(m))
       call LevG%encap%axpy(delG(m), -1.0_pfdp, LevG%pQ(m))
       call LevF%interpolate(delGF(m), delG(m), LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,tG(m))
    end do
    
    ! interpolate corrections
    call pf_apply_mat(LevF%Q, 1.0_pfdp, LevF%tmat, delGF, LevF%encap, .false.)

    Finterp_loc = .FALSE.
    if(present(Finterp)) then
       if  (Finterp)  then
          Finterp_loc = .TRUE.
       end if
    end if

    if (Finterp_loc) then
       !!  Interpolating F
       do p = 1,size(LevG%F(1,:))
          do m = 1, LevG%nnodes
             call LevG%encap%setval(delG(m),   0.0_pfdp)
             call LevF%encap%setval(delGF(m),  0.0_pfdp)
          end do
          do m = 1, LevG%nnodes
            call LevG%encap%copy(delG(m), LevG%F(m,p))
            call LevG%encap%axpy(delG(m), -1.0_pfdp, LevG%pF(m,p))

            call LevF%interpolate(delGF(m), delG(m), LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,tG(m))
         end do

!           do n = 1, LevF%nnodes
!               do m = 1, LevG%nnodes
!                  call LevF%encap%axpy(LevF%F(n,p), LevF%tmat(n,m), delGF(m))
!               end do
!            end do
!!$          
         ! interpolate corrections 
          call pf_apply_mat(LevF%F(:,p), 1.0_pfdp, LevF%tmat, delGF, LevF%encap, .false.)
!!$
       end do !  Loop on npieces

       !  Do interpolation of qSDC(1) to update initial condition
       !  Not needed if we are interpolating solutions too
       !  call LevG%encap%setval(delG(1),   0.0_pfdp)
       !  call LevG%encap%setval(delGF(1),  0.0_pfdp)
       !  call LevG%encap%copy(delG(1), LevG%Q(1))
       !  call LevG%encap%axpy(delG(1), -1.0_pfdp, LevG%pQ(1))
       ! call LevF%interpolate(delGF(1), delG(1), LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,tG(1))

       ! This updates all solutions with jump in initial data       
       ! interpolate corrections
       ! do n = 1, LevF%nnodes
       !   call LevF%encap%axpy(LevF%Q(n), LevF%tmat(n,1), delGF(1))
       ! end do

       ! recompute solution from new F, this is the old solution plus Residual
!       call pf_integrate(pf,LevF, dt)
!       do m = 1, LevF%nnodes-1
!          call LevF%encap%axpy(LevF%Q(m+1), 1.0_pfdp, LevF%R(m))
!      end do
      
       ! recompute fs (for debugging)
       !       do m = 1, LevF%nnodes
       !   call LevF%sweeper%evaluate(LevF, tF(m), m)
       !  end do
    else
       ! recompute fs
       do m = 1, LevF%nnodes
          call LevF%sweeper%evaluate(LevF, tF(m), m)
       end do
    end if  !  Feval

    !  Reset qend so that it is up to date
    call LevF%encap%copy(LevF%qend, LevF%Q(LevF%nnodes))
    ! destroy workspaces
    do m = 1, LevG%nnodes
       call LevG%encap%destroy(delG(m))
       call LevG%encap%destroy(delGF(m))
    end do

    call end_timer(pf, TINTERPOLATE + LevF%level - 1)
    call call_hooks(pf, LevF%level, PF_POST_INTERP_ALL)
  end subroutine interpolate_time_space

  subroutine interpolate_q0(pf, LevF, LevG)
    !  Use to update the fine initial condition from increment

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: LevF, LevG

    type(c_ptr) ::    delG, delF
    type(c_ptr) ::    q0F,q0G

    call call_hooks(pf, LevF%level, PF_PRE_INTERP_Q0)
    call start_timer(pf, TINTERPOLATE + LevF%level - 1)
  
    ! create workspaces
    call LevG%encap%create(q0G,  LevF%level, SDC_KIND_SOL_NO_FEVAL, &
         LevG%nvars, LevG%shape, LevG%levelctx, LevG%encap%encapctx)
    call LevF%encap%create(q0F,  LevF%level, SDC_KIND_SOL_NO_FEVAL, &
         LevF%nvars, LevF%shape, LevF%levelctx, LevF%encap%encapctx)
    call LevG%encap%create(delG, LevG%level, SDC_KIND_CORRECTION, &
         LevG%nvars, LevG%shape, LevG%levelctx, LevG%encap%encapctx)
    call LevF%encap%create(delF, LevF%level, SDC_KIND_CORRECTION, &
         LevF%nvars, LevF%shape, LevF%levelctx, LevF%encap%encapctx)

    ! needed for amr
    call LevF%encap%setval(q0F,  0.0_pfdp)
    call LevG%encap%setval(q0G,  0.0_pfdp)
    call LevG%encap%setval(delG, 0.0_pfdp)
    call LevF%encap%setval(delF, 0.0_pfdp)

    call LevG%encap%unpack(q0G, LevG%q0)
    call LevF%encap%unpack(q0F, LevF%q0)

    call LevF%restrict(q0F, delG, LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,pf%state%t0)
    call LevG%encap%axpy(delG, -1.0_pfdp, q0G)

    call LevF%interpolate(delF, delG, LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,pf%state%t0)
    call LevF%encap%axpy(q0F, -1.0_pfdp, delF)

    call LevF%encap%pack(LevF%q0, q0F)
    call LevG%encap%destroy(delG)
    call LevF%encap%destroy(delF)
    call LevF%encap%destroy(q0F)
    call LevG%encap%destroy(q0G)

    call end_timer(pf, TINTERPOLATE + LevF%level - 1)
    call call_hooks(pf, LevF%level, PF_POST_INTERP_Q0)
  end subroutine interpolate_q0

end module pf_mod_interpolate
