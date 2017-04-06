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
  subroutine interpolate_time_space(pf, t0, dt, LevF, LevG, F_INTERP)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    class(pf_level_t),  intent(inout) :: LevF, LevG
    logical,           intent(in), optional :: F_INTERP !<  Flag, if true, then do interp on f not sol

    integer    :: m, p
    real(pfdp) :: tF(LevF%nnodes)
    real(pfdp) :: tG(LevG%nnodes)
    logical :: F_INTERP_LOC

    class(pf_encap_t), allocatable :: delG(:)    !  Coarse in time and space
    class(pf_encap_t), allocatable :: delGF(:)   !  Coarse in time but fine in space

    call call_hooks(pf, LevF%level, PF_PRE_INTERP_ALL)
    call start_timer(pf, TINTERPOLATE + LevF%level - 1)

    ! create workspaces
    call LevG%ulevel%factory%create_array(delG,  LevG%nnodes, LevG%level, SDC_KIND_CORRECTION, LevG%nvars, LevG%shape)
    call LevF%ulevel%factory%create_array(delGF, LevG%nnodes, LevF%level, SDC_KIND_CORRECTION, LevF%nvars, LevF%shape)

    ! set time at coarse and fine nodes
    tG = t0 + dt*LevG%nodes
    tF = t0 + dt*LevF%nodes

    ! needed for amr
    do m = 1, LevG%nnodes
       call delG(m)%setval(0.0_pfdp)
       call delGF(m)%setval(0.0_pfdp)
    end do

    !!  Interpolating q
    do m = 1, LevG%nnodes
       call delG(m)%copy(LevG%Q(m))
       call delG(m)%axpy(-1.0_pfdp, LevG%pQ(m))
       call LevF%ulevel%interpolate(LevF, LevG, delGF(m), delG(m), tG(m))
    end do

    ! interpolate corrections
    call pf_apply_mat(LevF%Q, 1.0_pfdp, LevF%tmat, delGF, .false.)

    F_INTERP_LOC = .FALSE.
    if(present(F_INTERP)) then
       if  (F_INTERP)  then
          F_INTERP_LOC = .TRUE.
       end if
    end if

    if (F_INTERP_LOC) then
       !!  Interpolating F
       do p = 1,size(LevG%F(1,:))
          do m = 1, LevG%nnodes
             call delG(m)%setval(0.0_pfdp)
             call delGF(m)%setval(0.0_pfdp)
          end do
          do m = 1, LevG%nnodes
            call delG(m)%copy(LevG%F(m,p))
            call delG(m)%axpy(-1.0_pfdp, LevG%pF(m,p))

            call LevF%ulevel%interpolate(LevF, LevG, delGF(m), delG(m), tG(m))
         end do

         ! interpolate corrections  in time
          call pf_apply_mat(LevF%F(:,p), 1.0_pfdp, LevF%tmat, delGF, .false.)

       end do !  Loop on npieces

       !  Do interpolation of qSDC(1) to update initial condition
       !  Not needed if we are interpolating solutions too
       !  call LevG%encap%setval(delG(1),   0.0_pfdp)
       !  call LevG%encap%setval(delGF(1),  0.0_pfdp)
       !  call LevG%encap%copy(delG(1), LevG%Q(1))
       !  call LevG%encap%axpy(delG(1), -1.0_pfdp, LevG%pQ(1))
       ! call LevF%ulevel%interpolate(delGF(1), delG(1), LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,tG(1))

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
       !   call LevF%ulevel%sweeper%evaluate(LevF, tF(m), m)
       !  end do
    else
       ! recompute fs
!       do m = 1, LevF%nnodes
!          call LevF%ulevel%sweeper%evaluate(LevF, tF(m), m)
!       end do
       call LevF%ulevel%sweeper%evaluate_all(LevF, tF)
    end if  !  Feval

    !  Reset qend so that it is up to date
    call LevF%qend%copy(LevF%Q(LevF%nnodes))

    call LevG%ulevel%factory%destroy_array(delG,  LevG%nnodes, LevG%level, SDC_KIND_CORRECTION, LevG%nvars, LevG%shape)
    call LevF%ulevel%factory%destroy_array(delGF, LevG%nnodes, LevF%level, SDC_KIND_CORRECTION, LevF%nvars, LevF%shape)

    call end_timer(pf, TINTERPOLATE + LevF%level - 1)
    call call_hooks(pf, LevF%level, PF_POST_INTERP_ALL)
  end subroutine interpolate_time_space

  subroutine interpolate_q0(pf, LevF, LevG)
    !  Use to update the fine initial condition from increment

    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: LevF, LevG

    class(pf_encap_t), allocatable ::    delG, delF
    class(pf_encap_t), allocatable ::    q0F,q0G

    call call_hooks(pf, LevF%level, PF_PRE_INTERP_Q0)
    call start_timer(pf, TINTERPOLATE + LevF%level - 1)

    ! create workspaces
    call LevG%ulevel%factory%create_single(q0G,  LevF%level, SDC_KIND_SOL_NO_FEVAL, LevG%nvars, LevG%shape)
    call LevF%ulevel%factory%create_single(q0F,  LevF%level, SDC_KIND_SOL_NO_FEVAL, LevF%nvars, LevF%shape)
    call LevG%ulevel%factory%create_single(delG, LevG%level, SDC_KIND_CORRECTION,   LevG%nvars, LevG%shape)
    call LevF%ulevel%factory%create_single(delF, LevF%level, SDC_KIND_CORRECTION,   LevF%nvars, LevF%shape)

    ! needed for amr
    call q0F%setval(0.0_pfdp)
    call q0G%setval(0.0_pfdp)
    call delG%setval(0.0_pfdp)
    call delF%setval(0.0_pfdp)

    call q0G%unpack(LevG%q0)
    call q0F%unpack(LevF%q0)

    call LevF%ulevel%restrict(LevF, LevG, q0F, delG, pf%state%t0)
    call delG%axpy(-1.0_pfdp, q0G)

    call LevF%ulevel%interpolate(LevF, levG, delF, delG, pf%state%t0)
    call q0F%axpy(-1.0_pfdp, delF)
    call q0F%pack(LevF%q0)

    call end_timer(pf, TINTERPOLATE + LevF%level - 1)
    call call_hooks(pf, LevF%level, PF_POST_INTERP_Q0)

    ! destroy workspaces
    call LevG%ulevel%factory%destroy_single(q0G,  LevF%level, SDC_KIND_SOL_NO_FEVAL, LevG%nvars, LevG%shape) 
    call LevF%ulevel%factory%destroy_single(q0F,  LevF%level, SDC_KIND_SOL_NO_FEVAL, LevF%nvars, LevF%shape) 
    call LevG%ulevel%factory%destroy_single(delG, LevG%level, SDC_KIND_CORRECTION,   LevG%nvars, LevG%shape) 
    call LevF%ulevel%factory%destroy_single(delF, LevF%level, SDC_KIND_CORRECTION,   LevF%nvars, LevF%shape) 

  end subroutine interpolate_q0

end module pf_mod_interpolate
