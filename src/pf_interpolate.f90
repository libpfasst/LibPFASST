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

  ! Interpolate (in time and space) level_index-1 to level_index
  !
  ! Interpolation is done by interpolating increments.  The fine
  ! function values are re-evaluated after interpolation.
  subroutine interpolate_time_space(pf, t0, dt, level_index, F_INTERP)
    type(pf_pfasst_t), intent(inout),target :: pf
    real(pfdp),        intent(in)    :: t0
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: level_index    !< defines which level to interpolate
    logical,           intent(in)    :: F_INTERP !<  Flag, if true, then do interp on f not sol
    !  Local variables
    class(pf_level_t), pointer :: c_lev_ptr
    class(pf_level_t), pointer :: f_lev_ptr

    integer    :: m, p
    real(pfdp), allocatable :: c_times(:)
    real(pfdp), allocatable :: f_times(:)

    class(pf_encap_t), allocatable :: delG(:)    !  Coarse in time and space
    class(pf_encap_t), allocatable :: delGF(:)   !  Coarse in time but fine in space

    f_lev_ptr => pf%levels(level_index)
    c_lev_ptr => pf%levels(level_index-1)

    call call_hooks(pf, level_index, PF_PRE_INTERP_ALL)
    call start_timer(pf, TINTERPOLATE + level_index - 1)

    ! create workspaces
    call c_lev_ptr%ulevel%factory%create_array(delG,  c_lev_ptr%nnodes, &
      c_lev_ptr%level, SDC_KIND_CORRECTION, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%create_array(delGF, c_lev_ptr%nnodes, &
      f_lev_ptr%level, SDC_KIND_CORRECTION, f_lev_ptr%nvars, f_lev_ptr%shape)

    ! set time at coarse and fine nodes
    allocate(c_times(c_lev_ptr%nnodes))
    allocate(f_times(f_lev_ptr%nnodes))

    c_times = t0 + dt*c_lev_ptr%nodes
    f_times = t0 + dt*f_lev_ptr%nodes

    ! needed for amr
    do m = 1, c_lev_ptr%nnodes
       call delG(m)%setval(0.0_pfdp)
       call delGF(m)%setval(0.0_pfdp)
    end do

    !!  Interpolating q
    do m = 1, c_lev_ptr%nnodes
       call delG(m)%copy(c_lev_ptr%Q(m))
       call delG(m)%axpy(-1.0_pfdp, c_lev_ptr%pQ(m))
       call f_lev_ptr%ulevel%interpolate(f_lev_ptr,c_lev_ptr, delGF(m), delG(m), c_times(m))
    end do

    ! interpolate corrections
    call pf_apply_mat(f_lev_ptr%Q, 1.0_pfdp, f_lev_ptr%tmat, delGF, .false.)

    if (F_INTERP) then
       !!  Interpolating F
       do p = 1,size(c_lev_ptr%F(1,:))
          do m = 1, c_lev_ptr%nnodes
             call delG(m)%setval(0.0_pfdp)
             call delGF(m)%setval(0.0_pfdp)
          end do
          do m = 1, c_lev_ptr%nnodes
            call delG(m)%copy(c_lev_ptr%F(m,p))
            call delG(m)%axpy(-1.0_pfdp, c_lev_ptr%pF(m,p))

            call f_lev_ptr%ulevel%interpolate(f_lev_ptr, c_lev_ptr, delGF(m), delG(m), c_times(m))
         end do

         ! interpolate corrections  in time
          call pf_apply_mat(f_lev_ptr%F(:,p), 1.0_pfdp, f_lev_ptr%tmat, delGF, .false.)

       end do !  Loop on npieces

       !  Do interpolation of qSDC(1) to update initial condition
       !  Not needed if we are interpolating solutions too
       !  call c_lev_ptr%encap%setval(delG(1),   0.0_pfdp)
       !  call c_lev_ptr%encap%setval(delGF(1),  0.0_pfdp)
       !  call c_lev_ptr%encap%copy(delG(1), c_lev_ptr%Q(1))
       !  call c_lev_ptr%encap%axpy(delG(1), -1.0_pfdp, c_lev_ptr%pQ(1))
       ! call f_lev_ptr%ulevel%interpolate(delGF(1), delG(1), f_lev_ptr%level, f_lev_ptr%levelctx, c_lev_ptr%level, c_lev_ptr%levelctx,c_times(1))

       ! This updates all solutions with jump in initial data
       ! interpolate corrections
       ! do n = 1, f_lev_ptr%nnodes
       !   call f_lev_ptr%encap%axpy(f_lev_ptr%Q(n), f_lev_ptr%tmat(n,1), delGF(1))
       ! end do

       ! recompute solution from new F, this is the old solution plus Residual
       !       call pf_integrate(pf,f_lev_ptr, dt)
       !       do m = 1, f_lev_ptr%nnodes-1
       !          call f_lev_ptr%encap%axpy(f_lev_ptr%Q(m+1), 1.0_pfdp, f_lev_ptr%R(m))
       !      end do
    else
       ! recompute fs
       call f_lev_ptr%ulevel%sweeper%evaluate_all(f_lev_ptr, f_times)
    end if  !  Feval

    !  Reset qend so that it is up to date
    call f_lev_ptr%qend%copy(f_lev_ptr%Q(f_lev_ptr%nnodes))

    call c_lev_ptr%ulevel%factory%destroy_array(delG,  c_lev_ptr%nnodes, &
      c_lev_ptr%level, SDC_KIND_CORRECTION, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_array(delGF, c_lev_ptr%nnodes, &
      f_lev_ptr%level, SDC_KIND_CORRECTION, f_lev_ptr%nvars, f_lev_ptr%shape)

    call end_timer(pf, TINTERPOLATE + f_lev_ptr%level - 1)
    call call_hooks(pf, f_lev_ptr%level, PF_POST_INTERP_ALL)
  end subroutine interpolate_time_space

  subroutine interpolate_q0(pf, f_lev_ptr, c_lev_ptr)
    !  Use to update the fine initial condition from increment

    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: f_lev_ptr, c_lev_ptr

    class(pf_encap_t), allocatable ::    delG, delF
    class(pf_encap_t), allocatable ::    q0F,q0G

    call call_hooks(pf, f_lev_ptr%level, PF_PRE_INTERP_Q0)
    call start_timer(pf, TINTERPOLATE + f_lev_ptr%level - 1)

    ! create workspaces
    call c_lev_ptr%ulevel%factory%create_single(q0G,  f_lev_ptr%level, SDC_KIND_SOL_NO_FEVAL, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%create_single(q0F,  f_lev_ptr%level, SDC_KIND_SOL_NO_FEVAL, f_lev_ptr%nvars, f_lev_ptr%shape)
    call c_lev_ptr%ulevel%factory%create_single(delG, c_lev_ptr%level, SDC_KIND_CORRECTION,   c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%create_single(delF, f_lev_ptr%level, SDC_KIND_CORRECTION,   f_lev_ptr%nvars, f_lev_ptr%shape)

    ! needed for amr
    call q0F%setval(0.0_pfdp)
    call q0G%setval(0.0_pfdp)
    call delG%setval(0.0_pfdp)
    call delF%setval(0.0_pfdp)

    call q0G%unpack(c_lev_ptr%q0)
    call q0F%unpack(f_lev_ptr%q0)

    call f_lev_ptr%ulevel%restrict(f_lev_ptr, c_lev_ptr, q0F, delG, pf%state%t0)
    call delG%axpy(-1.0_pfdp, q0G)

    call f_lev_ptr%ulevel%interpolate(f_lev_ptr, c_lev_ptr, delF, delG, pf%state%t0)
    call q0F%axpy(-1.0_pfdp, delF)
    call q0F%pack(f_lev_ptr%q0)

    call end_timer(pf, TINTERPOLATE + f_lev_ptr%level - 1)
    call call_hooks(pf, f_lev_ptr%level, PF_POST_INTERP_Q0)

    ! destroy workspaces
    call c_lev_ptr%ulevel%factory%destroy_single(q0G,  f_lev_ptr%level, SDC_KIND_SOL_NO_FEVAL, c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_single(q0F,  f_lev_ptr%level, SDC_KIND_SOL_NO_FEVAL, f_lev_ptr%nvars, f_lev_ptr%shape)
    call c_lev_ptr%ulevel%factory%destroy_single(delG, c_lev_ptr%level, SDC_KIND_CORRECTION,   c_lev_ptr%nvars, c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_single(delF, f_lev_ptr%level, SDC_KIND_CORRECTION,   f_lev_ptr%nvars, f_lev_ptr%shape)

  end subroutine interpolate_q0

end module pf_mod_interpolate
