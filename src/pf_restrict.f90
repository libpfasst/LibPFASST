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

! Restriction and FAS routines.

!
! Notes:
!
!   2013-04-30 - Matthew Emmett
!
!     The pf_residual subroutine is now called after each SDC sweep,
!     and it computes the '0 to node' integrals and stores them in
!     'F%I' while it is computing the full SDC residual.  Furthermore,
!     these 'F%I' integrals also contain the appropriate tau corrections.
!
!     This means that when computing FAS corrections: the fine
!     integral term is already done for us, and it is already FAS
!     corrected, so we dont't have to "bring down fas corrections"
!     from finer levels.
!
!
!   2013-04-17 - Matthew Emmett
!
!     Time restriction was switched from point injection to polynomial
!     interpolation (ie, using the 'rmat's in each level) so that we
!     can use proper nodes for each level.
!
!     To recover point injection (ie, use copy instead of axpy)
!     properly we should really do some masking trickery with the
!     restriction matrices (rmat).  XXX.
!
!     Finally, perhaps the workspaces should be preallocated (along
!     with interpolation workspaces...).  XXX
!

module pf_mod_restrict
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_hooks
  implicit none
contains


  !
  ! Restrict (in time and space) qF to qG.
  !
  subroutine restrict_sdc(LevF, LevG, qF, qG, integral,tF)
    use pf_mod_utils, only: pf_apply_mat

    type(pf_level_t), intent(inout) :: LevF, LevG
    type(c_ptr),      intent(inout) :: qF(:), qG(:)
    logical,          intent(in)    :: integral
    real(pfdp),       intent(in) :: tF(:)

    type(c_ptr), pointer :: qFr(:)
    integer :: m

    if (integral) then

       allocate(qFr(LevF%nnodes-1))

       do m = 1, LevF%nnodes-1
          call LevG%encap%create(qFr(m), LevG%level, SDC_KIND_INTEGRAL, &
               LevG%nvars, LevG%shape, LevG%levelctx, LevG%encap%encapctx)
          call LevF%restrict(qF(m), qFr(m), LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,tF(m))
       end do

       ! when restricting '0 to node' integral terms, skip the first
       ! entry since it is zero
       call pf_apply_mat(qG, 1.d0, LevF%rmat(2:,2:), qFr, LevG%encap)
    else

       allocate(qFr(LevF%nnodes))

       do m = 1, LevF%nnodes
          call LevG%encap%create(qFr(m), LevG%level, SDC_KIND_SOL_NO_FEVAL, &
               LevG%nvars, LevG%shape, LevG%levelctx, LevG%encap%encapctx)
          call LevF%restrict(qF(m), qFr(m), LevF%level, LevF%levelctx, LevG%level, LevG%levelctx,tF(m))
       end do

       call pf_apply_mat(qG, 1.d0, LevF%rmat, qFr, LevG%encap)

    end if

    do m = 1, size(qFr)
       call LevG%encap%destroy(qFr(m))
    end do
    deallocate(qFr)

  end subroutine restrict_sdc


  !
  ! Restrict (in time and space) F to G and set G's FAS correction.
  !
  ! The coarse function values are re-evaluated after restriction.
  ! Note that even if the number of variables and nodes is the same,
  ! we should still compute the FAS correction since the function
  ! evaluations may be different.
  !
  subroutine restrict_time_space_fas(pf, t0, dt, LevF, LevG)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: LevF, LevG

    integer    :: m
    logical    :: old_way
    real(pfdp) :: tG(LevG%nnodes)
    real(pfdp) :: tF(LevF%nnodes)
    type(c_ptr) :: &
         tmpG(LevG%nnodes), &    ! coarse integral of coarse function values
         tmpF(LevF%nnodes), &    ! fine integral of fine function values
         tmpFr(LevG%nnodes)      ! coarse integral of restricted fine function values

    call call_hooks(pf, LevF%level, PF_PRE_RESTRICT_ALL)
    call start_timer(pf, TRESTRICT + LevF%level - 1)

    !
    ! create workspaces
    !
    do m = 1, LevG%nnodes
       call LevG%encap%create(tmpG(m), LevG%level, SDC_KIND_INTEGRAL, &
            LevG%nvars, LevG%shape, LevG%levelctx, LevG%encap%encapctx)
       call LevG%encap%create(tmpFr(m), LevG%level, SDC_KIND_INTEGRAL, &
            LevG%nvars, LevG%shape, LevG%levelctx, LevG%encap%encapctx)
    end do

    do m = 1, LevF%nnodes
       call LevF%encap%create(tmpF(m), LevF%level, SDC_KIND_INTEGRAL, &
            LevF%nvars, LevF%shape, LevF%levelctx, LevF%encap%encapctx)
    end do


    !
    ! restrict q's and recompute f's
    !
    tG = t0 + dt*LevG%nodes
    tF = t0 + dt*LevF%nodes

    call restrict_sdc(LevF, LevG, LevF%Q, LevG%Q, .false.,tF)

!    do m = 1, LevG%nnodes
!       call LevG%sweeper%evaluate(LevG, tG(m), m)
!    end do
     call LevG%sweeper%evaluate_all(LevG, tG)


    !
    ! fas correction
    !
    do m = 1, LevG%nnodes-1
       call LevG%encap%setval(LevG%tau(m), 0.0_pfdp)
       call LevG%encap%setval(LevG%tauQ(m), 0.0_pfdp)
    end do
    if (pf%state%iter >= pf%taui0)  then
       ! compute '0 to node' integral on the coarse level
       call LevG%sweeper%integrate(LevG, LevG%Q, LevG%F, dt, tmpG)
       !MMQ       do m = 2, LevG%nnodes-1
       !   call LevG%encap%axpy(tmpG(m), 1.0_pfdp, tmpG(m-1))
       !end do
       
       ! compute '0 to node' integral on the fine level
       call LevF%sweeper%integrate(LevF, LevF%Q, LevF%F, dt, LevF%I)

       !  put tau in
       !MMQ do m = 2, LevF%nnodes-1
       !   call LevF%encap%axpy(LevF%I(m), 1.0_pfdp, LevF%I(m-1))
       !end do
       if (associated(LevF%tauQ)) then
          do m = 1, LevF%nnodes-1
             call LevF%encap%axpy(LevF%I(m), 1.0_pfdp, LevF%tauQ(m))
          end do
       end if
       
       ! restrict '0 to node' integral on the fine level (which was
       ! computed during the last call to pf_residual)
       call restrict_sdc(LevF, LevG, LevF%I, tmpFr, .true.,tF)
       
       ! compute 'node to node' tau correction
       call LevG%encap%axpy(LevG%tau(1),  1.0_pfdp, tmpFr(1))
       call LevG%encap%axpy(LevG%tau(1), -1.0_pfdp, tmpG(1))
       
       do m = 2, LevG%nnodes-1
          call LevG%encap%axpy(LevG%tau(m),  1.0_pfdp, tmpFr(m))
          call LevG%encap%axpy(LevG%tau(m), -1.0_pfdp, tmpFr(m-1))
          
          call LevG%encap%axpy(LevG%tau(m), -1.0_pfdp, tmpG(m))
          call LevG%encap%axpy(LevG%tau(m),  1.0_pfdp, tmpG(m-1))
       end do
      ! compute '0 to node' tau correction
       do m = 1, LevG%nnodes-1
          call LevG%encap%axpy(LevG%tauQ(m),  1.0_pfdp, tmpFr(m))
          call LevG%encap%axpy(LevG%tauQ(m), -1.0_pfdp, tmpG(m))
       end do
    end if
 
       !
    ! tidy
    !

    do m = 1, LevG%nnodes
       call LevG%encap%destroy(tmpG(m))
       call LevG%encap%destroy(tmpFr(m))
    end do

    do m = 1, LevF%nnodes
       call LevF%encap%destroy(tmpF(m))
    end do

    call end_timer(pf, TRESTRICT + LevF%level - 1)
    call call_hooks(pf, LevF%level, PF_POST_RESTRICT_ALL)

  end subroutine restrict_time_space_fas

end module pf_mod_restrict

