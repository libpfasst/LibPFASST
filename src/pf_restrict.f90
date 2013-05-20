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
  subroutine restrict_sdc(F, G, qF, qG, integral)
    use pf_mod_utils, only: pf_apply_mat

    type(pf_level_t), intent(inout) :: F, G
    type(c_ptr),      intent(inout) :: qF(:), qG(:)
    logical,          intent(in)    :: integral

    type(c_ptr), pointer :: qFr(:)
    integer :: m

    if (integral) then

       allocate(qFr(F%nnodes-1))

       do m = 1, F%nnodes-1
          call G%encap%create(qFr(m), G%level, SDC_KIND_INTEGRAL, &
               G%nvars, G%shape, G%ctx, G%encap%ctx)
          call F%restrict(qF(m), qFr(m), F%level, F%ctx, G%level, G%ctx)
       end do

       ! when restricting '0 to node' integral terms, skip the first
       ! entry since it is zero
       call pf_apply_mat(qG, 1.d0, F%rmat(2:,2:), qFr, G%encap)

    else

       allocate(qFr(F%nnodes))

       do m = 1, F%nnodes
          call G%encap%create(qFr(m), G%level, SDC_KIND_SOL_NO_FEVAL, &
               G%nvars, G%shape, G%ctx, G%encap%ctx)
          call F%restrict(qF(m), qFr(m), F%level, F%ctx, G%level, G%ctx)
       end do

       call pf_apply_mat(qG, 1.d0, F%rmat, qFr, G%encap)

    end if

    do m = 1, size(qFr)
       call G%encap%destroy(qFr(m))
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
  subroutine restrict_time_space_fas(pf, t0, dt, F, G)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: F, G

    integer    :: m
    real(pfdp) :: tm(G%nnodes)
    type(c_ptr) :: &
         tmpG(G%nnodes), &    ! coarse integral of coarse function values
         tmpF(F%nnodes), &    ! fine integral of fine function values
         tmpFr(G%nnodes)      ! coarse integral of restricted fine function values

    call call_hooks(pf, F%level, PF_PRE_RESTRICT_ALL)
    call start_timer(pf, TRESTRICT + F%level - 1)

    !
    ! create workspaces
    !
    do m = 1, G%nnodes
       call G%encap%create(tmpG(m), G%level, SDC_KIND_INTEGRAL, &
            G%nvars, G%shape, G%ctx, G%encap%ctx)
       call G%encap%create(tmpFr(m), G%level, SDC_KIND_INTEGRAL, &
            G%nvars, G%shape, G%ctx, G%encap%ctx)
    end do

    do m = 1, F%nnodes
       call F%encap%create(tmpF(m), F%level, SDC_KIND_INTEGRAL, &
            F%nvars, F%shape, F%ctx, F%encap%ctx)
    end do


    !
    ! restrict q's and recompute f's
    !
    tm = t0 + dt*G%nodes

    call restrict_sdc(F, G, F%Q, G%Q, .false.)

    do m = 1, G%nnodes
       call G%sweeper%evaluate(G, tm(m), m)
    end do


    !
    ! fas correction
    !

    do m = 1, G%nnodes-1
       call G%encap%setval(G%tau(m), 0.0_pfdp)
    end do

    ! compute '0 to node' integral on the coarse level

    call G%sweeper%integrate(G, G%Q, G%F, dt, tmpG)
    do m = 2, G%nnodes-1
       call G%encap%axpy(tmpG(m), 1.0_pfdp, tmpG(m-1))
    end do

    ! restrict '0 to node' integral on the fine level (which was
    ! computed during the last call to pf_residual)

    call restrict_sdc(F, G, F%I, tmpFr, .true.)

    ! compute 'node to node' tau correction

    call G%encap%axpy(G%tau(1),  1.0_pfdp, tmpFr(1))
    call G%encap%axpy(G%tau(1), -1.0_pfdp, tmpG(1))

    do m = 2, G%nnodes-1
       call G%encap%axpy(G%tau(m),  1.0_pfdp, tmpFr(m))
       call G%encap%axpy(G%tau(m), -1.0_pfdp, tmpFr(m-1))

       call G%encap%axpy(G%tau(m), -1.0_pfdp, tmpG(m))
       call G%encap%axpy(G%tau(m),  1.0_pfdp, tmpG(m-1))
    end do


    !
    ! tidy
    !

    do m = 1, G%nnodes
       call G%encap%destroy(tmpG(m))
       call G%encap%destroy(tmpFr(m))
    end do

    do m = 1, F%nnodes
       call F%encap%destroy(tmpF(m))
    end do

    call end_timer(pf, TRESTRICT + F%level - 1)
    call call_hooks(pf, F%level, PF_POST_RESTRICT_ALL)

  end subroutine restrict_time_space_fas

end module pf_mod_restrict

