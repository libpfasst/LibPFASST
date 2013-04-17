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
!   2013-04-17 - Matthew Emmett
!
!     Time restriction was switched from point injection polynomial
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
  implicit none
contains

  ! Restrict (in time and space) qF to qG
  subroutine restrict_sdc(F, G, qF, qG)
    use pf_mod_utils

    type(pf_level_t),  intent(inout) :: F, G
    type(c_ptr),       intent(inout) :: qF(F%nnodes), qG(G%nnodes)

    type(c_ptr) :: qFr(F%nnodes)
    integer :: m, n

    ! create workspaces
    do m = 1, F%nnodes
       call G%encap%create(qFr(m), G%level, SDC_KIND_SOL_NO_FEVAL, &
            G%nvars, G%shape, G%ctx, G%encap%ctx)
    end do

    ! restrict and apply restriction matrix
    do m = 1, F%nnodes
       call F%restrict(qF(m), qFr(m), F%level, F%ctx, G%level, G%ctx)
    end do

    call pf_apply_mat(qG, 1.0_pfdp, F%rmat, qFr, F%encap)

    ! tidy
    do m = 1, F%nnodes
       call G%encap%destroy(qFr(m))
    end do

  end subroutine restrict_sdc


  ! Restrict (in time and space) F to G and set G's FAS correction.
  !
  ! The coarse function values are re-evaluated after restriction.
  subroutine restrict_time_space_fas(pf, t0, dt, F, G)
    use pf_mod_dtype
    use pf_mod_utils
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: F, G

    integer    :: m, mc, trat
    real(pfdp) :: tm(G%nnodes)
    type(c_ptr) :: &
         tmpG(G%nnodes), &    ! coarse integral of coarse function values
         tmpF(F%nnodes), &    ! fine integral of fine function values
         GofF(G%nnodes)       ! coarse integral of restricted fine function values

    call start_timer(pf, TRESTRICT + F%level - 1)

    trat = (F%nnodes - 1) / (G%nnodes - 1)

    ! note: even if the number of variables and nodes is the same, we
    ! should still compute the fas correction since the function
    ! evaluations may be different

    !!!! create workspaces
    do m = 1, G%nnodes
       call G%encap%create(tmpG(m), G%level, SDC_KIND_INTEGRAL, G%nvars, G%shape, G%ctx, G%encap%ctx)
       call G%encap%create(GofF(m), G%level, SDC_KIND_INTEGRAL, G%nvars, G%shape, G%ctx, G%encap%ctx)
    end do

    do m = 1, F%nnodes
       call F%encap%create(tmpF(m), F%level, SDC_KIND_INTEGRAL, F%nvars, F%shape, F%ctx, F%encap%ctx)
    end do

    !!!! restrict q's and recompute f's
    tm = t0 + dt*G%nodes

    call restrict_sdc(F, G, F%Q, G%Q)

    do m = 1, G%nnodes
       call G%sweeper%evaluate(G, tm(m), m)
    end do

    !!!! bring down fas correction from level above
    do m = 1, G%nnodes-1
       call G%encap%setval(G%tau(m), 0.0_pfdp)
    end do

    if (associated(F%tau)) then
       ! convert from 'node to node' to '0 to node'
       call F%encap%setval(tmpF(1), 0.0_pfdp)
       do m = 2, F%nnodes
          call F%encap%copy(tmpF(m), tmpF(m-1))
          call F%encap%axpy(tmpF(m), 1.0_pfdp, F%tau(m-1))
       end do

       call restrict_sdc(F, G, tmpF, tmpG)

       ! convert from '0 to node' to 'node to node'
       call G%encap%copy(G%tau(1), tmpG(2))
       do m = 2, G%nnodes-1
          call G%encap%copy(G%tau(m), tmpG(m+1))
          call G%encap%axpy(G%tau(m), -1.0_pfdp, tmpG(m))
       end do
    end if

    !!!! fas correction

    ! compute '0 to node' integrals (Q mat), and set the value at the
    ! first node to 0.0 (which is the integral from 0 to 0)

    call G%encap%setval(tmpG(1), 0.0_pfdp)
    call F%encap%setval(tmpF(1), 0.0_pfdp)

    call pf_integrate(G, dt, G%Qmat, tmpG(2:))
    call pf_integrate(F, dt, F%Qmat, tmpF(2:))

    ! now restrict the fine integral (using rmat) and compute the
    ! 'node to node' tau correction

    call restrict_sdc(F, G, tmpF, GofF)

    do m = 1, G%nnodes-1
       call G%encap%axpy(G%tau(m),  1.0_pfdp, GofF(m+1))
       call G%encap%axpy(G%tau(m), -1.0_pfdp, GofF(m))

       call G%encap%axpy(G%tau(m), -1.0_pfdp, tmpG(m+1))
       call G%encap%axpy(G%tau(m),  1.0_pfdp, tmpG(m))
    end do

    !!!! destroy workspaces
    do m = 1, G%nnodes
       call G%encap%destroy(tmpG(m))
       call G%encap%destroy(GofF(m))
    end do

    do m = 1, F%nnodes
       call F%encap%destroy(tmpF(m))
    end do

    call end_timer(pf, TRESTRICT + F%level - 1)

  end subroutine restrict_time_space_fas

end module pf_mod_restrict

