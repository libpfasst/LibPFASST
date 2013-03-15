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

module pf_mod_restrict
  implicit none
contains

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
         CofG(G%nnodes-1), &    ! coarse integral of coarse function values
         FofF(F%nnodes-1), &    ! fine integral of fine function values
         CofF(G%nnodes-1), &    ! coarse integral of restricted fine function values
         tmp

    call start_timer(pf, TRESTRICT + F%level - 1)

    trat = (F%nnodes - 1) / (G%nnodes - 1)

    ! note: even if the number of variables and nodes is the same, we
    ! should still compute the fas correction since the function
    ! evaluations may be different

    !!!! create workspaces
    do m = 1, G%nnodes-1
       call G%encap%create(CofG(m), G%level, .false., G%nvars, G%shape, G%ctx)
       call G%encap%create(CofF(m), G%level, .false., G%nvars, G%shape, G%ctx)
    end do

    call G%encap%create(tmp, G%level, .false., G%nvars, G%shape, G%ctx)

    do m = 1, F%nnodes-1
       call F%encap%create(FofF(m), F%level, .false., F%nvars, F%shape, F%ctx)
    end do

    !!!! restrict qs and recompute fs
    tm = t0 + dt*G%nodes
    do m = 1, G%nnodes
       ! XXX: use rmat here...
       call F%restrict(F%qSDC(trat*(m-1)+1), G%qSDC(m), F%level, F%ctx, G%level, G%ctx)
       call G%sweeper%evaluate(G, tm(m), m)
    end do

    !!!! bring down fas correction from level above
    do m = 1, G%nnodes-1
       call G%encap%setval(G%tau(m), 0.0_pfdp)
    end do

    if (associated(F%tau)) then
       ! restrict fine fas corrections and sum between coarse nodes
       call G%encap%setval(tmp, 0.0_pfdp) ! needed for amr
       do m = 1, F%nnodes-1
          mc = int(ceiling(1.0_pfdp*m/trat))
          call F%restrict(F%tau(m), tmp, F%level, F%ctx, G%level, G%ctx)
          call G%encap%axpy(G%tau(mc), 1.0_pfdp, tmp)
       end do
    end if

    !!!! fas correction
    call G%sweeper%integrate(G, G%qSDC, G%fSDC, dt, CofG)
    call F%sweeper%integrate(F, F%qSDC, F%fSDC, dt, FofF)

    do m = 1, G%nnodes-1
       call G%encap%setval(CofF(m), 0.0_pfdp)
    end do

    ! restrict fine function values and sum between coarse nodes
    call G%encap%setval(tmp, 0.0_pfdp)
    do m = 1, F%nnodes-1
       mc = int(ceiling(1.0_pfdp*m/trat))
       call F%restrict(FofF(m), tmp, F%level, F%ctx, G%level, G%ctx)
       call G%encap%axpy(CofF(mc), 1.0_pfdp, tmp)
    end do

    do m = 1, G%nnodes-1
       call G%encap%axpy(G%tau(m),  1.0_pfdp, CofF(m))
       call G%encap%axpy(G%tau(m), -1.0_pfdp, CofG(m))
    end do

    !!!! destroy workspaces
    do m = 1, G%nnodes-1
       call G%encap%destroy(CofG(m))
       call G%encap%destroy(CofF(m))
    end do

    call G%encap%destroy(tmp)

    do m = 1, F%nnodes-1
       call F%encap%destroy(FofF(m))
    end do

    call end_timer(pf, TRESTRICT + F%level - 1)

  end subroutine restrict_time_space_fas

end module pf_mod_restrict

