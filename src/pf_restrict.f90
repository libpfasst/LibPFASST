!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
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
    use pf_mod_sweep
    use pf_mod_timer
    use transfer, only: restrict

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: F, G

    integer    :: m, mc, trat
    real(pfdp) :: tm(G%nnodes)
    type(pf_encap_t) :: &
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
       call create(CofG(m), G%level, .false., G%nvars, G%shape, G%ctx)
       call create(CofF(m), G%level, .false., G%nvars, G%shape, G%ctx)
    end do

    call create(tmp, G%level, .false., G%nvars, G%shape, G%ctx)

    do m = 1, F%nnodes-1
       call create(FofF(m), F%level, .false., F%nvars, F%shape, F%ctx)
    end do

    !!!! restrict qs and recompute fs
    tm = t0 + dt*G%nodes
    do m = 1, G%nnodes
       ! XXX: use rmat here...
       call restrict(F%qSDC(trat*(m-1)+1), G%qSDC(m), F%level, F%ctx, G%level, G%ctx)
       call sdceval(tm(m), m, G)
    end do

    !!!! bring down fas correction from level above
    do m = 1, G%nnodes-1
       call setval(G%tau(m), 0.0_pfdp)
    end do

    if (associated(F%tau)) then
       ! restrict fine fas corrections and sum between coarse nodes
       call setval(tmp, 0.0_pfdp) ! needed for amr
       do m = 1, F%nnodes-1
          mc = int(ceiling(1.0_pfdp*m/trat))
          call restrict(F%tau(m), tmp, F%level, F%ctx, G%level, G%ctx)
          call axpy(G%tau(mc), 1.0_pfdp, tmp)
       end do
    end if

    !!!! fas correction
    call sdc_integrate(G%qSDC, G%fSDC, dt, G, CofG)
    call sdc_integrate(F%qSDC, F%fSDC, dt, F, FofF)

    do m = 1, G%nnodes-1
       call setval(CofF(m), 0.0_pfdp)
    end do

    ! restrict fine function values and sum between coarse nodes
    call setval(tmp, 0.0_pfdp)
    do m = 1, F%nnodes-1
       mc = int(ceiling(1.0_pfdp*m/trat))
       call restrict(FofF(m), tmp, F%level, F%ctx, G%level, G%ctx)
       call axpy(CofF(mc), 1.0_pfdp, tmp)
    end do

    do m = 1, G%nnodes-1
       call axpy(G%tau(m),  1.0_pfdp, CofF(m))
       call axpy(G%tau(m), -1.0_pfdp, CofG(m))
    end do

    !!!! destroy workspaces
    do m = 1, G%nnodes-1
       call destroy(CofG(m))
       call destroy(CofF(m))
    end do

    call destroy(tmp)

    do m = 1, F%nnodes-1
       call destroy(FofF(m))
    end do

    call end_timer(pf, TRESTRICT + F%level - 1)

  end subroutine restrict_time_space_fas

end module pf_mod_restrict

