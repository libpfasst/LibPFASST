!
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
!

! ------------------------------------------------------ module `feval`: start
! PFASST: RHS routines
module feval
	! For, e.g., C-pointers
	use iso_c_binding 

	! For LIBPFASST
	use pfasst

	! Problem parameters
	use probin
	! Encapsulation of data structure
	use encap

	! Variables starting with `i`, `j`, `k`, `l`, `m`, or `n` represent integers
	implicit none

contains
	! ------------------------------------------------------ subroutine `f1eval`: start
	! Explicit time-integration
	subroutine f1eval(Cptr2_sol_dat, t, lvl, Cptr2_ctx, Cptr2_f1_dat)
		! C-pointer to solution `sol`
		type(c_ptr), intent(in), value :: Cptr2_sol_dat
		real(pfdp), intent(in) :: t
		integer(c_int), intent(in) :: lvl
		type(c_ptr), intent(in), value :: Cptr2_ctx
		! C-pointer to `f1(t, sol)`
		type(c_ptr), intent(in), value :: Cptr2_f1_dat

		! Let `f1[i](t, sol) = -2*i*sol[i](t, x, y, z)`

		real(pfdp), pointer :: sol_dat_u(:, :, :, :)
		integer :: sol_dat_nfields, sol_dat_nx, sol_dat_ny, sol_dat_nz

		real(pfdp), pointer :: f1_dat_u(:, :, :, :)
		integer :: sol_dat_field

		sol_dat_u => get_u(Cptr2_sol_dat)
		sol_dat_nfields = get_nfields(Cptr2_sol_dat)
		sol_dat_nx = get_nx(Cptr2_sol_dat)
		sol_dat_ny = get_ny(Cptr2_sol_dat)
		sol_dat_nz = get_nz(Cptr2_sol_dat)

		f1_dat_u => get_u(Cptr2_f1_dat)

		! Dirichlet values everywhere
		f1_dat_u = 0.0_pfdp

		! Bulk values
		do sol_dat_field = 1, sol_dat_nfields
			f1_dat_u(sol_dat_field, 2:sol_dat_nx-1, 2:sol_dat_ny-1, 2:sol_dat_nz-1) = &
				-2.0_pfdp*dble(sol_dat_field)*sol_dat_u(sol_dat_field, 2:sol_dat_nx-1, 2:sol_dat_ny-1, 2:sol_dat_nz-1)
		end do
	end subroutine f1eval
	! ------------------------------------------------------ subroutine `f1eval`: stop

	! ------------------------------------------------------ subroutine `f2eval`: start
	! Implicit time-integration
	subroutine f2eval(Cptr2_sol_dat, t, lvl, Cptr2_ctx, Cptr2_f2_dat)
		! C-pointer to solution `sol`
		type(c_ptr), intent(in), value :: Cptr2_sol_dat
		real(pfdp), intent(in) :: t
		integer(c_int), intent(in) :: lvl
		type(c_ptr), intent(in), value :: Cptr2_ctx
		! C-pointer to `f2(t, sol)`
		type(c_ptr), intent(in), value :: Cptr2_f2_dat

		! Let `f2[i](t, sol) = -4*i*sol[i](t, x, y, z)`

		real(pfdp), pointer :: sol_dat_u(:, :, :, :)
		integer :: sol_dat_nfields, sol_dat_nx, sol_dat_ny, sol_dat_nz

		real(pfdp), pointer :: f2_dat_u(:, :, :, :)
		integer :: sol_dat_field

		sol_dat_u => get_u(Cptr2_sol_dat)
		sol_dat_nfields = get_nfields(Cptr2_sol_dat)
		sol_dat_nx = get_nx(Cptr2_sol_dat)
		sol_dat_ny = get_ny(Cptr2_sol_dat)
		sol_dat_nz = get_nz(Cptr2_sol_dat)

		f2_dat_u => get_u(Cptr2_f2_dat)

		! Dirichlet values everywhere
		f2_dat_u = 0.0_pfdp

		! Bulk values
		do sol_dat_field = 1, sol_dat_nfields
			f2_dat_u(sol_dat_field, 2:sol_dat_nx-1, 2:sol_dat_ny-1, 2:sol_dat_nz-1) = &
				-4.0_pfdp*dble(sol_dat_field)*sol_dat_u(sol_dat_field, 2:sol_dat_nx-1, 2:sol_dat_ny-1, 2:sol_dat_nz-1)
		end do
	end subroutine f2eval
	! ------------------------------------------------------ subroutine `f2eval`: stop

	! ------------------------------------------------------ subroutine `f2comp`: start
	! Source evaluation for implicit time-integration
	subroutine f2comp(Cptr2_sol_dat, t, dt, Cptr2_rhs_dat, lvl, Cptr2_ctx, Cptr2_f2_dat)
		! C-pointer to `sol`
		type(c_ptr), intent(in), value :: Cptr2_sol_dat
		real(pfdp), intent(in) :: t, dt
		! C-pointer to RHS `rhs` in `sol-dt*f2(t, sol) = rhs(t, sol)`
		type(c_ptr), intent(in), value :: Cptr2_rhs_dat
		integer(c_int), intent(in) :: lvl
		type(c_ptr), intent(in), value :: Cptr2_ctx
		! C-pointer to `f2(t, sol)`
		type(c_ptr), intent(in), value :: Cptr2_f2_dat

		! Solve `sol-dt*f2(t, sol) = rhs(t, sol)` for `f2(t, sol)`

		integer :: sol_dat_nfields, sol_dat_nx, sol_dat_ny, sol_dat_nz
		real(pfdp), pointer :: sol_dat_u(:, :, :, :)

		real(pfdp), pointer :: rhs_dat_u(:, :, :, :)
		real(pfdp), pointer :: f2_dat_u(:, :, :, :)
		integer :: sol_dat_field

		sol_dat_u => get_u(Cptr2_sol_dat)
		sol_dat_nfields = get_nfields(Cptr2_sol_dat)
		sol_dat_nx = get_nx(Cptr2_sol_dat)
		sol_dat_ny = get_ny(Cptr2_sol_dat)
		sol_dat_nz = get_nz(Cptr2_sol_dat)

		rhs_dat_u => get_u(Cptr2_rhs_dat)

		f2_dat_u => get_u(Cptr2_f2_dat)

		! Dirichlet values everywhere
		f2_dat_u = 0.0_pfdp

		! Bulk values
		do sol_dat_field = 1, sol_dat_nfields
			! Solve `u[i](t(n+1), x, y, z)-dt*f2[i](t(n+1), sol) = rhs[i](t(n), x, y, z)` for `u[i](t(n+1), x, y, z)` given `f2[i](t(n+1), sol)`
			sol_dat_u(sol_dat_field, 2:sol_dat_nx-1, 2:sol_dat_ny-1, 2:sol_dat_nz-1) = &
				rhs_dat_u(sol_dat_field, 2:sol_dat_nx-1, 2:sol_dat_ny-1, 2:sol_dat_nz-1)/(1.0_pfdp+4.0_pfdp*dble(sol_dat_field)*dt)
		end do

		! We have `f2[i](t(n+1), sol) = -4*i*sol[i](t(n+1), x, y, z)`
		call f2eval(Cptr2_sol_dat, t, lvl, Cptr2_ctx, Cptr2_f2_dat)
	end subroutine f2comp
	! ------------------------------------------------------ subroutine `f2comp`: stop
end module feval
! ------------------------------------------------------ module `feval`: stop
