!
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
!

! ------------------------------------------------------ module *transfer`: start
! PFASST: Restriction and interpolation routines
module transfer
	! For, e.g., C-pointers
	use iso_c_binding

	! For LIBPFASST
	use pfasst

	! Problem parameters
	use probin
	! Encapsulation of data structure
	use encap
	! Function evaluations
	use feval

	! Variables starting with *i*, *j*, *k*, *l*, *m*, or *n* represent integers
	implicit none

contains
	! ------------------------------------------------------ subroutine `interpolate`: start
	! Coarse to fine spatial grid interpolation
	subroutine interpolate(Cptr2_fine_dat, Cptr2_coarse_dat, fine_lvl, Cptr2_fine_ctx, coarse_lvl, Cptr2_coarse_ctx, t)
		type(c_ptr), intent(in), value :: Cptr2_fine_dat, Cptr2_coarse_dat
		integer, intent(in) :: fine_lvl
		type(c_ptr), intent(in), value :: Cptr2_fine_ctx
		integer, intent(in) :: coarse_lvl
		type(c_ptr), intent(in), value :: Cptr2_coarse_ctx
		real(pfdp), intent(in) :: t

		! Second-order interpolation from coarse to fine grid data

		! Determine number of grid points
		integer :: coarse_dat_nx, coarse_dat_ny, coarse_dat_nz
		integer :: fine_dat_nx, fine_dat_ny, fine_dat_nz

		! Get data array associated to input data encapsulation
		real(pfdp), pointer :: coarse_dat_u(:, :, :, :)

		! Define data array associated to output data encapsulation
		real(pfdp), pointer :: fine_dat_u(:, :, :, :)

		! Define data array associated to output data encapsulation
		coarse_dat_u => get_u(Cptr2_coarse_dat)

		! Get data array associated to input data encapsulation
		fine_dat_u => get_u(Cptr2_fine_dat)

		! Get number of coarse gird points
		coarse_dat_nx = get_nx(Cptr2_coarse_dat)
		coarse_dat_ny = get_ny(Cptr2_coarse_dat)
		coarse_dat_nz = get_nz(Cptr2_coarse_dat)

		! Get number of fine gird points
		fine_dat_nx = get_nx(Cptr2_fine_dat)
		fine_dat_ny = get_ny(Cptr2_fine_dat)
		fine_dat_nz = get_nz(Cptr2_fine_dat)

!		if(coarse_dat_nx == fine_dat_nx) then
!			! Coarse and fine grid may be identical
			fine_dat_u = coarse_dat_u
!		else
!			! Evaluate `fine_dat_u` given `coarse_dat_u`
!			! Point injection
!			fine_dat_u(:, 1:fine_dat_nx:2, 1:fine_dat_ny:2, 1:fine_dat_nz:2) = &
!				coarse_dat_u(:, 1:coarse_dat_nx, 1:coarse_dat_ny, 1:coarse_dat_nz)
!	
!			! Linear `x`-interpolation
!			fine_dat_u(:, 2:fine_dat_nx-1:2, 1:fine_dat_ny:2, 1:fine_dat_nz:2) = 0.5_pfdp*(&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 1:coarse_dat_ny, 1:coarse_dat_nz)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 1:coarse_dat_ny, 1:coarse_dat_nz))
!	
!			! Linear `y`-interpolation
!			fine_dat_u(:, 1:fine_dat_nx:2, 2:fine_dat_ny-1:2, 1:fine_dat_nz:2) = 0.5_pfdp*(&
!				coarse_dat_u(:, 1:coarse_dat_nx, 1:coarse_dat_ny-1, 1:coarse_dat_nz)+&
!				coarse_dat_u(:, 1:coarse_dat_nx, 2:coarse_dat_ny, 1:coarse_dat_nz))
!	
!			! Linear `z`-interpolation
!			fine_dat_u(:, 1:fine_dat_nx:2, 1:fine_dat_ny:2, 2:fine_dat_nz-1:2) = 0.5_pfdp*(&
!				coarse_dat_u(:, 1:coarse_dat_nx, 1:coarse_dat_ny, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 1:coarse_dat_nx, 1:coarse_dat_ny, 2:coarse_dat_nz))
!	
!			! Linear `xy`-interpolation
!			fine_dat_u(:, 2:fine_dat_nx-1:2, 2:fine_dat_ny-1:2, 1:fine_dat_nz:2) = 0.25_pfdp*(&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 1:coarse_dat_ny-1, 1:coarse_dat_nz)+&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 2:coarse_dat_ny, 1:coarse_dat_nz)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 1:coarse_dat_ny-1, 1:coarse_dat_nz)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 2:coarse_dat_ny, 1:coarse_dat_nz))
!	
!			! Linear `xz`-interpolation
!			fine_dat_u(:, 2:fine_dat_nx-1:2, 1:fine_dat_ny:2, 2:fine_dat_nz-1:2) = 0.25_pfdp*(&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 1:coarse_dat_ny, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 1:coarse_dat_ny, 2:coarse_dat_nz)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 1:coarse_dat_ny, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 1:coarse_dat_ny, 2:coarse_dat_nz))
!	
!			! Linear `yz`-interpolation
!			fine_dat_u(:, 1:fine_dat_nx:2, 2:fine_dat_ny-1:2, 2:fine_dat_nz-1:2) = 0.25_pfdp*(&
!				coarse_dat_u(:, 1:coarse_dat_nx, 1:coarse_dat_ny-1, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 1:coarse_dat_nx, 1:coarse_dat_ny-1, 2:coarse_dat_nz)+&
!				coarse_dat_u(:, 1:coarse_dat_nx, 2:coarse_dat_ny, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 1:coarse_dat_nx, 2:coarse_dat_ny, 2:coarse_dat_nz))
!	
!			! Linear `xyz`-interpolation
!			fine_dat_u(:, 2:fine_dat_nx-1:2, 2:fine_dat_ny-1:2, 2:fine_dat_nz-1:2) = 0.125_pfdp*(&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 1:coarse_dat_ny-1, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 1:coarse_dat_ny-1, 2:coarse_dat_nz)+&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 2:coarse_dat_ny, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 1:coarse_dat_nx-1, 2:coarse_dat_ny, 2:coarse_dat_nz)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 1:coarse_dat_ny-1, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 1:coarse_dat_ny-1, 2:coarse_dat_nz)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 2:coarse_dat_ny, 1:coarse_dat_nz-1)+&
!				coarse_dat_u(:, 2:coarse_dat_nx, 2:coarse_dat_ny, 2:coarse_dat_nz))
!		end if
	end subroutine interpolate
	! ------------------------------------------------------ subroutine `interpolate`: stop

	! ------------------------------------------------------ subroutine `restrict`: start
	! Fine to coarse spatial grid restriction
	subroutine restrict(Cptr2_fine_dat, Cptr2_coarse_dat, fine_lvl, Cptr2_fine_ctx, coarse_lvl, Cptr2_coarse_ctx, t)
		type(c_ptr), intent(in), value :: Cptr2_fine_dat, Cptr2_coarse_dat
		integer, intent(in) :: fine_lvl
		type(c_ptr), intent(in), value :: Cptr2_fine_ctx
		integer, intent(in) :: coarse_lvl
		type(c_ptr), intent(in), value :: Cptr2_coarse_ctx
		real(pfdp), intent(in) :: t

		! Point injection from coarse to fine grid data

		! Determine number of grid points in x-, y-, and z-direction
		integer :: fine_dat_nx, fine_dat_ny, fine_dat_nz
		integer :: coarse_dat_nx, coarse_dat_ny, coarse_dat_nz

		! Get data array associated to input data encapsulation
		real(pfdp), pointer :: fine_dat_u(:, :, :, :)

		! Define data array associated to output data encapsulation
		real(pfdp), pointer :: coarse_dat_u(:, :, :, :)

		! Get data array associated to input data encapsulation
		fine_dat_u => get_u(Cptr2_fine_dat)

		! Define data array associated to output data encapsulation
		coarse_dat_u => get_u(Cptr2_coarse_dat)

		! Get number of fine gird points
		fine_dat_nx = get_nx(Cptr2_fine_dat)
		fine_dat_ny = get_ny(Cptr2_fine_dat)
		fine_dat_nz = get_nz(Cptr2_fine_dat)

		! Get number of coarse gird points
		coarse_dat_nx = get_nx(Cptr2_coarse_dat)
		coarse_dat_ny = get_ny(Cptr2_coarse_dat)
		coarse_dat_nz = get_nz(Cptr2_coarse_dat)

!		if(fine_dat_nx == coarse_dat_nx) then
!			! Coarse and fine grid may be identical
			coarse_dat_u = fine_dat_u
!		else
!			! Evaluate `coarse_dat_u` given `fine_dat_u`
!			coarse_dat_u(:, 1:coarse_dat_nx, 1:coarse_dat_ny, 1:coarse_dat_nz) = &
!				fine_dat_u(:, 1:fine_dat_nx:2, 1:fine_dat_ny:2, 1:fine_dat_nz:2)
!		end if
	end subroutine restrict
	! ------------------------------------------------------ subroutine `restrict`: stop
end module transfer
! ------------------------------------------------------ module *transfer`: stop
