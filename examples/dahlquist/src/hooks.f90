!
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
!

! ------------------------------------------------------ module `hooks`: start
! PFASST: Post-processing routines
module hooks
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
	! ------------------------------------------------------ subroutine `print_err`: start
	! Print error after each sweep
	subroutine print_err(pf, lvl, state, ctx)
		type(pf_pfasst_t), intent(inout) :: pf
		type(pf_level_t), intent(inout) :: lvl
		type(pf_state_t), intent(in) :: state
		type(c_ptr), intent(in) :: ctx

		! Time
		real(pfdp) :: t_start, t_stop

		! PFASST
		type(c_ptr) :: Cptr2_pfasst_dat
		real(pfdp), pointer :: pfasst_dat_u(:, :, :, :)

		! Exact
		type(c_ptr) :: Cptr2_exact_start_dat
		real(pfdp), pointer :: exact_start_dat_u(:, :, :, :)
		real(pfdp) :: exact_start_dat_x, exact_start_dat_y, exact_start_dat_z
		integer :: i_exact_start_dat_field, i_exact_start_dat_x, i_exact_start_dat_y, i_exact_start_dat_z

		type(c_ptr) :: Cptr2_exact_stop_dat
		real(pfdp), pointer :: exact_stop_dat_u(:, :, :, :)
		real(pfdp) :: exact_stop_dat_x, exact_stop_dat_y, exact_stop_dat_z
		integer :: i_exact_stop_dat_field, i_exact_stop_dat_x, i_exact_stop_dat_y, i_exact_stop_dat_z

		real(pfdp), parameter :: pi_pfdp = 3.141592653589793_pfdp
		real(pfdp), parameter :: x_char_pfdp = 2.0_pfdp
		real(pfdp), parameter :: y_char_pfdp = 3.0_pfdp
		real(pfdp), parameter :: z_char_pfdp = 4.0_pfdp

		! Error
		real(pfdp) :: err_start, err_stop

		! Time
		t_start = state%t0
		t_stop = state%t0+state%dt

		! PFASST
		Cptr2_pfasst_dat = lvl%Q(lvl%nnodes)
		pfasst_dat_u => get_u(Cptr2_pfasst_dat)

		! Exact
		call create_dat(Cptr2_exact_start_dat, pf%nlevels, 1, pf%levels(pf%nlevels)%nvars, pf%levels(pf%nlevels)%shape, ctx)
		exact_start_dat_u => get_u(Cptr2_exact_start_dat)
		do i_exact_start_dat_field = 1, nfields(pf%nlevels)
			do i_exact_start_dat_x = 1, nx(pf%nlevels)
				exact_start_dat_x = (i_exact_start_dat_x-1)*x_char_pfdp/(nx(pf%nlevels)-1.0_pfdp)
				do i_exact_start_dat_y = 1, ny(pf%nlevels)
					exact_start_dat_y = (i_exact_start_dat_y-1)*y_char_pfdp/(ny(pf%nlevels)-1.0_pfdp)
					do i_exact_start_dat_z = 1, nz(pf%nlevels)
						exact_start_dat_z = (i_exact_start_dat_z-1)*z_char_pfdp/(nz(pf%nlevels)-1.0_pfdp)
						exact_start_dat_u(i_exact_start_dat_field, i_exact_start_dat_x, i_exact_start_dat_y, i_exact_start_dat_z) = &
							exp(-6.0_pfdp*i_exact_start_dat_field*t_start)*&
							sin(pi_pfdp*exact_start_dat_x/x_char_pfdp)*&
							sin(pi_pfdp*exact_start_dat_y/y_char_pfdp)*&
							sin(pi_pfdp*exact_start_dat_z/z_char_pfdp)/(6.0_pfdp*i_exact_start_dat_field)
					end do
				end do
			end do
		end do

		call create_dat(Cptr2_exact_stop_dat, pf%nlevels, 1, pf%levels(pf%nlevels)%nvars, pf%levels(pf%nlevels)%shape, ctx)
		exact_stop_dat_u => get_u(Cptr2_exact_stop_dat)
		do i_exact_stop_dat_field = 1, nfields(pf%nlevels)
			do i_exact_stop_dat_x = 1, nx(pf%nlevels)
				exact_stop_dat_x = (i_exact_stop_dat_x-1)*x_char_pfdp/(nx(pf%nlevels)-1.0_pfdp)
				do i_exact_stop_dat_y = 1, ny(pf%nlevels)
					exact_stop_dat_y = (i_exact_stop_dat_y-1)*y_char_pfdp/(ny(pf%nlevels)-1.0_pfdp)
					do i_exact_stop_dat_z = 1, nz(pf%nlevels)
						exact_stop_dat_z = (i_exact_stop_dat_z-1)*z_char_pfdp/(nz(pf%nlevels)-1.0_pfdp)
						exact_stop_dat_u(i_exact_stop_dat_field, i_exact_stop_dat_x, i_exact_stop_dat_y, i_exact_stop_dat_z) = &
							exp(-6.0_pfdp*i_exact_stop_dat_field*t_stop)*&
							sin(pi_pfdp*exact_stop_dat_x/x_char_pfdp)*&
							sin(pi_pfdp*exact_stop_dat_y/y_char_pfdp)*&
							sin(pi_pfdp*exact_stop_dat_z/z_char_pfdp)/(6.0_pfdp*i_exact_stop_dat_field)
					end do
				end do
			end do
		end do

		! Error
		err_start = maxval(abs(pfasst_dat_u-exact_start_dat_u))
		err_stop = maxval(abs(pfasst_dat_u-exact_stop_dat_u))

		! Write maximum norm difference
		if(state%iter == 1) then
			write(*, '(a, 1x, a, 1x, a, 1x, a, 1x, a, 1x, a)') '# step', '# iter', '# t_start', '# err_start', '# t_stop', '# err_stop'
		end if
		write(*, '(1i0, 1x, 1i0, 1x, 1e10.4, 1x, 1e10.4, 1x, 1e10.4)') state%step, state%iter, t_start, t_stop, err_stop
	end subroutine
	! ------------------------------------------------------ subroutine `print_err`: stop
end module hooks
! ------------------------------------------------------ module `hooks`: stop
