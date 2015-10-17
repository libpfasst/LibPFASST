!
! Copyright (c) 2015, Michael Minion and Andreas Kreienbuehl. All rights reserved.
!

! ------------------------------------------------------ module *feval*: start
! PFASST: RHS routines for compressible Navier-Stokes equations example (i.e. Cart)
module feval
	! For, e.g., C-pointers
	use iso_c_binding 

	! For LIBPFASST
	use pfasst

	! Problem parameters
	use probin
	! Encapsulation of data structure
	use encap

	! Variables starting with *i*, *j*, *k*, *l*, *m*, or *n* represent integers
	implicit none

contains
	! ------------------------------------------------------ subroutine *f1eval*: start
	! Explicit time-integration
	subroutine f1eval(Cptr2_u_dat, t, lvl, Cptr2_levctx, Cptr2_f1_dat)
		type(c_ptr), intent(in), value :: Cptr2_u_dat   ! pointer to solution u
		real(pfdp), intent(in) :: t
		integer(c_int), intent(in) :: lvl
		type(c_ptr), intent(in), value :: Cptr2_levctx   ! pointer to level context
		type(c_ptr), intent(in), value :: Cptr2_f1_dat   ! pointer to solution f(u)

		! Consider $\dot{u} = -f1(u)$ 

		! Determine number of grid points in x-, y-, and z-direction
		integer :: nx, ny, nz

		! Get data array associated to input data encapsulation
		real(pfdp), pointer :: u_dat_arr(:, :, :, :)

		! Define data array associated to output data encapsulation
		real(pfdp), pointer :: f1_dat_arr(:, :, :, :)

		real(pfdp)  :: tau  !local variable for Andreas
                integer :: k


		! Determine number of grid points in x-, y-, and z-direction
		nx = get_nx(Cptr2_u_dat)
		ny = get_ny(Cptr2_u_dat)
		nz = get_nz(Cptr2_u_dat)

		! Get data array associated to input data encapsulation
		u_dat_arr => get_arr(Cptr2_f1_dat)

		! Define data array associated to output data encapsulation
		f1_dat_arr => get_arr(Cptr2_f1_dat)

		! Evaluate *f1_dat_arr* given *u_dat_arr*

		! Set homogeneous Dirichlet boundary values
                f1_dat_arr=0.0_pfdp

		! Set bulk values and, therefore, have *f1_dat_arr = -tau*u_dat_arr*
                do k = 1,5
                   tau = dble(k)
                   f1_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1) = &
			-tau*u_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1)
                end do
	end subroutine f1eval
	! ------------------------------------------------------ subroutine *f1eval*: stop

	! ------------------------------------------------------ subroutine *f2eval*: start
	! Implicit time-integration
	subroutine f2eval(Cptr2_u_dat, t, lvl, Cptr2_levctx, Cptr2_f2_dat)
		type(c_ptr), intent(in), value :: Cptr2_u_dat   ! pointer to solution u
		real(pfdp), intent(in) :: t                     ! time
		integer(c_int), intent(in) :: lvl               ! level index
		type(c_ptr), intent(in), value :: Cptr2_levctx  ! pointer to level context
		type(c_ptr), intent(in), value :: Cptr2_f2_dat  ! pointer to f2(u)

		! Define data array associated to Cptr arguments
		real(pfdp), pointer :: u_dat_arr(:, :, :, :)
		real(pfdp), pointer :: f2_dat_arr(:, :, :, :)
		real(pfdp)  :: tau     !local variable for Andreas
                integer :: k           !  loop counter
		integer :: nx, ny,nz   ! size of the grids

		! Get fortran pointers 
		u_dat_arr => get_arr(Cptr2_u_dat)
		f2_dat_arr => get_arr(Cptr2_f2_dat)

		nx = get_nx(Cptr2_u_dat)
		ny = get_ny(Cptr2_u_dat)
		nz = get_nz(Cptr2_u_dat)

		! Set variables
                f2_dat_arr=0.0_pfdp
                do k = 1,5
                   tau = dble(k)
                   f2_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1) = &
			-2.0_pfdp*tau*u_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1)
                end do

	end subroutine f2eval
	! ------------------------------------------------------ subroutine *f2eval*: stop

	! ------------------------------------------------------ subroutine *f2comp*: start
	! Solve system of equations for implicit time-step
	subroutine f2comp(Cptr2_u_dat, t, dt, Cptr2_rhs_dat, lvl, Cptr2_levctx, Cptr2_f2_dat)
		type(c_ptr), intent(in), value :: Cptr2_u_dat   ! pointer to solution f(u)
		real(pfdp), intent(in) :: t, dt
		type(c_ptr), intent(in), value :: Cptr2_rhs_dat ! pointer to rhs
		integer(c_int), intent(in) :: lvl  
		type(c_ptr), intent(in), value :: Cptr2_levctx  ! pointer to level context
		type(c_ptr), intent(in), value :: Cptr2_f2_dat  !  solution f(u)

                ! Solve the equation u-dt*f_2(u,t) = rhs

		! Get data array associated to input data encapsulation
		real(pfdp), pointer :: u_dat_arr(:, :, :, :)
		real(pfdp), pointer :: rhs_dat_arr(:, :, :, :)
		real(pfdp), pointer :: f2_dat_arr(:, :, :, :)
            
		real(pfdp)  :: tau  !local variable for Andreas
		integer :: nx, ny,nz   ! size of the grids
                integer :: k           !  loop counter

		u_dat_arr => get_arr(Cptr2_u_dat)
		f2_dat_arr => get_arr(Cptr2_f2_dat)
		rhs_dat_arr => get_arr(Cptr2_rhs_dat)

		nx = get_nx(Cptr2_u_dat)
		ny = get_ny(Cptr2_u_dat)
		nz = get_nz(Cptr2_u_dat)

                do k = 1,5
                   tau = dble(k)
                    u_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1) = &
			rhs_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1)/(1.0_pfdp+2.0_pfdp*tau)
                    f2_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1) = &
			-2.0_pfdp*tau*u_dat_arr(k, 2:nx-1, 2:ny-1, 2:nz-1)
                end do

              !  f2_dat_arr=(u_dat_arr-rhs_dat_arr)/dt  !!!  seems easy, but not recommended

	end subroutine f2comp
	! ------------------------------------------------------ subroutine *f2comp*: stop
end module feval
! ------------------------------------------------------ module *feval*: stop
