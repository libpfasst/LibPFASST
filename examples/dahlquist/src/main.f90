!
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
!

! ------------------------------------------------------ program `main`: start
! PFASST: Main program for Dahlquist equation
program main
	! For, e.g., C-pointers
	use iso_c_binding 

	! For LIBPFASST
	use pfasst
	! For MPI
	use pf_mod_mpi, only: MPI_COMM_WORLD

	! Problem parameters
	use probin
	! Encapsulation of data structure
	use encap
	! Function evaluations
	use feval
	! Restriction and interpolation
	use transfer
	! Diagnostic routines
	use hooks

	! Variables starting with `i`, `j`, `k`, `l`, `m`, or `n` represent integers
	implicit none

	! Input file name
	character(len = 512) :: ifname

	! MPI error
	integer :: err

	! MPI rank
	integer :: rank

	! LIBPFASST data container defined in `libpfasst/src/pf_dtype_t.f90`
	type(pf_pfasst_t) :: pf

	! Communicator defined in `libpfasst/src/pf_dtype_t.f90`
	type(pf_comm_t) :: comm

	! Encapsulated data defined in `pfasst-cart/src/encap.f90` (here, ultimately, a pointer destination)
	type(pf_encap_t), target :: pf_encap

	! Timing
	real t_start, t_stop

	! Counter for levels
	integer :: lvl

	! Initial data
	type(c_ptr) :: Cptr2_ini_dat
	type(dat), pointer :: Fptr2_ini_dat
	real(pfdp), pointer :: ini_dat_u(:, :, :, :)
	real(pfdp) :: ini_dat_x, ini_dat_y, ini_dat_z
	integer :: i_ini_dat_field, i_ini_dat_x, i_ini_dat_y, i_ini_dat_z
	real(pfdp), parameter :: pi_pfdp = 3.141592653589793_pfdp
	real(pfdp), parameter :: x_char_pfdp = 2.0_pfdp
	real(pfdp), parameter :: y_char_pfdp = 3.0_pfdp
	real(pfdp), parameter :: z_char_pfdp = 4.0_pfdp

	! Final data
	type(c_ptr) :: Cptr2_end_dat

	!
	! Initialization of input
	!

	if(command_argument_count() >= 1) then
		 call get_command_argument(1, value = ifname)
	else
		 ifname = 'probin.nml'
	end if

	call probin_init(ifname)

	!
	! Initialization of MPI
	!

	call mpi_init(err)
	if(err .ne. 0) stop 'Cannot initialize MPI'
	call mpi_comm_rank(MPI_COMM_WORLD, rank, err)

	!
	! Initialization of LIBPFASST
	!

	! See `libpfasst/src/pf_mpi.f90`
	call pf_mpi_create(comm, MPI_COMM_WORLD)

	! See `libpfasst/src/pf_pfasst.f90`
	call pf_pfasst_create(pf, comm, fname = ifname)

	!
	! Initialization of user data encapsulation
	!

	call create_encap(pf_encap)

	! Set level parameters
	do lvl = 1, pf%nlevels
		! Create `shape` array
		allocate(pf%levels(lvl)%shape(4))

		! Set `shape` values
		pf%levels(lvl)%shape(1) = nfields(lvl)
		pf%levels(lvl)%shape(2) = nx(lvl)
		pf%levels(lvl)%shape(3) = ny(lvl)
		pf%levels(lvl)%shape(4) = nz(lvl)

		! Total number of variables in solution
		pf%levels(lvl)%nvars = nfields(lvl)*nx(lvl)*ny(lvl)*nz(lvl)

		! Number of nodes on each level
		pf%levels(lvl)%nnodes = nnodes(lvl)

		! Number of sweeps at each level
		pf%levels(lvl)%nsweeps = nsweeps(lvl)

		! Number of sweeps at each level during prediction
		pf%levels(lvl)%nsweeps_pred = nsweeps_pred(lvl)

		! Assign interpolation routine
		pf%levels(lvl)%interpolate => interpolate

		! Assing restriction routine
		pf%levels(lvl)%restrict => restrict

		! Assign instance of PFASST encapsulation
		pf%levels(lvl)%encap => pf_encap

		! Create solver/sweeper
		call pf_explicitQ_create(pf%levels(lvl)%sweeper, f1eval)
		!call pf_imexQ_create(pf%levels(lvl)%sweeper, f1eval, f2eval, f2comp)
	end do

	!
	! Define initial data
	!

	call c_f_pointer(Cptr2_ini_dat, Fptr2_ini_dat)
	call create_dat(Cptr2_ini_dat, pf%nlevels, 1, pf%levels(pf%nlevels)%nvars, pf%levels(pf%nlevels)%shape, pf%levels(pf%nlevels)%ctx)
	ini_dat_u => get_u(Cptr2_ini_dat)

	do i_ini_dat_field = 1, nfields(pf%nlevels)
		do i_ini_dat_x = 1, nx(pf%nlevels)
			ini_dat_x = (i_ini_dat_x-1)*x_char_pfdp/(nx(pf%nlevels)-1.0_pfdp)
			do i_ini_dat_y = 1, ny(pf%nlevels)
				ini_dat_y = (i_ini_dat_y-1)*y_char_pfdp/(ny(pf%nlevels)-1.0_pfdp)
				do i_ini_dat_z = 1, nz(pf%nlevels)
					ini_dat_z = (i_ini_dat_z-1)*z_char_pfdp/(nz(pf%nlevels)-1.0_pfdp)
					ini_dat_u(i_ini_dat_field, i_ini_dat_x, i_ini_dat_y, i_ini_dat_z) = &
						sin(pi_pfdp*ini_dat_x/x_char_pfdp)*&
						sin(pi_pfdp*ini_dat_y/y_char_pfdp)*&
						sin(pi_pfdp*ini_dat_z/z_char_pfdp)/(6.0_pfdp*i_ini_dat_field)
				end do
			end do
		end do
	end do

	!
	! Create encapsulation for final solution
	!

	call create_dat(Cptr2_end_dat, pf%nlevels, 1, pf%levels(pf%nlevels)%nvars, pf%levels(pf%nlevels)%shape, pf%levels(pf%nlevels)%ctx)

	!
	! PFASST communicator setup
	!

	! See `libpfasst/src/pf_mpi.f90`
	call pf_mpi_setup(comm, pf)

	!
	! PFASST object setup (for time-interpolation matrices)
	!

	! See `libpfasst/src/pf_pfasst.f90`
	call pf_pfasst_setup(pf)

	!
	! Add hook(s)
	!

	call pf_add_hook(pf, -1, PF_POST_SWEEP, print_err)

	!
	! Determine number of time-steps per processor
	!
	
	! Perform one time-step per processor and temporal subinterval if `nsteps < 0`
	if(nsteps < 0) then
		nsteps = comm%nproc
	end if

	!
	! All processes shall start the integration simultaneously
	!

	call mpi_barrier(MPI_COMM_WORLD, err)

	!
	! Start ticking
	!

	call cpu_time(t_start)

	!
	! Run PFASST
	!

	! Execute PFASST solver (see `libpfasst/src/pf_parallel.f90` for STD output)
	call pf_pfasst_run(pf, Cptr2_ini_dat, dt, 0.0_pfdp, nsteps, Cptr2_end_dat)

	!
	! Stop ticking
	!

	call cpu_time(t_stop)
end program main
! ------------------------------------------------------ program `main`: stop
