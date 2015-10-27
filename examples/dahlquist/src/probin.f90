!
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
!

! ------------------------------------------------------ module `probin`: start
! PFASST: Input parameters
module probin
	! For, e.g., C-pointers
	use iso_c_binding

	! For LIBPFASST
	use pfasst

	! Variables starting with `i`, `j`, `k`, `l`, `m`, or `n` represent integers
	implicit none

	!
	! PFASST specific parameters
	!

	! See (1) `libpfasst/src/pf_dtype.f90` and (2) `libpfasst/src/pf_options.f90`

	integer, save :: nlevels
	integer, save :: niters
	integer, save :: qtype
	integer, save :: window
	real(pfdp), save :: abs_res_tol
	real(pfdp), save :: rel_res_tol
	logical, save :: Pipeline_g
	logical, save :: PFASST_pred
	logical, save :: echo_timings
	integer, save :: taui0
	character(len = 512), save :: outdir

	! List of variable names (name "pf_params" seems to be a must)
	namelist /pf_params/ nlevels, niters, qtype, window, abs_res_tol, rel_res_tol, Pipeline_g, PFASST_pred, echo_timings, taui0, outdir

	!
	! My (problem) specific parameters
	!

	! Maximum number of levels
	integer, parameter :: MAXLEVS = 3

	! Number of field variables
	integer, save :: nfields(MAXLEVS)
	! Number of grid points on each level
	integer, save :: nx(MAXLEVS)
	! Number of grid points on each level
	integer, save :: ny(MAXLEVS)
	! Number of grid points on each level
	integer, save :: nz(MAXLEVS)
	! Number of nodes on each level
	integer, save :: nnodes(MAXLEVS)
	! Sweeps at each levels
	integer, save :: nsweeps(MAXLEVS)
	! Sweeps during prediction
	integer, save :: nsweeps_pred(MAXLEVS)
	! Decides if interpolation of function values is done
	integer, save :: Finterp
	! Equal 1 if IMEX time stepping
	integer, save :: doIMEX
	! Number of time steps to take
	integer, save :: nsteps
	! Time step
	real(pfdp), save :: dt
	! Final time of run
	real(pfdp), save :: Tfin

	! List of variable names (name "my_params" seems to be up for grabs)
	namelist /my_params/ nfields, nx, ny, nz, nnodes, nsweeps, nsweeps_pred, Finterp, doIMEX, nsteps, dt, Tfin

contains
	! ------------------------------------------------------ subroutine `probin_init`: start
	! Parameter initialization
	subroutine probin_init(ifname)
		character(len = *), intent(in) :: ifname

		! I/O unit specifier (see [unit](http://h21007.www2.hp.com/portal/download/files/unprot/fortran/docs/lrm/lrm0355.htm#unit_spec))	
		integer, parameter :: un = 9
		! Read `i`-th command line argument (see [get_command_argument](https://gcc.gnu.org/onlinedocs/gfortran/GET_005fCOMMAND_005fARGUMENT.html))
		integer :: ird
		! Command line argument
		character(len = 512) :: arg
		! Parsed command line argument
		character(len = 512) :: parsed_arg
		! I/O status variable
		integer :: ios
		! I/O error message
		character(len = 512) :: err_msg

		!
		! PFASST specific parameters
		!

		! See (1) `libpfasst/src/pf_dtype.f90` and (2) `libpfasst/src/pf_options.f90`

		nlevels = 1
		niters = 5
		qtype = SDC_GAUSS_LOBATTO
		window = PF_WINDOW_BLOCK
		abs_res_tol = 0.0_pfdp
		rel_res_tol = 0.0_pfdp
		Pipeline_g = .false.
		PFASST_pred = .false.
		echo_timings = .false.
		taui0 = -999999
		outdir = ''

		!
		! My (problem) specific parameters
		!

		! Number of field variables
		nfields = 5
		! Number of grid points on each level
    nx = 3
		! Number of grid points on each level
    ny = 3
		! Number of grid points on each level
    nz = 3
		! Number of nodes on each level
		nnodes = 5
		! Sweeps at each levels
		nsweeps = 1
		! Sweeps during prediction
		nsweeps_pred = 1
		! Decides if interpolation of function values is done
    Finterp = 0
		! Equal 1 if IMEX time stepping
		doIMEX = 0
		! Number of time steps to take
    nsteps = 4
		! Time step
    dt = 0.25_pfdp
		! Final time of run
    Tfin = 1.0_pfdp

		! Read from `params.nml` file
		open(unit = un, file = ifname, status = 'old', action = 'read')
			read(unit = un, nml = pf_params)
			read(unit = un, nml = my_params)
		close(unit = un)

		! Read command line arguments to possibly overwrite content of `params.nml`
		ird = 0
		do
			call get_command_argument(ird, arg)
			if(len_trim(arg) == 0) exit
			if(ird > 0) then
				parsed_arg = "&PF_PARAMS "//trim(arg)//" /"    
				read(parsed_arg, nml = pf_params, iostat = ios, iomsg = err_msg)
				! Read namelist
				parsed_arg = "&MY_PARAMS "//trim(arg)//" /"    
				read(parsed_arg, nml = my_params, iostat = ios, iomsg = err_msg)
			end if
			ird = ird+1
		end do

		! Time-domain and step-size
    if(Tfin .gt. 0.0_pfdp) then
       dt = Tfin/nsteps
    else
       Tfin = nsteps*dt
    endif
	end subroutine probin_init
	! ------------------------------------------------------ subroutine `probin_init`: stop
end module probin
! ------------------------------------------------------ module `probin`: stop
