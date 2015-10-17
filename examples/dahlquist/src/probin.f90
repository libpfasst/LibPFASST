!
! Copyright (c) 2015, Michael Minion and Andreas Kreienbuehl. All rights reserved.
!

! ------------------------------------------------------ module *probin*: start
! PFASST: Input parameters (initialized via *probin_init*)
module probin
	! For, e.g., C-pointers
	use iso_c_binding

	! For LIBPFASST
	use pfasst

	! Variables starting with *i*, *j*, *k*, *l*, *m*, or *n* represent integers
	implicit none

	! PFASST

	! See (1) *libpfasst/src/pf_dtype.f90* and (2) *libpfasst/src/pf_options.f90*
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
	character(len = 256), save :: outdir

!	! Maximum number of levels
	integer, parameter :: MAXLEVS = 3
!	! Decides if interpolation of function values is done
	integer, save :: Finterp
!	! Number of dimensions
!	integer, save :: ndim
!	! Number of fake processors
!	integer, save :: nfake
!	! Number of nodes on each level
	integer, save :: nnodes(MAXLEVS)
!	! Define which problem to run
!	integer, save :: nprob
!	! To store the difference stencil type
!	integer, save :: nstencil
!	! Boundary condition
!	integer, save :: nbcnr
!	! Number of time steps to take
	integer, save :: nsteps
!	! Number of field variables
	integer, save :: nflds
!	! Number of grid points on each level
	integer, save :: nx(MAXLEVS)
!	! Number of grid points on each level
	integer, save :: ny(MAXLEVS)
!	! Number of grid points on each level
	integer, save :: nz(MAXLEVS)
!	! Equal 1 if IMEX time stepping
	integer, save :: doIMEX
!	! Time step
	real(pfdp), save :: dt
!	! Final time of run
	real(pfdp), save :: Tfin
!	! Sweeps at each levels
	integer, save :: nsweeps(MAXLEVS)
!	! Sweeps during predictor
	integer, save :: nsweeps_pred(MAXLEVS)



		namelist /params/ Finterp, nnodes, nsweeps, nsweeps_pred, nx, ny, nz, Tfin, nsteps, dt
	! MAIN (i.e. "Cart", see *tar/pfasst/cartInterface.f90* and *tar/pfasst/inviscidRHSunified.f90*)

!	! Number of field variables stored
!	integer, save :: nq
!	! Number of field variables to compute residuals for
!	integer, save :: nvar
!	! Ratio of specific heats
!	real(pfdp), save :: gamma
!	! Grid speeds in three directions
!	real(pfdp), save :: timeMetric(3)
!	! Coordinate 1 spacing
!	real(pfdp), save :: dx
!	! Coordinate 2 spacing
!	real(pfdp), save :: dy
!	! coordinate 3 spacing
!	real(pfdp), save :: dz
!	! Coordinate 1 dimension
!	integer, save :: jmax
!	! Coordinate 2 dimension
!	integer, save :: kmax
!	! Coordinate 3 dimension
!	integer, save :: lmax
!	! Order of physical flux
!	integer, save :: flux_order
!	! Order of dissipation
!	integer, save :: diss_order
!	! Scaling for dissipation
!	real(pfdp), save :: efac
!	! Storage type
!	character(len = 256), save :: istor

contains
	! ------------------------------------------------------ subroutine *probin_init*: start
	! Parameter initialization
	subroutine probin_init(in_fname)
		character(len = *), intent(in) :: in_fname

    integer,parameter :: un=9
    integer :: ird, ios
    CHARACTER(len=64) :: arg
    CHARACTER(LEN=64) :: istring   ! For pulling off command line options
    CHARACTER(LEN=64) :: message   ! use for I/O error messages

		namelist /params/ Finterp, nnodes, nsweeps, nsweeps_pred, nx, ny, nz, Tfin, nsteps, dt, nflds


		! PFASST (i.e. "libpfasst")

		namelist /params/ nlevels, niters, qtype, window, abs_res_tol, rel_res_tol, Pipeline_g, PFASST_pred, echo_timings, taui0, outdir

		! See (1) *libpfasst/src/pf_dtype.f90* and (2) *libpfasst/src/pf_options.f90*
		nlevels = -1
		niters = 5
		qtype = SDC_GAUSS_LOBATTO
		window = PF_WINDOW_BLOCK
		abs_res_tol = 0.d0
		rel_res_tol = 0.d0
		Pipeline_g = .false.
		PFASST_pred = .false.
		echo_timings = .false.
		taui0 = -999999
		outdir = ''

!		namelist /params/ Finterp, ndim, nfake, nnodes, nprob, nstencil, nbcnr, nsteps, nx, ny, nz, doIMEX, dt, Tfin, nsweeps, nsweeps_pred

!		Finterp
!		ndim
!		nfake
!		nnodes
!		nprob
!		nstencil
!		nbcnr
!		nsteps
!		nx
!		ny
!		nz
!		doIMEX
!		dt
!		Tfin
!		nsweeps
!		nsweeps_pred

		! MAIN (i.e. "Cart")

!		namelist /params/ nq, nvar, gamma, timeMetric, dx, dy, dz, jmax, kmax, lmax, flux_order, diss_order, efac, istor

!		nq = 9
!		nvar = 5
!		gamma
!		timeMetric(3)
!		dx
!		dy
!		dz
!		jmax = 60
!		kmax = 60
!		lmax = 60
!		flux_order = 2
!		diss_order = 4
!		efac
!		istor = 'row'


    Finterp = 0   !  default is 0
			nnodes =      5
			nsweeps = 1
			nsweeps_pred = 1
    nx =   8
    ny =   8
    nz =   8

    Tfin = 1.0
    nsteps = 4

    dt     = 0.25


		nflds = 5

    open(unit=un, file = in_fname, status = 'old', action = 'read')
    read(unit=un, nml = params)
    close(unit=un)

    !  Read command line
    ird = 0
    DO
       CALL get_command_argument(ird, arg)
       IF (LEN_TRIM(arg) == 0) EXIT
       if (ird > 0) then
          istring="&PARAMS "//TRIM(arg)//" /"    
          READ(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
          istring="&JBPARAMS "//TRIM(arg)//" /"    
       end if
       ird = ird+1
    END DO


    if( Tfin .gt. 0.0_pfdp )then
       dt=Tfin/nsteps
    else
       Tfin = dt*nsteps
    endif

	end subroutine probin_init
	! ------------------------------------------------------ subroutine *probin_init*: stop
end module probin
! ------------------------------------------------------ module *probin*: stop
