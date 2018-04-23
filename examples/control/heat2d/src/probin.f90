module probin
  use pf_mod_dtype

  double precision, parameter :: pi     = 3.141592653589793d0
  double precision, parameter :: two_pi = 6.2831853071795862d0

  integer, parameter :: maxlevs = 3

  integer, save :: wtype

  character(len=64), save :: problem_type, window_type

  double precision, save :: v      ! advection velocity
  double precision, save :: Lx     ! domain size x
  double precision, save :: Ly     ! domain size y
  double precision, save :: nu     ! viscosity
  double precision, save :: t0     ! initial time for exact solution
  double precision, save :: sigma  ! initial condition parameter
  integer,          save :: kfreq  ! initial condition parameter
  double precision, save :: dt     ! time step
  double precision, save :: Tfin   ! Final time

  double precision, save :: abs_res_tol ! absolute residual tolerance
  double precision, save :: rel_res_tol ! relative residual tolerance
  
  ! optimization problem parameters
  double precision, save :: alpha        ! regularization parameter
  double precision, save :: tol_grad     ! gradient tolerance, stop optimization if gradient norm below
  double precision, save :: tol_obj      ! objective function tolerance, stop optimization if objective function value below
  integer,          save :: max_opt_iter ! maximum number of optimization iterations
  

  integer, save :: ndim            ! number of dimensions
  integer, save :: nnodes(maxlevs) ! number of nodes
  integer, save :: nvars(maxlevs)  ! number of grid points
  integer, save :: nprob           ! which problem
  integer, save :: nsteps          ! number of time steps
  logical, save :: Finterp
  integer, save :: spatial_order   ! spatial order for operators
  integer, save :: interp_order
  integer, save ::  mg_interp_order
  integer, save ::  do_spec
  integer, save ::  N_Vcycles
  integer, save ::  Nrelax
  integer, save ::  do_imex         ! set to 1 to use imex sweeper, otherwise misdc
  integer, save ::  warmstart       ! set to 1 to use previous solution as initialization in state solve
  integer, save :: do_mixed         ! set to 1 to sweep on state and adjoint simultaneously (but without communication in adjoint)
  integer, save :: nsweeps(maxlevs) !  Sweeps at each levels
  integer, save :: nsweeps_pred(maxlevs)   !  Sweeps during predictor

  logical, save :: solve_y
  
  character(len=32), save :: pfasst_nml
  character(len=20), save :: fbase   !  base name for run
  character(len=64), save :: foutbase   !  base name for output file
  character(len=44) :: foutfoo          !  temp used to create foutbase
  character(len=128), save :: logfile    !  file to use for output of optimization progress
  integer, save    :: poutmod

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/ Finterp, ndim, nnodes, nvars,nprob, nsteps
  namelist /params/ spatial_order,interp_order, mg_interp_order, do_spec, N_Vcycles,Nrelax
  namelist /params/ pfasst_nml,fbase ,poutmod
  namelist /params/  abs_res_tol, rel_res_tol 
  namelist /params/  v, nu, t0, dt, Tfin,sigma, kfreq, Lx, Ly, alpha, max_opt_iter
  namelist /params/  do_imex, warmstart, do_mixed, logfile, nsweeps, nsweeps_pred

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename
    integer :: i
    CHARACTER(len=32) :: arg
    integer :: un


    !
    ! defaults
    !

    Finterp = .FALSE.
    ndim     = 2
    nnodes  = [ 2, 3, 3 ]

    nsteps  = -1

    v       = 0.0_pfdp
    Lx      = 1._pfdp
    Ly      = 1._pfdp
    nu      = 1._pfdp
    sigma   = 0.004_pfdp
    kfreq   = 1
    t0      = 0.0_pfdp
    dt      = 0.01_pfdp
    Tfin    = 1.0_pfdp

    abs_res_tol = 0.0
    rel_res_tol = 0.0

    spatial_order=2
    interp_order = 2

    do_spec = 1
    N_Vcycles = 1
    Nrelax = 1
    mg_interp_order = 2
    
    do_imex = 1
    warmstart = 0
    do_mixed = 0   
    
    max_opt_iter = 100
    alpha = 0.05
    tol_grad = 1e-6
    tol_obj  = 1e-6
 
    poutmod = 1

    logfile = "progress.log"

    nsweeps =  [1, 1, 1]
    nsweeps_pred = [1, 1, 1]
    !
    ! read
    !
    !  Read in stuff from a file
    un = 9
    write(*,*) 'opening file ',TRIM(filename), '  for input'
    open(unit=un, file = filename, status = 'old', action = 'read')
    read(unit=un, nml = params)
    close(unit=un)
          
    i = 0
    DO
       CALL get_command_argument(i, arg)
       IF (LEN_TRIM(arg) == 0) EXIT
       if (i > 0) then
          istring="&PARAMS "//TRIM(arg)//" /"    
          READ(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
       end if
       i = i+1
    END DO

    !  Reset dt if Tfin is set
    if (Tfin .gt. 0.0) dt = Tfin/nsteps

    !
    ! init
    !

    select case (window_type)
    case ("block")
       wtype = PF_WINDOW_BLOCK
    case ("ring")
       wtype = PF_WINDOW_RING
    end select

    output = "Dat/pfasst_V/numpy"

  end subroutine probin_init

end module probin

