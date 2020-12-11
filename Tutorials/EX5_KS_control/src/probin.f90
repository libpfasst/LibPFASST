module probin
  use pf_mod_dtype

  double precision, parameter :: pi     = 3.141592653589793d0

  integer, parameter :: maxlevs = 3

  integer, save :: wtype

  character(len=64), save :: problem_type, window_type

  double precision, save :: Lx     ! domain size x
  double precision, save :: dt     ! time step
  double precision, save :: Tfin   ! Final time

  double precision, save :: abs_res_tol ! absolute residual tolerance
  double precision, save :: rel_res_tol ! relative residual tolerance
  
  ! optimization problem parameters
  double precision, save :: alpha        ! regularization parameter
  double precision, save :: tol_grad     ! gradient tolerance, stop optimization if gradient norm below
  double precision, save :: tol_obj      ! objective function tolerance, stop optimization if objective function value below
  integer,          save :: max_opt_iter ! maximum number of optimization iterations
  double precision, save :: max_step_size ! maximum/initial step size for line search (Armijo only; in WP this gets increased); use for experimenting with warm starts
  

  integer, save :: ndim            ! number of dimensions
  integer, save :: nnodes(maxlevs) ! number of nodes
  integer, save :: nvars(maxlevs)  ! number of grid points
  integer, save :: nsteps          ! number of time steps
  logical, save :: Finterp
  integer, save :: do_imex          ! set to 1 to use imex sweeper, otherwise misdc
  integer, save :: warmstart        ! set to 1 to use previous solution as initialization in state solve
  integer, save :: nsweeps(maxlevs) ! Sweeps at each levels
  integer, save :: nsweeps_pred(maxlevs)   ! Sweeps during predictor
  integer, save :: test_no          ! to distinguish individual test runs
  integer, save :: opt_method       ! 0: steepest descent, 1: Polak-Ribiere nonlinear conjugate gradient
                                    ! 2: Dai-Yuan ncg      3: Fletcher-Reeves ncg
  logical, save :: use_wolfe        ! whether to use strong Wolfe stepsize (if false, use Armijo)
  logical, save :: write_numpy      ! whether to output numpy arrays of control, state, gradient
                                    ! turn off for timing experiments
  
  character(len=32), save :: pfasst_nml !  file for setting PFASST options
  character(len=20), save :: fbase      !  base name for run
  character(len=64), save :: foutbase   !  base name for output file
  character(len=44) :: foutfoo          !  temp used to create foutbase
  character(len=128), save :: logfile   !  file to use for output of optimization progress
  integer, save    :: poutmod

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring     ! stores command line argument
  CHARACTER(LEN=255) :: message     ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/ Finterp, ndim, nnodes, nvars, nsteps
  namelist /params/ pfasst_nml,fbase ,poutmod, output
  namelist /params/ abs_res_tol, rel_res_tol, tol_grad, tol_obj
  namelist /params/ dt, Tfin, Lx, alpha, max_opt_iter
  namelist /params/ do_imex, warmstart, logfile, nsweeps, nsweeps_pred
  namelist /params/ max_step_size, test_no, opt_method, use_wolfe, write_numpy


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
    ndim     = 1
    nnodes  = [ 2, 3, 3 ]
    nsteps  = -1

    Lx      = 1._pfdp
    dt      = 0.01_pfdp
    Tfin    = 1.0_pfdp

    abs_res_tol = 0.0
    rel_res_tol = 0.0

    do_imex = 1
    warmstart = 0
    max_opt_iter = 100
    alpha = 0.05
    tol_grad = 1e-6
    tol_obj  = 1e-6
    max_step_size = 1.0

    opt_method = 0
    use_wolfe = .false.
    write_numpy = .true.

 
    poutmod = 1

    test_no=0
    logfile = "progress.log"

    nsweeps =  [1, 1, 1]
    nsweeps_pred = [1, 1, 1]
    !
    ! read
    !
    !  Read in stuff from a file
    un = 9
    write(*,*) 'opening file ',TRIM(filename), '  for input'
    open(unit=un, file = filename, status = 'old', action = 'read') ! do that with MPI_File_Open?
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

    output = "Dat/pfasst_V/numpy"

  end subroutine probin_init

end module probin

