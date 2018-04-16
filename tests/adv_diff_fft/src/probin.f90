module probin
  use pf_mod_dtype

  double precision, parameter :: pi     = 3.141592653589793d0
  double precision, parameter :: two_pi = 6.2831853071795862d0

  integer, parameter :: MAXLEVS = 4


  character(len=64), save :: problem_type, window_type

  double precision, save :: v      ! advection velocity
  double precision, save :: Lx     ! domain size
  double precision, save :: nu     ! viscosity
  double precision, save :: t00     ! initial time for exact solution
  double precision, save :: sigma  ! initial condition parameter
  integer,          save :: kfreq  ! initial condition parameter
  double precision, save :: dt     ! time step
  double precision, save :: Tfin   ! Final time
  integer, save    :: nsweeps(MAXLEVS) !  Sweeps at each levels
  integer, save    :: nsweeps_pred(MAXLEVS)   !  Sweeps during predictor


  integer, save :: nnodes(MAXLEVS) ! number of nodes
  integer, save :: nx(MAXLEVS)     ! number of grid points
  integer, save :: nprob           ! which problem
  integer, save :: nsteps          ! number of time steps
  logical, save :: Finterp
  logical, save :: use_LUq
  integer, save :: imex_stat

  character(len=32), save :: pfasst_nml

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/ Finterp, nnodes, nx,nprob, nsteps, dt, Tfin
  namelist /params/  pfasst_nml, nsweeps, nsweeps_pred
  namelist /params/  v, nu, t00, sigma, kfreq, use_LUq,imex_stat

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
    nnodes  = 3

    nsteps  = -1

    v       = 1.0_pfdp
    Lx      = 1._pfdp
    nu      = 0.1_pfdp
    sigma   = 0.004_pfdp
    kfreq   = 1
    t00      = 0.25_pfdp
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp

    use_LUq=.TRUE.
    imex_stat=2    !  Default is full IMEX
    
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

  end subroutine probin_init

end module probin
