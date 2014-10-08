module probin
  use pf_mod_dtype

  double precision, parameter :: pi     = 3.141592653589793d0
  double precision, parameter :: two_pi = 6.2831853071795862d0

  integer, parameter :: maxlevs = 3

  integer, save :: wtype

  character(len=64), save :: problem_type, window_type

  double precision, save :: v      ! advection velocity
  double precision, save :: Lx     ! domain size
  double precision, save :: nu     ! viscosity
  double precision, save :: t0     ! initial time for exact solution
  double precision, save :: sigma  ! initial condition parameter
  integer,          save :: kfreq  ! initial condition parameter
  double precision, save :: dt     ! time step
  double precision, save :: Tfin   ! Final time

  double precision, save :: abs_tol ! absolute residual tolerance
  double precision, save :: rel_tol ! relative residual tolerance

  integer, save :: ndim             ! number of dimensions
  integer, save :: nlevs           ! number of pfasst levels
  integer, save :: nnodes(maxlevs) ! number of nodes
  integer, save :: nvars(maxlevs)  ! number of grid points
  integer, save :: nprob           ! which problem
  integer, save :: nsteps          ! number of time steps
  integer, save :: niters          ! number of iterations
  logical, save :: Finterp
  integer, save :: spatial_order   ! spatial order for operators
  integer, save :: interp_order
  integer, save ::  mg_interp_order
  integer, save ::  do_spec
  integer, save ::  N_Vcycles
  integer, save ::  Nrelax

  character(len=32), save :: pfasst_nml
  character(len=20), save :: fbase   !  base name for run
  character(len=64), save :: foutbase   !  base name for output file
  character(len=44) :: foutfoo          !  temp used to create foutbase
  integer, save    :: poutmod

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/ Finterp, ndim,nlevs, nnodes, nvars,nprob, nsteps, niters
  namelist /params/ spatial_order,interp_order, mg_interp_order, do_spec, N_Vcycles,Nrelax
  namelist /params/ window_type, pfasst_nml,fbase ,poutmod
  namelist /params/  abs_tol, rel_tol 
  namelist /params/  v, nu, t0, dt, Tfin,sigma, kfreq

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename
    integer :: i
    CHARACTER(len=32) :: arg
    integer :: un


    !
    ! defaults
    !

    window_type  = "block"
    Finterp = .FALSE.
    ndim     = 1
    nlevs   = 3
    nnodes  = [ 2, 3, 3 ]

    niters  = 8
    nsteps  = -1

    v       = 0.0_pfdp
    Lx      = 1._pfdp
    nu      = 0.02_pfdp
    sigma   = 0.004_pfdp
    kfreq   = 1
    t0      = 0.25_pfdp
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp

    abs_res_tol = 0.0
    rel_res_tol = 0.0

    spatial_order=2
    interp_order = 2

    do_spec = 1
    N_Vcycles = 1
    Nrelax = 1
    mg_interp_order = 2
    
    poutmod = 1
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

  end subroutine probin_init

end module probin

