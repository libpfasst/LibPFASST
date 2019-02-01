! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin

  use pf_mod_dtype
  implicit none

  integer, parameter :: MAXLEVS = 4
  integer, save    :: Finterp         !  Decides if interpolation of function values is done
  integer, save    :: ndim            !  Number of dimensions
  integer, save    :: nfake           !  Number of fake processors
  integer, save    :: nnodes(MAXLEVS) !  number of nodes on each level
  integer, save    :: nprob           !  define which problem to run

  integer, save    :: nsteps          !  Number of time steps to take
  integer, save    :: N(MAXLEVS)      !  Size of problem on each level
  integer, save    :: magnus_order(MAXLEVS)    !  The order of the Magnus expansion to use
  real(pfdp), save :: dt              !  Time step
  real(pfdp), save :: Tfin            !  Final time of run
  real(pfdp), save :: alpha           !  energy of state A for 2x2 problem
  real(pfdp), save :: beta            !  energy of state B for 2x2 problem
  real(pfdp), save :: Znuc            !  Single integral cst
  real(pfdp), save :: E0              !  Double integral cst
  real(pfdp), save :: Xmax            !  Max distance from origin
  real(pfdp), save :: vab             !  coupling of states A, B for 2x2 problem
  real(pfdp), save :: exptol(MAXLEVS)
  integer, save    :: nsweeps(MAXLEVS) !  Sweeps at each levels
  integer, save    :: nsweeps_pred(MAXLEVS)   !  Sweeps during predictor
  ! character(len=64), save :: basis(MAXLEVS)
  character(len=64), save :: basis
  character(len=1024), save :: molecule
  character(len=1024), save :: exact_dir
  integer, save :: nparticles
  logical, save :: toda_periodic, save_solutions, use_sdc

  !  Output
  character(len=64), save :: fbase
  integer, save    :: poutmod

  namelist /params/ Finterp, ndim, nfake, nnodes, nprob, nsteps, N, dt, Tfin
  namelist /params/ nsweeps,nsweeps_pred, use_sdc
  namelist /params/ fbase, poutmod, exptol, save_solutions, nparticles, toda_periodic
  namelist /params/ alpha, beta, vab, magnus_order, basis, molecule, exact_dir,E0,Znuc

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename

    integer,parameter :: un=9
    integer :: ird,ios
    CHARACTER(len=32) :: arg
    CHARACTER(LEN=64) :: istring   ! For pulling off command line options
    CHARACTER(LEN=64) :: message   ! use for I/O error messages

    ! initialize the runtime parameters

    ! Defaults
    Finterp = 0
    ndim = 2
    nfake = 0
    nprob = 3
    nsteps = -1
    N = -1
    poutmod = 1
    fbase = 'output'

    dt = 0.01_pfdp
    Tfin = 1.0_pfdp

    alpha = -1.0_pfdp
    beta  = -0.5_pfdp
    vab   = 0.25_pfdp
    E0=0.2
    Znuc=2.0
    Xmax=1.0        
    basis = "sto3g"
    molecule = "H 0 0 0; H 0 0 1.414"
    exact_dir = '' ! NO DEFAULT INITIAL DMAT
    magnus_order(:) = 3
    exptol(:) = 1.0d-20
    save_solutions = .true.
    toda_periodic = .false.
    use_sdc = .false.
    nparticles = 3

    nnodes = 3         !  Set them all to three

    open(unit=un, file = filename, status = 'old', action = 'read')
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
       end if
       ird = ird+1
    END DO

    if( Tfin .gt. 0.0_pfdp )then
       dt=Tfin/nsteps
    else
       Tfin = dt*nsteps
    endif

    do ird = 1, MAXLEVS
      if (nnodes(ird) < 2) then
        magnus_order(ird) = 1
      else
      endif
    enddo
  end subroutine probin_init

end module probin
