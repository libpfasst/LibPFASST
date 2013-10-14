module probin
  use pf_mod_dtype

  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: two_pi = 6.2831853071795862d0

  integer, parameter :: maxlevs = 3

  integer, parameter :: PROB_AD    = 11
  integer, parameter :: PROB_HEAT  = 12
  integer, parameter :: PROB_VB    = 13
  integer, parameter :: PROB_WAVE  = 14
  integer, parameter :: PROB_KS    = 15
  integer, parameter :: PROB_SHEAR = 21

  integer, save :: problem
  integer, save :: wtype

  character(len=64), save :: problem_type, window_type

  double precision, save :: v      ! advection velocity (PROB_AD only)
  double precision, save :: Lx     ! domain size
  double precision, save :: nu     ! viscosity
  double precision, save :: t0     ! initial time for exact solution (PROB_AD only)
  double precision, save :: sigma  ! initial condition parameter
  double precision, save :: delta  ! initial condition parameter
  double precision, save :: rho    ! initial condition parameter
  double precision, save :: dt     ! time step

  double precision, save :: abs_tol ! absolute residual tolerance
  double precision, save :: rel_tol ! relative residual tolerance

  integer, save :: dim             ! number of dimensions
  integer, save :: nlevs           ! number of pfasst levels
  integer, save :: nnodes(maxlevs) ! number of nodes
  integer, save :: nvars(maxlevs)  ! number of grid points
  integer, save :: nsteps          ! number of time steps
  integer, save :: niters          ! number of iterations


  character(len=64), save :: output ! directory name for output
  

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename

    integer :: un

    namelist /prbin/ &
         problem_type, window_type, output, abs_tol, rel_tol, &
         v, nu, t0, dt, sigma, rho, delta, &
         nlevs, nnodes, nvars, nsteps, niters


    !
    ! defaults
    !

    problem_type = "ad"
    window_type  = "block"
    output       = ""

    nlevs   = 2
    nnodes  = [ 2, 3, 5 ]
    nvars   = [ 32, 64, 128 ]
    niters  = 8
    nsteps  = -1

    v       = 1.d0
    Lx      = 1.d0
    nu      = 0.02d0
    sigma   = 0.004d0
    t0      = 0.25d0
    dt      = 0.01d0
    delta   = 1.0d0
    rho     = 1.0d0

    abs_tol = 0.d0
    rel_tol = 0.d0


    !
    ! read
    !

    un = 66
    open(unit=un, file=filename, status='old', action='read')
    read(unit=un, nml=prbin)
    close(unit=un)


    !
    ! init
    !

    select case (problem_type)
    case ("ad")
       problem = PROB_AD
       dim     = 1
    case ("heat")
       problem = PROB_HEAT
       dim     = 1
    case ("burgers")
       problem = PROB_VB
       dim     = 1
    case ("wave")
       problem = PROB_WAVE
       dim     = 2
    case ("ks")
       problem = PROB_KS
       nu      = -1.d0
       Lx      = 100.d0
       dim     = 1
    case ("shear")
       problem = PROB_SHEAR
       dim     = 2
    end select

    select case (window_type)
    case ("block")
       wtype = PF_WINDOW_BLOCK
    case ("ring")
       wtype = PF_WINDOW_RING
    end select

  end subroutine probin_init

end module probin

