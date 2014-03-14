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
         dim, window_type, output, abs_tol, rel_tol, &
         v, nu, t0, dt, sigma, &
         nlevs, nnodes, nvars, nsteps, niters

    !
    ! defaults
    !

    window_type  = "block"
    output       = ""

    dim     = 1
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

    select case (window_type)
    case ("block")
       wtype = PF_WINDOW_BLOCK
    case ("ring")
       wtype = PF_WINDOW_RING
    end select

  end subroutine probin_init

end module probin

