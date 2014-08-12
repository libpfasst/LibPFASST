module probin
  use pf_mod_dtype

  integer, parameter :: maxlevs = 3

  integer, save           :: nlevs, nsteps, niters, nskip
  integer, save           :: nx(maxlevs), nvars(maxlevs), nnodes(maxlevs)
  double precision, save  :: nu, dt, sigma
  character(len=64), save :: input, output

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename

    integer :: un

    namelist /prbin/ input, output, nu, dt, nlevs, nx, nsteps, niters, nnodes, sigma, nskip

    !
    ! defaults
    !

    nlevs   = 2
    nnodes  = [ 2, 3, 5 ]
    nx      = [ 16, 32, 64 ]
    niters  = 8
    nsteps  = -1
    nskip   = 1

    nu      = 7.6d-4
    dt      = 0.0001d0
    sigma   = 99.0d0

    !
    ! read
    !

    un = 66
    open(unit=un, file=filename, status='old', action='read')
    read(unit=un, nml=prbin)
    close(unit=un)

    nvars  = 2 * 3 * nx**3

    if (nsteps < 0) stop "nsteps must be set"

  end subroutine probin_init

end module probin

