!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem
module probin
  use pf_mod_dtype


  character(len=64), save :: problem_type

  double precision, save :: v      ! advection velocity
  double precision, save :: nu     ! viscosity
  double precision, save :: t00     ! initial time for exact solution
  double precision, save :: sigma  ! initial condition parameter
  integer,          save :: kfreq  ! initial condition parameter
  double precision, save :: dt     ! time step
  double precision, save :: Tfin   ! Final time

  integer, save :: nx(PF_MAXLEVS)     ! number of grid points
  integer, save :: nprob           ! which problem
  integer, save :: nsteps          ! number of time steps
  integer, save :: nsteps_rk       ! number of time steps for rk
  integer, save :: imex_stat       ! type of imex splitting

  character(len=32), save :: pfasst_nml

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/  nx,nprob, nsteps,nsteps_rk, dt, Tfin
  namelist /params/  pfasst_nml, v, nu, t00, sigma, kfreq,imex_stat

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename
    integer :: i
    character(len=32) :: arg
    integer :: un

    !> set defaults
    nsteps  = -1
    nsteps_rk  = -1

    v       = 1.0_pfdp
    nu      = 0.01_pfdp
    kfreq   = 1.0_pfdp
    t00      = 0.08_pfdp
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp
    nprob = 0  !  0: Gaussian, 1: Sin wave
    imex_stat=2    !  Default is full IMEX
    

    !>  Read in stuff from input file
    un = 9
    write(*,*) 'opening file ',TRIM(filename), '  for input'
    open(unit=un, file = filename, status = 'old', action = 'read')
    read(unit=un, nml = params)
    close(unit=un)
          
    !>  Read the command line
    i = 0
    do
       call get_command_argument(i, arg)
       if (LEN_TRIM(arg) == 0) EXIT
       if (i > 0) then
          istring="&PARAMS "//TRIM(arg)//" /"    
          READ(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
       end if
       i = i+1
    end do

    !  Reset dt if Tfin is set
    if (Tfin .gt. 0.0) dt = Tfin/dble(nsteps)
  end subroutine probin_init

end module probin
