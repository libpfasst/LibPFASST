!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem
module probin
  use pf_mod_dtype


  character(len=64), save :: problem_type

  double precision, save :: lam1
  double precision, save :: lam2
  double precision, save :: dt     ! time step
  double precision, save :: Tfin   ! Final time

  integer, save :: nsteps          ! number of time steps

  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/  lam1,lam2, dt, Tfin, nsteps


contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename
    integer :: i
    character(len=32) :: arg
    integer :: un

    !> set defaults
    nsteps  = -1
    nsteps_rk  = -1

    lam1       = 1.0_pfdp
    lam2       = -2.0_pfdp

    dt      = 0.01_pfdp
    Tfin    = 1.0_pfdp

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
