!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem
module probin
  use pf_mod_dtype


  character(len=64), save :: problem_type

  real(pfdp), save :: lam1
  real(pfdp), save :: lam2
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time

  integer, save :: nsteps          ! number of time steps
  character(len=32), save :: pfasst_nml  ! file for reading pfasst parameters
  
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/  lam1,lam2, dt, Tfin, nsteps,pfasst_nml


contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename
    integer :: i
    character(len=32) :: arg
    integer :: un

    !> set defaults
    nsteps  = -1

    lam1       = 1.0_pfdp
    lam2       = -2.0_pfdp

    dt      = 0.01_pfdp
    Tfin    = 1.0_pfdp
    pfasst_nml=filename
    
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

  subroutine print_loc_options(pf, un_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: un_opt
    integer :: un = 6

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt

    !  Print out the local parameters
    write(un,*) '=================================================='
    write(un,*) ' '
    write(un,*) 'Local Variables'
    write(un,*) '----------------'
    write(un,*) 'nsteps: ', nsteps, '! Number of steps'
    write(un,*) 'Dt:     ', Dt, '! Time step size'
    write(un,*) 'Tfin:   ', Tfin,   '! Final time of run'
    write(un,*) 'lam1:   ', lam1, '! explicit constant'
    write(un,*) 'lam1:   ', lam2, '! implicit constant'    


    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  
  
  
end module probin
