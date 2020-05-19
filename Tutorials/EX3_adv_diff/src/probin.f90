!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem
module probin
  use pf_mod_dtype


  character(len=64), save :: problem_type

  real(pfdp), save :: v      ! advection velocity
  real(pfdp), save :: nu     ! viscosity
  real(pfdp), save :: kfreq  ! initial condition parameter
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time
  real(pfdp), save :: Lx     ! Domain size

  integer, save :: nx(PF_MAXLEVS)     ! number of grid points
  integer, save :: nsteps          ! number of time steps
  integer, save :: imex_stat       ! type of imex splitting
  integer, save :: ic_type         ! specifies the initial condition
  integer, save :: fd_ord          ! specifies the spatial order (0 for spectral, 2, or 4)
  integer, save :: interp_ord      ! specifies the spatial order of interpolation

  character(len=128), save :: pfasst_nml

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/  nx, nsteps, dt, Tfin
  namelist /params/  pfasst_nml, v, nu, kfreq,Lx,imex_stat,ic_type,fd_ord,interp_ord

contains

  subroutine probin_init(pf_fname)
    character(len=*), intent(inout) :: pf_fname
    integer :: i   !  loop variable
    integer :: un  !  file read unit
    character(len=128) :: arg  !  command line argument
    character(128)    :: probin_fname   !<  file name for input parameters
    
    !> Set the name of the input file
    probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
    
    !> set defaults
    nsteps  = -1
    v       = 1.0_pfdp
    nu      = 0.01_pfdp
    kfreq   = 1.0_pfdp
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp
    Lx      = 1.0_pfdp
    imex_stat=2    !  Default is full IMEX
    ic_type=1      !  Default is a sine wave    
    fd_ord=0       !  Default is spectral accuracy
    interp_ord=0       !  Default is spectral accuracy
    pfasst_nml=probin_fname

    !>  Read in stuff from input file
    un = 9
    write(*,*) 'opening file ',TRIM(probin_fname), '  for input'
    open(unit=un, file = probin_fname, status = 'old', action = 'read')
    read(unit=un, nml = params)
    close(unit=un)
          
    !>  Read the command line
    i = 0
    do
       call get_command_argument(i, arg)
       if (LEN_TRIM(arg) == 0) EXIT
       if (i > 0) then
          istring="&PARAMS "//TRIM(arg)//" /"    
          read(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
       end if
       i = i+1
    end do

    !  Reset dt if Tfin is set
    if (Tfin .gt. 0.0) dt = Tfin/dble(nsteps)

    !  Return the name of the file from which to read PFASST parameters
    pf_fname=pfasst_nml

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
    write(un,*) 'nx:     ',  nx(1:pf%nlevels), '! grid size per level'
    write(un,*) 'v:      ',  v, '! advection constant'
    write(un,*) 'nu:     ', nu, '! diffusion constant'
    write(un,*) 'fd_ord: ', fd_ord, '! spatial order of accuracy'
    write(un,*) 'interp_ord: ', interp_ord, '! spatial order of interpolation'
    select case (imex_stat)
    case (0)  
       write(un,*) 'imex_stat:', imex_stat, '! Fully explicit'
    case (1)  
       write(un,*) 'imex_stat:', imex_stat, '! Fully implicit'
    case (2)  
       write(un,*) 'imex_stat:', imex_stat, '! Implicit/Explicit'
    case DEFAULT
       print *,'Bad case for imex_stat in probin ', imex_stat
       call exit(0)
    write(un,*) 'ic_type:    ', ic_type, '! which initial condition'
    end select

    
    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  

end module probin
