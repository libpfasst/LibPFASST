!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem
module probin
  use pf_mod_dtype


  character(len=64), save :: problem_type
  real(pfdp),  save :: Lx(2) ! Domain size
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time
  

  integer, save :: nx(PF_MAXLEVS)     ! number of grid points
  integer, save :: ic_type         ! which initial condition
  integer, save :: eq_type         ! which equation to solve
  integer, save :: nsteps          ! number of time steps
  integer, save :: nsteps_rk       ! number of time steps for rk
  integer, save :: splitting       ! type of imex splitting
  !  parameters for advection diffusion
  real(pfdp), save :: nu      ! viscosity
  real(pfdp), save :: kappa    ! thermal conductivity
  real(pfdp), save :: grav     ! gravity
  real(pfdp), save :: t00     ! initial time for exact ad solution
  real(pfdp), save :: sigma   ! initial condition parameter for ad solution
  real(pfdp),  save :: kfreq(2)  ! initial condition parameter for ad solution
  real(pfdp),  save :: v0(2)     ! constant velocity

  character(len=32), save :: pfasst_nml

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/  nx,ic_type, eq_type, nsteps,nsteps_rk, dt, Tfin
  namelist /params/  pfasst_nml,  nu, kappa, grav,t00, sigma, v0, kfreq, Lx,splitting

contains

  subroutine probin_init(pf_fname)
    character(len=*), intent(inout) :: pf_fname
    integer :: i   !  loop variable
    integer :: un  !  file read unit
    character(len=32) :: arg  !  command line argument
    character(128)    :: probin_fname   !<  file name for input parameters

    !> Set the name of the input file
    probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)

    !> set defaults
    nsteps  = -1
    nsteps_rk  = -1

    nu      = 0.01_pfdp
    kappa    = 0.02_pfdp
    grav    = 0.0_pfdp
    v0   = [1.0_pfdp,1.0_pfdp]
    kfreq   = [1.0_pfdp,1.0_pfdp]
    Lx      = [two_pi,two_pi]
    t00      = 0.08_pfdp
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp
    eq_type = 1  !  1: AD 2: Burgers 3: NLS 4: Kdv, 5: KS
    ic_type = 1  !  0: Gaussian, 1: Sin wave
    splitting = 1  !  Default is fully exponential (exact)
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
    write(un,*) 'v0:     ',  v0,'! constant advection velocity'
    write(un,*) 'nu:     ', nu, '! viscosity'
    write(un,*) 'kapp:   ', kappa, '! thermal diffusion constant'
    write(un,*) 'grav:   ', grav,  '! gravity constant'
    write(un,*) 'Lx:     ', Lx, '! domain size'
!    select case (splitting)
!!$    case (1)  
!!$       write(un,*) 'splitting:', splitting, '! Exponential u_xx, Explicit u_x'
!!$    case (2)  
!!$       write(un,*) 'splitting:', splitting, '! Exponential u_xx and u_x'
!!$    case (3)  
!!$       write(un,*) 'splitting:', splitting, '! Exponential u_x, Explicit u_xx'
!!$    case DEFAULT
!!$       print *,'Bad case for splitting in probin ', splitting
!!$       call exit(0)
!!$    end select


    select case (ic_type)
    case (0)  
       write(un,*) 'Periodic Gaussian initial conditions with t00=',t00
    case (1)  
       write(un,*) 'Sine initial conditions with kfreq=',kfreq
    case (2)  
       write(un,*) 'Taylor-Green vortex test'
    case DEFAULT
       print *,'Bad case for ic_type in probin ', ic_type
       call exit(0)
    end select
    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  

end module probin
