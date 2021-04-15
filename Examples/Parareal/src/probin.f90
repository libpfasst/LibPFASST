!
! This file is part of LIBPFASST.
!
!>  Module for reading local parameters for the problem
module probin
  use pfasst
  use pf_mod_dim, only: echo_dim

  character(len=64), save :: problem_type
  integer,  save :: Ndim   ! Number of dimesions
  real(pfdp),  save,allocatable :: dom_size(:)    ! Domain size
  real(pfdp),  save :: Lx,Ly,Lz    ! Components of domain size
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time
  integer, save, allocatable ::  grid_size(:)  !  Will hold the size of the grid
  real(pfdp),  save,allocatable :: kfreq(:)  ! wave numbers
  integer, save :: nx(PF_MAXLEVS)     ! number of grid points per level
  integer, save :: ic_type         ! which initial condition
  integer, save :: eq_type         ! which equation to solve
  integer, save :: nsteps          ! number of time steps
  integer, save :: nsteps_rk(PF_MAXLEVS)       ! number of time steps for rk
  integer, save :: rk_order(PF_MAXLEVS)         ! number of time steps for rk
  integer, save :: splitting       ! type of imex splitting
  !  parameters for advection diffusion
  real(pfdp), save :: lam1,lam2 ! coefficients for Dahlquist
  real(pfdp), save :: a,b,c   ! advection velocities
  real(pfdp), save :: nu      ! viscosity
  real(pfdp), save :: t00     ! initial time for exact ad solution
  real(pfdp), save :: sigma   ! initial condition parameter for ad solution
  real(pfdp),  save :: kfreqx  ! initial condition parameter for ad solution
  real(pfdp),  save :: kfreqy  ! initial condition parameter for ad solution
  real(pfdp),  save :: kfreqz  ! initial condition parameter for ad solution
  !  parameters for kdv
  real(pfdp), save :: beta    ! scaling factor
  !  parameters for burgers term
  real(pfdp), save :: gamma    ! scaling factor


  character(len=32), save :: pfasst_nml

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat 
  namelist /params/  nx,ic_type, eq_type, nsteps,nsteps_rk,rk_order, dt, Tfin
  namelist /params/  pfasst_nml, lam1,lam2,a,b,c, nu, t00, sigma, beta, gamma, splitting
  namelist /params/  kfreqx,kfreqy,kfreqz,Lx,Ly,Lz

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

    Ndim=echo_dim()  !  in dim module
    !  Allocate grid_size to have correct dimension
    allocate(grid_size(Ndim))
    allocate(dom_size(Ndim))
    allocate(kfreq(Ndim))
    !> set defaults
    nsteps  = -1
    nsteps_rk  = -1
    rk_order  = 2

    lam1    = -1.0_pfdp
    lam2    = 0.5_pfdp
    a       = 0.0_pfdp
    b       = 1.0_pfdp
    c       = 1.0_pfdp
    nu      = 0.01_pfdp
    kfreqx   = 1.0_pfdp
    kfreqy   = 1.0_pfdp
    kfreqz   = 1.0_pfdp
    beta    = 3.0000_pfdp
    gamma   = 1.0000_pfdp
    Lx      = two_pi
    Ly      = two_pi
    Lz      = two_pi
    t00      = 0.08_pfdp
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp
    eq_type = 1  !  1: AD 2: Burgers 3: NLS 4: Kdv, 5: KS 7:ZD
    ic_type = 1  !  0: Gaussian, 1: Sin wave
    splitting = 2  !  Default is fully exponential (exact)
    pfasst_nml=probin_fname

    !>  Read in stuff from input file
    un = 9
    ! write(*,*) 'opening file ',TRIM(probin_fname), '  for input'
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

    !  Load the domain size and kfreq
    dom_size(1)=Lx
    kfreq(1)=kfreqx
    if (Ndim .gt. 1) then
       dom_size(2)=Ly
       kfreq(2)=kfreqy
    end if
    if (Ndim .gt. 2) then
       dom_size(3)=Lz
       kfreq(3)=kfreqz
    end if
       
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
    write(un,*) 'Ndim:   ', Ndim, '! Number of dimensions'
    write(un,*) 'nsteps: ', nsteps, '! Number of steps'
    write(un,*) 'nsteps_rk: ', nsteps_rk(1:pf%nlevels), '! Number of rk substeps'
    write(un,*) 'rk_order: ', rk_order(1:pf%nlevels), '! Order of rk substeps'        
    write(un,*) 'Dt:     ', Dt, '! Time step size'
    write(un,*) 'Tfin:   ', Tfin,   '! Final time of run'
    write(un,*) 'nx:     ',  nx(1:pf%nlevels), '! grid size per level'
    write(un,*) 'Domain:     ', Lx,Ly,Lz, '! domain size'
    write(un,*) 'splitting:', splitting, '! 1,2, or 3'

    select case (eq_type)
    case (0)  
       write(un,*) 'Solving the linear Dahlquist equation'
       write(un,*) 'lam1,lam2', lam1,lam2, '! Dahlquist constants'
    case (1)  
       write(un,*) 'Solving the linear advection diffusion equation'
       write(un,*) 'adevct coef: ',  a,b,c, '! advection velocities'
       write(un,*) 'nu:     ', nu, '! diffusion constant'
       select case (ic_type)
       case (0)  
          write(un,*) 'Periodic Gaussian initial conditions with t00=',t00
       case (1)  
          write(un,*) 'Sine initial conditions with kfreq=',kfreq
       case DEFAULT
          call pf_stop(__FILE__,__LINE__,'Bad case  for ic_type ',ic_type)
       end select
    case (2)  
       write(un,*) 'Solving Burgers equation'
    case (3)  
       write(un,*) 'Solving the nonlinear Schroedinger equation'
    case (4)  
       write(un,*) 'Solving the KDV equation'
    case (5)  
       write(un,*) 'Solving the Kuramoto-Shivashinsky equation'
    case (6)  
       write(un,*) 'Solving the Reaction-diffusion equation'
    case (7)  
       write(un,*) 'Solving the Zero dispersion optics'
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case  for eq_type ',eq_type)
    end select

    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  

end module probin
