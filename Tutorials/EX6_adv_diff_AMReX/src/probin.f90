!
! This file is part of LIBPFASST.
!
!>  Module for reading local parameters for the problem
module probin
  use pfasst
  use pf_mod_dim, only: echo_dim
  use pf_mod_mpi
  
  character(len=64), save :: problem_type
  integer,  save :: Ndim   ! Number of dimesions
  integer,  save :: nspace   ! Number of process in space
  integer,  save :: ntime   ! Number of process in time
  integer,  save :: ndiag   ! Number of process for diag sweeper
  real(pfdp),  save,allocatable :: dom_size(:)    ! Domain size
  real(pfdp),  save :: Lx,Ly,Lz    ! Components of domain size
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time
  integer, save, allocatable ::  grid_size(:)  !  Will hold the size of the grid
  
  real(pfdp),  save,allocatable :: kfreq(:)     ! wave numbers
  integer, save :: nx(PF_MAXLEVS)               ! number of grid points per level x-Dim
  integer, save :: ny(PF_MAXLEVS)               ! number of grid points per level y-Dim
  integer, save :: nz(PF_MAXLEVS)               ! number of grid points per level z-Dim
  integer, save :: ic_type                      ! which initial condition
  integer, save :: eq_type                      ! which equation to solve
  integer, save :: nsteps                       ! number of time steps
  integer, save :: nsteps_rk(PF_MAXLEVS)        ! number of time steps for rk
  integer, save :: rk_order(PF_MAXLEVS)         ! number of time steps for rk
  integer, save :: splitting                    ! type of imex splitting
  !  parameters for advection diffusion
  real(pfdp), save :: lam1,lam2                 ! coefficients for Dahlquist
  real(pfdp), save :: vx,vy,vz                     ! advection velocities
  real(pfdp), save :: nu                        ! viscosity
  real(pfdp), save :: t00                       ! initial time for exact ad solution
  real(pfdp), save :: sigma                     ! initial condition parameter for ad solution
  real(pfdp),  save :: kfreqx                   ! initial condition parameter for ad solution
  real(pfdp),  save :: kfreqy                   ! initial condition parameter for ad solution
  real(pfdp),  save :: kfreqz                   ! initial condition parameter for ad solution
  ! parameters for AMReX
  integer, save :: max_grid_size    ! maximum number of grid points per dimenson per box (AMReX defaults: 2D -> 128, 3D -> 32)
  integer, save :: linop_maxorder   ! maximum order of linear operator
  integer, save, allocatable :: linop_bc_lo(:)      ! boundary conditions for MLMG Laplacian
  integer, save, allocatable :: linop_bc_hi(:)      ! boundary conditions for MLMG Laplacian
  integer, save, allocatable :: geom_bc_lo(:)       ! geometrical boundary conditions - needed for interpolation of results
  integer, save, allocatable :: geom_bc_hi(:)       ! geometrical boundary conditions - needed for interpolation of results
  integer, save :: mlmg_max_iter    ! maximum number of MLMG iterations
  real(pfdp), save :: mlmg_tol_rel  ! relative tolerance for MLMG solver
  real(pfdp), save :: mlmg_tol_abs  ! absolute tolerance for MLMG solver
  !  parameters for kdv
  real(pfdp), save :: beta    ! scaling factor
  !  parameters for burgers term
  real(pfdp), save :: gamma    ! scaling factor


  character(len=32), save :: pfasst_nml

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat 
  namelist /params/  nx,ny,nz,ic_type, eq_type, nsteps,nsteps_rk,rk_order, dt, Tfin
  namelist /params/  pfasst_nml, lam1,lam2,vx,vy,vz, nu, t00, sigma, beta, gamma, splitting
  namelist /params/  kfreqx,kfreqy,kfreqz,Lx,Ly,Lz, max_grid_size
  namelist /params/  linop_maxorder, linop_bc_lo, linop_bc_hi, geom_bc_lo, geom_bc_hi
  namelist /params/  mlmg_max_iter, mlmg_tol_rel, mlmg_tol_abs, nspace, ntime, ndiag

contains

  subroutine probin_init(pf_fname)
    character(len=*), intent(inout) :: pf_fname
    integer :: i   !  loop variable
    integer :: un  !  file read unit
    character(len=32) :: arg  !  command line argument
    character(128)    :: probin_fname   !<  file name for input parameters
    integer :: gSize, gRank, err

    !> Set the name of the input file
    probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)

    Ndim=echo_dim()  !  in dim module
    !  Allocate grid_size to have correct dimension
    allocate(grid_size(Ndim+1))  !  We have the ncomp first, then the grid size (ncomp,nx,ny,nz)
    allocate(dom_size(Ndim))
    allocate(kfreq(Ndim))
    allocate(linop_bc_lo(Ndim))
    allocate(linop_bc_hi(Ndim))
    allocate(geom_bc_lo(Ndim))
    allocate(geom_bc_hi(Ndim))
    
    !> set defaults
    nsteps  = -1
    nsteps_rk  = -1
    rk_order  = 2

    lam1          = 0.0_pfdp  !-1.0_pfdp
    lam2          = 0.0_pfdp  !0.5_pfdp
    vx            = 1.0_pfdp
    vy            = 1.0_pfdp
    vz            = 0.0_pfdp
    nu            = 0.01_pfdp
    kfreqx        = 1.0_pfdp
    kfreqy        = 0.0_pfdp
    kfreqz        = 0.0_pfdp
    beta          = 3.0000_pfdp
    gamma         = 1.0000_pfdp
    Lx            = 1.0_pfdp
    Ly            = 1.0_pfdp
    Lz            = 1.0_pfdp
    max_grid_size = 32         ! AMReX default is 128 for 2D and 32 for 3D
    t00           = 0.08_pfdp
    dt            = 0.01_pfdp
    Tfin          = 0.0_pfdp
    eq_type       = 1  !  1: AD 2: Burgers 3: NLS 4: Kdv, 5: KS
    ic_type       = 1  !  0: Gaussian, 1: Sin wave
    splitting     = 1  !  Default is IMEX
    pfasst_nml    = probin_fname

    !> MPI-communicator stuff
    ! default to time only paralleliation
    call mpi_comm_size(MPI_COMM_WORLD, gSize, err)
    call mpi_comm_rank(MPI_COMM_WORLD, gRank, err)
    ntime = 1
    ndiag = 1
    nspace = gSize      

    !>  Read in stuff from input file
    un = 9
    if (gRank .eq. 0) write(*,*) 'opening file ',TRIM(probin_fname), '  for input'
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
    write(un,*) 'nSpace: ', nspace, '! Number of proccesses in space'
    write(un,*) 'nTime:  ', ntime, '! Number of proccesses in time'
    write(un,*) 'nDiag:  ', ndiag, '! Number of proccesses for diagonal sweeper'
    write(un,*) 'nsteps: ', nsteps, '! Number of steps'
    write(un,*) 'nsteps_rk: ', nsteps_rk(1:pf%nlevels), '! Number of rk substeps'
    write(un,*) 'rk_order: ', rk_order(1:pf%nlevels), '! Order of rk substeps'        
    write(un,*) 'Dt:     ', Dt, '! Time step size'
    write(un,*) 'Tfin:   ', Tfin,   '! Final time of run'
    write(un,*) 'nx:     ',  nx(1:pf%nlevels), '! grid size per level'
    write(un,*) 'ny:     ',  ny(1:pf%nlevels), '! grid size per level'
    write(un,*) 'nz:     ',  nz(1:pf%nlevels), '! grid size per level'
    write(un,*) 'max_grid_size:   ',  max_grid_size, '! max nodes per dimension per box'
    write(un,*) 'Domain-size:     ', dom_size, '! domain size'
    if (Ndim == 1) then
      write(un,*) 'Domain:     ', Lx, '! domain size'
    else if (Ndim == 2) then
      write(un,*) 'Domain:     ', Lx,Ly, '! domain size'
    else if (Ndim == 3) then
      write(un,*) 'Domain:     ', Lx,Ly,Lz, '! domain size'
    else if (Ndim > 3) then
      call pf_stop(__FILE__,__LINE__,'Dimensionality > 3 !')
    end if
    write(un,*) 'splitting:', splitting, '! 0 fully implicit, 1 treats u_xx implicitly and u_x explicitly, 2 fully explicit'

    select case (eq_type)
    case (0)  
       write(un,*) 'Solving the linear Dahlquist equation'
       write(un,*) 'lam1,lam2', lam1,lam2, '! Dahlquist constants'
    case (1)  
       write(un,*) 'Solving the linear advection diffusion equation'
       select case (Ndim)
       case (1)
         write(un,*) 'adevct coef: ',  vx, '! advection velocities'
       case (2)
         write(un,*) 'adevct coef: ',  vx,vy, '! advection velocities'
       case (3)
         write(un,*) 'adevct coef: ',  vx,vy,vz, '! advection velocities'
       end select
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
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case  for eq_type ',eq_type)
    end select

    write(un,*) 'AMReX parameters:'
    write(un,*) 'linop_maxorder: ', linop_maxorder
    write(un,*) 'linop_bc_lo: ', linop_bc_lo
    write(un,*) 'linop_bc_hi: ', linop_bc_hi
    write(un,*) 'geom_bc_lo: ', geom_bc_lo
    write(un,*) 'geom_bc_hi: ', geom_bc_hi
    write(un,*) 'mlmg_max_iter: ', mlmg_max_iter
    write(un,*) 'mlmg_tol_rel: ', mlmg_tol_rel
    write(un,*) 'mlmg_tol_abs: ', mlmg_tol_abs

    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  

end module probin
