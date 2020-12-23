!
! This file is part of LIBPFASST.
!
module probin
  use pfasst
  use pf_mod_mpi

  !  The namlist for local variables
  integer, save :: num_grid_points, nspace, ntime, space_dim, max_space_v_cycles
  real(pfdp), save :: init_cond
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time
  integer, save :: nsteps    ! number of time steps
  character(len=128), save :: pfasst_nml  ! file for reading pfasst parameters

  integer, save :: nsteps_rk(PF_MAXLEVS)
  integer, save :: rk_order
  logical, save :: use_mgrit
  integer, save :: mgrit_n_coarse
  integer, save :: mgrit_refine_factor

  integer, save :: imex_stat
  integer, save :: ark_stat

  namelist /params/  space_dim, num_grid_points, init_cond, nspace, ntime, dt, Tfin, nsteps, pfasst_nml, max_space_v_cycles
  namelist /params/  nsteps_rk, use_mgrit, mgrit_n_coarse, mgrit_refine_factor, rk_order, imex_stat, ark_stat

contains
  
  subroutine probin_init(pf_fname)
    character(len=*), intent(inout) :: pf_fname

    !  Some local variables for reading
    character(len=128) :: arg  !  command line argument
    character(len=128)    :: probin_fname   !<  file name for input parameters
    character(len=256) :: istring           ! stores command line argument
    character(len=1024) :: message          ! use for I/O error messages
    integer :: ios,iostat
    integer :: i   !  loop variable
    integer :: un  !  file read unit
    integer :: nproc, rank, error

    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

    !> Set the name of the input file
    probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)

    !> set defaults
    nsteps  = -1

    nspace = 1
    ntime = nproc
    num_grid_points = 3
    init_cond = 50.0
    space_dim = 2
    max_space_v_cycles = 100

    dt      = 0.01_pfdp
    Tfin    = 1.0_pfdp
    pfasst_nml=probin_fname

    mgrit_n_coarse = 10
    mgrit_refine_factor = 2
    nsteps_rk = -1
    use_mgrit = .true.
    rk_order = 1
    nsteps_rk = -1
    
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
          READ(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
       end if
       i = i+1
    end do

    !  Reset dt if Tfin is set
    if (Tfin .gt. 0.0) dt = Tfin/dble(nsteps)

    !  Return the name of the file from which to read PFASST parameters
    pf_fname=pfasst_nml
  end subroutine probin_init

  !>  Subroutine to output run parameters 
  subroutine print_loc_options(pf, un_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: un_opt
    integer :: un = 6
    
    integer :: nproc, rank, error

    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt

    !>  Output the PFASST options with the LibPFASST routine
    call pf_print_options(pf,un_opt=un)


    if (rank == 0) then
       !  Print out the local parameters
       write(un,*) '=================================================='
       write(un,*) ' '
       write(un,*) 'Local Variables'
       write(un,*) '----------------'
       write(un,*) 'nsteps: ', nsteps, '! Number of steps'
       write(un,*) 'Dt:     ', Dt, '! Time step size'
       write(un,*) 'Tfin:   ', Tfin,   '! Final time of run'
       write(un,*) 'num spatial grid points (on each side of square domain) per processor:   ', num_grid_points
       write(un,*) 'num spatial procs per temporal proc:   ',   nspace
       write(un,*) 'num temporal procs:   ',   ntime  
       write(un,*) 'Constant initial condition:   ', init_cond
       write(un,*) 'Number of spacial dimensions:   ', space_dim
       write(un,*) 'Number of Hypre V-cycles:   ', max_space_v_cycles
   
   
       write(un,*) 'PFASST parameters read from input file ', pfasst_nml
       write(un,*) '=================================================='
    end if
  end subroutine print_loc_options
  
  
  
end module probin
