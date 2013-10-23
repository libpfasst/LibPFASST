  type :: pf_pfasst_t
     integer :: nlevels = -1            ! number of pfasst levels
     integer :: niters  = 5             ! number of iterations
     integer :: rank    = -1            ! rank of current processor
     integer :: qtype   = SDC_GAUSS_LOBATTO
     integer :: ctype   = SDC_CYCLE_V

     real(pfdp) :: abs_res_tol = 0.d0
     real(pfdp) :: rel_res_tol = 0.d0

     integer :: window = PF_WINDOW_BLOCK
     
  type :: pf_level_t
     integer     :: nvars = -1          ! number of variables (dofs)
     integer     :: nnodes = -1         ! number of sdc nodes
     integer     :: nsweeps = 1         ! number of sdc sweeps to perform
     integer     :: level = -1          ! level number (1 is the coarsest)
     logical     :: Finterp = .false.   ! interpolate functions instead of solutions



  type, bind(c) :: pf_state_t
     integer(c_int) :: block
     integer(c_int) :: cycle
     integer(c_int) :: hook
     integer(c_int) :: first        ! rank of first processor in time block
     integer(c_int) :: iter
     integer(c_int) :: last         ! rank of last processor in time block
     integer(c_int) :: level
     integer(c_int) :: nmoved       ! how many processors behind me have moved
     integer(c_int) :: nsteps
     integer(c_int) :: proc
     integer(c_int) :: pstatus      ! previous rank's status
     integer(c_int) :: step
     integer(c_int) :: status       ! status (iterating, converged etc)

     real(c_double) :: dt
     real(c_double) :: res
     real(c_double) :: t0
  end type pf_state_t


  namelist /params/ Finterp       !  Interpolate function values too?
  namelist /params/ nnodes        !  Number of nodes in each level
  namelist /params/ nfake         !  number of fake processors
  namelist /params/ niters        !  number of pfasst iterations
  namelist /params/ nlevs         !  number of PFASST levels
  namelist /params/ nprob         !  define which problem to run
  namelist /params/ nsteps        !  Number of time steps to take
  namelist /params/ qtype         !  type of quadrature nodes
  namelist /params/ poutmod       !  controls how often output comes

  namelist /params/ dt            !  time step
  namelist /params/ Tfin          !  End time of run
  namelist /params/ Htol          !  Tolerance for stopping


  namelist /params/ fbase          !  Base name for output
  namelist /params/ fnml          !  Base name for output
