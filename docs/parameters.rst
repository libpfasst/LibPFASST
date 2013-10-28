
Parameters and variables
========================

The libpfasst library has many parameters which control the 
behavior of the PFASST algorithm and can be changed by the 
user.  This section lists all the parameters and describes
their function, location, default values, 
and how they are set.  Below is also an explanation of how to set 
certain parameters to acheive particular variants of PFASST.

Types of parameters
-------------------

*  libpfasst static parameters:  are hard coded and cannot be changed at run time
*  libfpasst mandatory paramters: must be reassigned at run time by either specification in and input file or the command line.  The use of default values will result in program termination.
*  libfpasst optional paramters: can  be reassigned at run time by either specification in and input file or the command line.  The use of default values will result in default execution.
* local mandatory parameters:  must be passed in a call to pf_run_pfasst
* local optional parameters:   specified by the user application and unrelated to the workings of libpfasst



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



File input
----------

The default input file is "probin.nml" wherein the namelist 
PARAMS (defined  locally in probin.f90) can be specified.  Alternatively,
a different input file can be specified on the command line by adding
the file name directly after the executable.  The alternative input
file must be specified first before any command line parameter specifications 
(see next section).

Command line input
------------------

All the variables in the namelist PARAMS (defined  locally in probin.f90) can
be modified by simply specifying their value on the command line.  There is 
only one caveat to this in that any parameters must be specified after the
(optional) input file specification.  For example

mpirun -n 20 main.exe  myinput.nml niters=10

would set the  input file to "myinput.nml" and then over-ride any
specified value of niters with the value 10. 
