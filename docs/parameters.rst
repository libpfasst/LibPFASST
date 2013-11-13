
Parameters and variables
========================

The libpfasst library has many parameters which control the 
behavior of the PFASST algorithm and can be changed by the 
user.  This section lists all the parameters and describes
their function, location, and default values. Most of the 
parameters assume a default value that can be changed by
specifying the value either in an input file or on the 
command line.  Some of the parameters must be changed from
their default value or PFASST will not execute.  

Following these lists is an explanation of how to set 
parameters through input files or the command line, and 
how to choose certain parameters to acheive particular 
variants of PFASST.

Types of parameters
-------------------

*  libpfasst static parameters:  are hard coded and cannot be changed at run time
*  pfasst_t mandatory parameters: must be reassigned at run time 
    the use of default values will result in program termination.
*  pfasst_t optional parameters: can  be reassigned at run time,  
   the use of default values will result in default execution.
*  level_t mandatory  parameters: must be reassigned at run time,  
   the use of default values will result in default execution.
*  level_t optional parameters: can  be reassigned at run time, 
   the use of default values will result in default execution.
*  local mandatory parameters:  must be passed in a call to pf_run_pfasst 
*  local optional parameters:   specified by the user application and unrelated to the workings of libpfasst


libfpasst static parameters:
---------------------------

The parameters at the top of the file src/pf_dtype.f90 are all set at compile time and can't be changed at runtime.
The only parameter here of interest to the user is 

  integer, parameter :: pfdp = c_double

which controls the precision of all floating point numbers (or at least all those using pfdp in the declaration)

pfasst_t mandatory parameters:
---------------------------

The parameters defined in type  pf_pfasst_t the file src/pf_dtype.f90 are all given a default value.  Currently
on the variable nlevels is given a problematic default.  Hence setting this variable on the command line or in
an initialization file is mandatory

pfasst_t optional parameters:
---------------------------

The remaining variables in the specification  of  pf_pfasst_t  are given default values as below
  type :: pf_pfasst_t

     integer :: niters  = 5             ! number of iterations
     integer :: rank    = -1            ! rank of current processor
     integer :: qtype   = SDC_GAUSS_LOBATTO
     integer :: ctype   = SDC_CYCLE_V

     real(pfdp) :: abs_res_tol = 0.d0
     real(pfdp) :: rel_res_tol = 0.d0

     integer :: window = PF_WINDOW_BLOCK
 
     logical :: Pipeline_G =  .false.
     logical :: PFASST_pred = .false.

     ! timing
     logical    :: echo_timings  = .false.

These value can be changed as desired.

level_t mandatory parameters:
---------------------------

In the specification of pf_level_t, the first two variables nvar and nnodes
are given default values that will cause the program to abort.  These variables
are typically set per level in the main.f90 routine.

  type :: pf_level_t
     integer     :: nvars = -1          ! number of variables (dofs)
     integer     :: nnodes = -1         ! number of sdc nodes

level_t optional parameters:
---------------------------

In the specification of pf_level_t, the first  variables nsweeps and Finterp
are  default values.  These can be changed per level as the levels are created
in main.f90

     integer     :: nsweeps = 1         ! number of sdc sweeps to perform
     logical     :: Finterp = .false.   ! interpolate functions instead of solutions



local mandatory parameters: 
---------------------------

In the call to run pfasst

    pf_pfasst_run(pf, q0, dt, tend, nsteps, qend)

The variables q0, dt, and tend  must be included.  The variable
nsteps is optional, if it is not included, then nsteps is set to  

       pf%state%nsteps = ceiling(1.0*tend/dt)

qend is also optional and returns the final solution.



File input for user variables
-----------------------------

The usual default input file is "probin.nml" wherein the namelist 
PARAMS (defined  locally in probin.f90) can be specified.  Alternatively,
a different input file can be specified on the command line by adding
the file name directly after the executable.  The alternative input
file must be specified first before any command line parameter specifications 
(see next section).

File input for pfasst  variables
--------------------------------

The pfasst parameters are specified in a namelist PF_PARAMS defined
in routine pf_read_opts in src/pf_options.f90.  This routine is called
from pf_pfasst_create in pf_pfasst.f90 (which is typically 
called from main.f90).  If no file is specified in the call to pf_pfasst_create,
then no file is read.  Typically the main.f90 routine specifies this
input file (the default being probin.nml), and this file can be changed 
by specifying   the value of

  pfasst_nml = 'probin.nml'

either in the local input file or the command line. 

Command line input
------------------

All the variables in the namelist PARAMS (defined  locally in probin.f90) 
and PF_PARAMS can be modified by simply specifying their value on the command line.  There is 
only one caveat to this in that any parameters must be specified after the
(optional) input file specification.  For example

mpirun -n 20 main.exe  myinput.nml niters=10

would set the  input file to "myinput.nml" and then over-ride any
specified value of niters with the value 10. Command line options
over-ride input files.





Variables for the predictor
--------------------------

The two variables Pipeline_G and PFASST_pred  determine how the
predictor works.  The different combinations of these variables
and the parameter Nsweeps on the coarsest level great some subtle
differences in how the predictor performs.

Some cases:
1. If PFASST_pred is false and Pipeline_G is false, then 
the predictor is a serial application of SDC with Nsweeps.
This can be done without communication wherein every processor
mimics the behavior of the processors previous to it in time.

2. If PFASST_pred is false and Pipeline_G is true and Nsweeps
is one, then the predictor is a serial application of SDC with 1
sweep.  As above, there is no communication necessary.

3. If PFASST_pred is false and Pipeline_G is true and Nsweeps
is greater than one,  then the predictor is a version of pipelined
SDC. There is no communication necessary until the second sweep on
the each processor is done.  After that, each processor must recieve
a new initial value.

4. If PFASST_pred is true, and Nsweeps equals one, then it doesn't
matter what Pipeline_G is.  No communication is necessary, and 
we simply reuse the function values from the previous iteration
in each SDC sweep.  Some care must be taken here as to how to 
interpret the variable t0 especially in light of time dependent
boundary conditions.  Currently t0 does not change in these
iterations, hence one should use caution using PFASST_pred = true with
time dependent boundary conditions.

5. If PFASST_pred is true, and Nsweeps is greater than  one and 
Pipeline_G is true, then the predictor will act like the normal
PFASST_pred with Nsweeps equals one, but more iterations will be
taken.  This choice is a bit strange.  No communication is needed
until each processor is doing the P+1st iteration, then new
initial data must be used and in all cases, previous f values are
used in the SDCsweeps.  The caveat about t0 is still valid.    

6. Finally, if PFASST_pred is true, and Nsweeps is greater than  one and 
Pipeline_G is false, then the predictor will act like the normal
PFASST_pred with Nsweeps equals one, but additional iterations
are taken before the initial conditions at each processor are reset.
This can be done without communication.
The caveat about t0 is still valid.    

How is this implemented?  There are two pieces to the initialization.
The first consists of the process of giving every processor an initial
value which is consistent with t0 at that processor.  This can be
done without communication in all cases.


