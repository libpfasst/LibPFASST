
Parameters and variables
========================

The LibPFASST library has many parameters that control the
behavior of the PFASST algorithm and can be changed by the
user.  This section lists all the parameters and describes
their function, location, and default values. Most of the
parameters assume a default value that can be changed by
specifying the value either in an input file or on the
command line.  Some of the parameters must be changed from
their default value or PFASST will not execute.

Following these lists is an explanation of how to set
parameters through input files or the command line and
how to choose certain parameters to achieve particular
variants of PFASST.

Types of parameters
-------------------

* LibPFASST static parameters:  are hard coded and cannot be changed at run time
* ``pf_pfasst_t`` mandatory parameters: must be reassigned at run time because 
  the use of default values will result in program termination.
* ``pf_pfasst_t`` optional parameters: can be reassigned at run time, but
  the use of default values will result in default execution.
* ``pf_level_t`` mandatory parameters: must be reassigned at run time because
  the use of default values will result in default execution.
* ``pf_level_t`` optional parameters: can  be reassigned at run time, but
  the use of default values will result in default execution.
* local mandatory parameters: must be passed in a call to ``pf_run_pfasst``
* local optional parameters: specified by the user application and
  unrelated to the workings of LibPFASST


pfasst static parameters
---------------------------

The parameters at the top of the file ``src/pf_dtype.f90`` are all set
at compile time and can't be changed at runtime.  The only parameter
here of interest to the user is

.. code-block:: fortran

  integer, parameter :: pfdp = selected_real_kind(15, 307)  

which controls the precision of all floating point numbers (or at
least all those using ``pfdp`` in the declaration).

Mandatory pfasst parameters
---------------------------

The parameters defined in type ``pf_pfasst_t`` in ``src/pf_dtype.f90``
are all given a default value.  Currently only the variable
``nlevels`` is given a problematic default.  Hence setting this
variable on the command line or in an initialization file is mandatory

.. code-block:: fortran

  !>  The main PFASST data type which includes pretty much everything
  type :: pf_pfasst_t
     !> === Mandatory pfasst parameters (must be set on command line or input file)  ===
     integer :: nlevels = -1             !! number of pfasst levels



Optional pfasst parameters
--------------------------

The remaining variables in the specification of ``pf_pfasst_t`` are given
default values as below:

.. code-block:: fortran

     !>  ===  Optional pfasst parameters ====
     integer :: niters  = 5             !! number of PFASST iterations to do
     integer :: qtype   = SDC_GAUSS_LOBATTO  !! type of nodes
     logical :: use_proper_nodes =  .false.
     logical :: use_composite_nodes = .false.
     logical :: use_no_left_q = .false.

     ! --  level dependent parameters
     integer :: nsweeps(PF_MAXLEVS) = 1       !!  number of sweeps at each levels
     integer :: nsweeps_pred(PF_MAXLEVS) =1   !!  number of sweeps during predictor
     integer :: nnodes(PF_MAXLEVS)=3          !! number of nodes

     ! --  tolerances
     real(pfdp) :: abs_res_tol = 0.d0   !!  absolute convergence tolerance
     real(pfdp) :: rel_res_tol = 0.d0   !!  relative convergence tolerance

     ! --  predictor options  (should be set before pfasst_run is called)
     logical :: PFASST_pred = .true.    !!  true if the PFASST type predictor is used
     logical :: pipeline_pred = .false. !!  true if coarse sweeps after burn in are pipelined  (if nsweeps_pred>1 on coarse level)
     integer :: nsweeps_burn =  1       !!  number of sdc sweeps to perform during coarse level burn in
     integer :: q0_style =  0           !!  q0 can take 3 values
                                        !!  0:  Only the q0 at t=0 is valid  (default)
                                        !!  1:  The q0 at each processor is valid
                                        !!  2:  q0 and all nodes at each processor is valid


     ! --  run options  (should be set before pfasst_run is called)
     logical :: Vcycle = .true.         !!  decides if Vcycles are done
     logical :: use_pysdc_V = .false.         !!  decides if Vcycles are done
     logical :: sweep_at_conv = .false. !!  decides if one final sweep after convergence is done
     logical :: Finterp = .false.    !!  True if transfer functions operate on rhs
     logical :: use_LUq = .true.     !!  True if LU type implicit matrix is used
     logical :: use_Sform = .false.  !!  True if Qmat type of stepping is used
     integer :: taui0 = -99          !! iteration cutoff for tau inclusion


     ! -- RK and Parareal options
     logical :: use_sdc_sweeper =.true.  !! decides if SDC sweeper is used (can be turned off for pure parareal)
     logical :: use_rk_stepper = .false. !! decides if RK steps are used instead of the sweeps
     integer :: nsteps_rk(PF_MAXLEVS)=3  !! number of runge-kutta steps per time step
     logical :: RK_pred = .false.        !!  true if the coarse level is initialized with Runge-Kutta instead of PFASST

     ! -- misc
     logical :: debug = .false.         !!  If true, debug diagnostics are printed

     ! -- controller for the results 
     logical :: save_residuals = .true.  !!  Will save residuals every time they are set
     logical :: save_delta_q0 = .true.   !!  Will save change in initial condition
     logical :: save_errors  = .true.    !!  Will save errors, but set_error must be called externally
     integer :: save_timings  = 2        !!  0=none, 1=total only, 2=all, 3=all and echo


Mandatory level parameters
--------------------------

In the specification of ``pf_level_t``, the  variable
``mpibuflen`` is given the default values -1.  This must be
changed to the length of the mpi buffer required by the user's data type
per level as the levels are  initialized.

.. code-block:: fortran

     !  ===Mandatory level parameters===
     integer  :: mpibuflen    = -1   !! size of solution in pfdp units


Optional level parameters
-------------------------
 In the specification of ``pf_level_t``, there are a set of parameters that are given values from the corresponding copies in the ``pf_pfasst_t``.  This redundancy is simply for convenience.  
.. code-block:: fortran

     !  level parameters set by the pfasst_t values
     integer  :: index        = -1   !! level number (1 is the coarsest)
     integer  :: nnodes       = -1   !! number of sdc nodes
     integer  :: nsteps_rk    = -1   !! number of rk steps to perform
     integer  :: nsweeps      = -1   !! number of sdc sweeps to perform
     integer  :: nsweeps_pred = -1      !! number of coarse sdc sweeps to perform predictor in predictor
     logical     :: Finterp = .false.   !! interpolate functions instead of solutions


Mandatory local parameters
--------------------------

In the call to run pfasst

.. code-block:: fortran

  pf_pfasst_run(pf, q0, dt, tend, nsteps, qend)

The variables ``q0``, ``dt``, and ``tend`` must be included.  The
variable ``nsteps`` is optional, if it is not included, then
``nsteps`` is set to

.. code-block:: fortran

  pf%state%nsteps = ceiling(tend/dt)

Finally, ``qend`` is also optional and returns the final solution.

..

  File input for user variables
  -----------------------------

  The  default input file is "probin.nml" wherein the namelist
  PARAMS (defined locally in probin.f90) can be specified.
  Alternatively, a different input file can be specified on the command
  line by adding the file name directly after the executable.  The
  alternative input file must be specified first before any command line
  parameter specifications (see next section).


File input for pfasst  variables
--------------------------------

The pfasst parameters are specified in a namelist ``PF_PARAMS``
defined in routine ``pf_read_opts`` in ``src/pf_options.f90``.  This
routine is called from ``pf_pfasst_create`` in ``pf_pfasst.f90``
(which is typically called when initializing PFASST).  If no file is
specified in the call to ``pf_pfasst_create``, then no file is read.
Typically the main routine specifies this input file (the default
being probin.nml), and this file can be changed by specifying the
value of

  pfasst_nml = "my_file.nml"

either in the local input file or the command line  pfasst_nml=\\\\"my_file.nml\\\\"


Command line input
------------------

All the variables in the namelist ``PF_PARAMS`` can be modified by
simply specifying their value on the command line.  There is only one
caveat to this in that any parameters must be specified after the
(optional) input file specification.  For example

.. code-block:: sh

  mpirun -n 20 main.exe  myinput.nml niters=10

would set the input file to "myinput.nml" and then over-ride any
specified value of niters with the value 10. Command line options
over-ride input files values.


Variables for the predictor
---------------------------

The two variables ``pipeline_pred`` and ``PFASST_pred`` determine how the
predictor works.  The different combinations of these variables and
the parameter nsweeps_pred on the coarsest level great some subtle
differences in how the predictor performs.

Some cases:

#. If PFASST_pred is false and pipeline_pred is false, then the predictor
   is a serial application of SDC with nsweeps_pred sweeps.  This can be done
   without communication wherein every processor mimics the behavior
   of the processors previous to it in time.

#. If PFASST_pred is false and pipeline_pred is true and nsweeps_pred is one,
   then the predictor is a serial application of SDC with 1 sweep.  As
   above, there is no communication necessary.

#. If PFASST_pred is false and pipeline_pred is true and nsweeps_pred is
   greater than one, then the predictor is a version of pipelined
   SDC. There is no communication necessary until the second sweep on
   each processor is done.  After that, each processor must
   recieve a new initial value before each new sweep.

#. If PFASST_pred is true, and nsweeps_pred equals one, then it doesn't
   matter what pipeline_pred is.  No communication is necessary, and we
   simply reuse the function values from the previous iteration in
   each SDC sweep.  Some care must be taken here as to how to
   interpret the variable t0 especially in light of time dependent
   boundary conditions.  Currently t0 does not change in these
   iterations, hence one should use caution using PFASST_pred = true
   with time dependent boundary conditions.

#. If PFASST_pred is true, and nsweeps_pred is greater than one and
   pipeline_pred is true, then the predictor will act like the normal
   PFASST_pred with nsweeps equals one, but more iterations will be
   taken.  This choice is a bit strange.  No communication is needed
   until each processor is doing the P+1st iteration, then new initial
   data must be used and in all cases, previous f values are used in
   the SDCsweeps.  The caveat about t0 is still valid.

#. Finally, if PFASST_pred is true, and nsweeps_pred is greater than one
   and pipeline_pred is false, then the predictor will act like the
   normal PFASST_pred with nsweeps equals one, but additional
   iterations are taken before the initial conditions at each
   processor are reset.  This can be done without communication.  The
   caveat about t0 is still valid.

How is this implemented?  There are two pieces to the initialization.
The first consists of the process of giving every processor an initial
value which is consistent with t0 at that processor.  This can be done
without communication in all cases.  Then,  additional sweeps are done if specified with communication if necessary.
