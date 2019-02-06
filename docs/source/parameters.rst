
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

* libpfasst static parameters:  hard coded and cannot be changed at run time.
* ``pf_pfasst_t`` mandatory parameters: must be reassigned at run time,
  the use of default values will result in program termination.
* ``pf_pfasst_t`` optional parameters: can be reassigned at run time,
  the use of default values will result in default execution.
* ``pf_level_t`` optional parameters: can  be reassigned at run time,
  the use of default values will result in default execution.
* ``pf_level_t`` mandatory  parameters:   must be reassigned at run time,
  the use of default values will result in program termination.
* local mandatory parameters: must be passed in a call to ``pf_run_pfasst``.
* local optional parameters: specified by the user application and
  unrelated to the workings of libpfasst.


Libfpasst static parameters
---------------------------

The parameters at the top of the file `src/pf_dtype.f90` are all set
at compile time and can't be changed at runtime.  The only parameters
here of interest to the user are

.. code-block:: fortran

  integer, parameter :: pfdp = selected_real_kind(15, 307)  !!  Defines double precision type for all real and complex variables
  integer, parameter :: pfqp = selected_real_kind(33, 4931) !!  Defines quad precision type for all real and complex variables
  
which controls the precision of all floating point numbers  using ``pfdp`` or ``pfqp`` in the declaration and assignment (this is highly encouraged).
The ``selected_real_kind`` function is a Fortran intrinsic which is designed to allow for differing definition of precision for different compilers.

In theory, one can run libpfasst in quad precision by changing the first line to 

.. code-block:: fortran
		
  integer, parameter :: pfdp = selected_real_kind(33, 4931)  !! For quad precision everywhere (use at your risk and see top of pf_mpi.f90)

Since this will change the size of the data being passed by mpi, one then must also change the parameter

.. code-block:: fortran

   integer, parameter :: myMPI_Datatype=MPI_REAL8

The tutorial examples EX1_Dahlquist and  EX2_Dahlquist can be run in quad precision and will demonstrate quad precision in the residual and errors, but not all aspects of libpfasst will
work as desired with quad precision.  For example, some of the Runge-Kutta coefficients are not known to quad precision.  Hence, quad precision should be used at your own risk.

   
Pfasst mandatory parameters
--------------------------
The first variable in the specification of ``pf_pfasst_t`` in `pf_dtype.f90` is the only variable that is mandatory to set.  

.. code-block:: fortran
		
     !>  Mandatory parameters (must be set on command line or input file)
     integer :: nlevels = -1             !! number of pfasst levels

It can be set either on the command line, or more usually in an input file using the  `PF_PARAMS`  namelist (see below).


Pfasst optional parameters
--------------------------

All the remaining variables in the specification of ``pf_pfasst_t`` in `pf_dtype.f90` are given default values as below:

.. code-block:: fortran

     !>  Optional parameters
     integer :: niters  = 5             !! number of PFASST iterations to do
     integer :: qtype   = SDC_GAUSS_LOBATTO  !! type of nodes
     
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
     logical :: Finterp = .false.    !!  True if transfer functions operate on rhs
     logical :: use_LUq = .true.     !!  True if LU type implicit matrix is used 
     integer :: taui0 = -999999     !! iteration cutoff for tau inclusion


     !> RK and Parareal options
     logical :: use_rk_stepper = .false. !! decides if RK steps are used instead of the sweeps
     integer :: nsteps_rk(PF_MAXLEVS)=3  !! number of runge-kutta nodes
     logical :: RK_pred = .false.        !!  true if the coarse level is initialized with Runge-Kutta instead of PFASST

     ! -- misc
     logical :: debug = .false.         !!  If true, debug diagnostics are printed
     logical :: save_results = .false.  !!  If true, results are output
     logical    :: echo_timings  = .false.    !!  If true, timings are output


These values can be changed as desired either on the command line or in an input file as described below.
Except for the predictor parameters, the meaning  of most of the parameters should be clear from the context and from reading the
description of the pfasst algorithm.  See the section on  the predictor for more discussion of the predictor parameters.

Level mandatory parameters
--------------------------
There is one level parameter that must be set on each level by the user, namely mpibuflen, which gives the size of the
solution that is communicated by MPI.  There is no way for libpfasst to know the value of this parameter, so the code
will terminate if it is not set before the call to pfasst_setup in the user's main.

.. code-block:: fortran
		
     !  Mandatory level parameter
     integer  :: mpibuflen    = -1   !! size of solution in pfdp units

     
Level optional parameters
-------------------------

In the specification of ``pf_level_t`` in `pf_dtype.f90`, the first six parameters are assigned values in
the subroutine ``pf_pfasst_create`` located in `pf_pfasst.f90`.
Except for the first, these can be changed per level as the levels are initialized in the users main routine

.. code-block:: fortran

  type :: pf_level_t

     integer  :: index        = -1   !! level number (1 is the coarsest)
     integer  :: nnodes       = -1   !! number of sdc nodes
     integer  :: nsteps_rk    = -1   !! number of rk steps to perform
     integer  :: nsweeps      = -1   !! number of sdc sweeps to perform
     integer  :: nsweeps_pred = -1      !! number of coarse sdc sweeps to perform predictor in predictor
     logical     :: Finterp = .false.   !! interpolate functions instead of solutions





Local mandatory parameters
--------------------------

In the call to run pfasst

.. code-block:: fortran

  pf_pfasst_run(pf, q0, dt, tend, nsteps, qend, flags)

The variables ``q0``, ``dt``, and ``tend`` must be included.  These correspond to the initial condition, the time step, and the end time of the run.

The variable ``nsteps`` is optional, if it is not included, then ``nsteps`` is set to

.. code-block:: fortran

  pf%state%nsteps = ceiling(tend/dt)

If it is included, then the value of ``tend``  passed into the routine is ignored and the final time of the simulation will be ``nsteps*dt``

The input paratmer ``qend`` is also optional and returns the final solution if desired.  Finally, an integer array ``flags`` can be passed if desired.

..

File input for user variables
-----------------------------

  The usual default input file for libpfasst examples is  `probin.nml` wherein the namelist
  PARAMS (defined locally in `probin.f90`) can be specified.
  Alternatively, a different input file can be specified on the command
  line by adding the file name directly after the executable.  The
  alternative input file must be specified first before any command line
  parameter specifications (see next section).  For a given application, there is no requirement that the program reads in any local parameters, and the style of `probin.f90` can be changed to anything else.
  It is necessary however to provide an input for pfasst variables described next.


File input for pfasst  variables
--------------------------------

The pfasst parameters are specified in a namelist ``PF_PARAMS``
defined in routine ``pf_read_opts`` in `pf_pfasst.f90`.  This
routine is called from ``pf_pfasst_create`` in `pf_pfasst.f90`
(which is typically called when initializing PFASST).  If no file is
specified in the call to ``pf_pfasst_create``, then no file is read.
Typically the main routine specifies this input file (the default
being probin.nml), and this file can be changed by specifying the
value of

  pfasst_nml = 'probin.nml'

either in the local input file  or the command line.

This is not completely transparent, so consider some cases:

*  I include an input file on the command line and it contains  a ``PF_PARAMS`` namelist: This is fine as long as ``PF_PARAMS`` has an ``nlevels`` entry
*  I include no input file on the command line:  Then the input file `probin.nml` will be read for the namelist and two possibilities exist.
   * `probin.nml` has a ``PF_PARAMS`` namelist including an ``nlevels`` entry
   * `probin.nml` has  an assignment of a different `pfasst_nml` input file in which 
   * `probin.nml` has  no namelist but  ``nlevels`` is specified on the command line




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
over-ride input files.


Variables for the predictor
---------------------------
Warning:  This section may not be current due to an increase in the possible ways the predictor is called. The interested reader
might look directly in the source code in the file `src/pf_parallel.f90`

The two variables ``pipeline_pred`` and ``PFASST_pred`` determine how the
predictor works.  The different combinations of these variables and
the parameter ``Nsweeps_pred`` create some subtle
differences in how the predictor performs.

Some cases:

#. If ``PFASST_pred`` is false and ``pipeline_pred`` is false, then the predictor
   is a serial application of SDC with Nsweeps.  This can be done
   without communication wherein every processor mimics the behavior
   of the processors previous to it in time.

#. If ``PFASST_pred`` is false and ``pipeline_pred`` is true and ``Nsweeps`` is one,
   then the predictor is a serial application of SDC with 1 sweep.  As
   above, there is no communication necessary.

#. If ``PFASST_pred`` is false and ``pipeline_pred`` is true and ``Nsweeps`` is
   greater than one, then the predictor is a version of pipelined
   SDC. There is no communication necessary until the second sweep on
   the each processor is done.  After that, each processor must
   recieve a new initial value.

#. If ``PFASST_pred`` is true, and ``Nsweeps`` equals one, then it doesn't
   matter what pipeline_pred is.  No communication is necessary, and we
   simply reuse the function values from the previous iteration in
   each SDC sweep.  Some care must be taken here as to how to
   interpret the variable t0 especially in light of time dependent
   boundary conditions.  Currently t0 does not change in these
   iterations, hence one should use caution using PFASST_pred = true
   with time dependent boundary conditions.

#. If ``PFASST_pred`` is true, and ``Nsweeps`` is greater than one and
   ``pipeline_pred`` is true, then the predictor will act like the normal
   ``PFASST_pred`` with ``Nsweeps`` equal one, but more iterations will be
   taken.  This choice is a bit strange.  No communication is needed
   until each processor is doing the P+1st iteration, then new initial
   data must be used and in all cases, previous f values are used in
   the SDCsweeps.  The caveat about t0 is still valid.

#. Finally, if ``PFASST_pred`` is true, and ``Nsweeps`` is greater than one
   and ``pipeline_pred`` is false, then the predictor will act like the
   normal PFASST_pred with ``Nsweeps`` equals one, but additional
   iterations are taken before the initial conditions at each
   processor are reset.  This can be done without communication.  The
   caveat about t0 is still valid.

