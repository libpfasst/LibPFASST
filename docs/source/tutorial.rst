Tutorial
========

Example 1
---------
The following material will walk the user through a couple of examples to demonstrate how to set up an application using
LibPFASST.


Once libpfasst has been successfully built, move to the directory  libpfasst/Tutorials/EX1_Dahlquist
This example solves the  scalar model problem or Dahlquist equation

.. math::

  y'  = \lambda y

  y(0) = 1

An implicit-explicit or IMEX  (also known as semi-implicit) splitting is used in this example, so the equation can be written 

.. math::

   y'  = \lambda_1 y + \lambda_2 y

Typing 'make' in the directory should compile the example creating an exectuable called 'main.exe'.  In the same directory, there are a few parameter files with the extension '.nml'.  You can run the example using one of these files, as in

  `$ ./main.exe sdc.nml`

Using your favorite editor, open the file sdc.nml.  There are two namelists here, 'PF_PARAMS' and 'PARAMS'.  The second of these is for local variables which in this simple example are just the values of lam1 and lam2, the simulation time and the number of steps.  The 'PF_PARAMS' variables are discussed below in the section `Parameters <parameters>`_. 
Any parameter in either list can be overwritten by adding it to the command line, as in

`$ ./main.exe sdc.nml lam1=3.0 niters=10`

The order of the command line parameters is not important except that they must come after the input file is specified.

To run an example that does the actual PFASST algorithm

`$ mpirun -n 32 ./main.exe multilevel.nml`

The main program is in ``src/main.f90`` and can be used as a template for building applications using LibPFASST.  The main routine here
only initializes MPI, calls a routine ``run_pfasst`` to run the PFASST algorithm, then closes MPI.

The first thing done in ``run_pfasst`` is to read in local problem parameters by calling  the subroutine ``probin_init`` located in ``src/probin.f90``.

.. code-block:: fortran
		
    !> Read problem parameters
    call probin_init(pf_fname)

This routine also returns the location of the namelist for the PFASST parameters, which in this case will match the location of the local parameters (a different file can be specified if desired).  After
the MPI based communicator is set up, the routine ``pf_pfasst_create`` is called to allocate the main PFASST structure called 'pf'.  Note that the filename for PFASST parameters is
based to this routine so that it can be read.

.. code-block:: fortran
		
    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)
    

The most important part of the initialization in ``run_pfasst`` is the loop over levels where the level, sweeper, and data encapsulation are specified.

.. code-block:: fortran

    !> Loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data constructor
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)

       !>  Allocate the sweeper at this level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !>  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,[1])
    end do

The abstract type defining a level must be extended with spatial restriction and
interpolation operators (in this example, these are the identity operators).  This is done in the file ``src/level.f90``.
Next, the encapsulation of the data type must be specified by assigning its factory.  In this example, the data type 'ndarray' provided by LibPFASST and corresponding to an N-dimensional array is used.
In the next Tutorial, a local encapsulation will be demonstrated.  Third, the abstract sweeper type must be extended and assigned to the level. This is done in the file ``src/sweeper.f90`` discussed below.  Finally, the level type carries a one-dimensional integer array to carry information about the size of the data.  This array must be set, and here the size is simply a one-dimensional array of length 1.

   
The local sweeper type needs to define
functions to evaluate each term in the IMEX splitting and
a routine to solve an implicit equation equivalent to a
backward-Euler step.  These routines are in ``src/sweeper.f90`` and are called
``f_eval`` and ``f_comp``.

After the levels are assigned, the rest of the PFASST structure can is made by calling

.. code-block:: fortran

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

Next, a hook is added that will echo residuals to the screen after every iteration.

.. code-block:: fortran
		
    !> add some hooks for output  (using a LibPFASST hook here)
    call pf_add_hook(pf, -1, PF_POST_ITERATION, pf_echo_residual)

After a routine to echo 
the run options to the screen, the initial conditions are set, and then the
routine to actually do the time stepping is called.

.. code-block:: fortran

    !> Do the PFASST time stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps)

The rest is just cleanup.

Example 2
---------
This example solves exactly the same equation as Example 1, but using more user generated code.
The main difference is that instead of using the LibPFASST data encapsulation 'ndarray', a local data
encapsulation called 'scalar_encap' is defined 

.. code-block:: fortran
		
       !>  Allocate the user specific data constructor
       allocate(scalar_factory::pf%levels(l)%ulevel%factory)

The relevant code for the factory is in ``src/encap.f90``.       

.. code-block:: fortran

  !>  Type to create and destroy the local data encapsulation
  type, extends(pf_factory_t) :: scalar_factory
   contains
     procedure :: create_single  => scalar_create_single
     procedure :: create_array  => scalar_create_array
     procedure :: destroy_single => scalar_destroy_single
     procedure :: destroy_array => scalar_destroy_array
  end type scalar_factory

The four required subroutines are in this case trivial since no data structures need to be allocated to make the encapsulation.  

To define a data encapsulation, the user must also provide
7 subroutines that define actions on the data set corresponding to the procedures in
``src/encap.f90``.       

.. code-block:: fortran

  contains
     procedure :: setval => scalar_setval
     procedure :: copy => scalar_copy
     procedure :: norm => scalar_norm
     procedure :: pack => scalar_pack
     procedure :: unpack => scalar_unpack
     procedure :: axpy => scalar_axpy
     procedure :: eprint => scalar_eprint
  end type scalar_encap

In this example, these are all trivial and should be self-explanatory from the code.  The last of these, eprint, is not typically needed  by LibPFASST but is included for convenience in debugging.

The sweeper assigned in this example is the same as in Example 1, but there are two additional routines defined in ``src/sweeper.f90``.

.. code-block:: fortran

     procedure :: initialize !  Overwrites imex sweeper initialize
     procedure :: destroy    !  Overwrites imex sweeper destroy

These two routines will be called instead of the base sweeper initialize and destroy in LibPFASST.  The point is that this then allows the user to add whatever things to the sweeper as necessary.  Here, there is nothing to do, but one must explicitly call the LibPFASST versions of these routines, as in

.. code-block:: fortran

    !>  Call the imex sweeper initialization
    call this%imex_initialize(pf,level_index)


Another difference in this example, is that a local hook is defined in the file ``src/hooks.f90`` to print the error to the screen.  It is assigned by


.. code-block:: fortran
		
    !>  Add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)


The user can construct custom hooks following this template.

Finally, note that in this example, an optional argument to return the solution at the final time, y_end is included in the call to ``pf_pfasst_run``

.. code-block:: fortran
		
    !> Do the PFASST time stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)
  

Example 3
---------
Please see the ``LibPFASST/test/adv_diff_fft`` directory included in LibPFASST
for a simple PDE application of LibPFASST.

This example solves a 1d linear advection diffusion equation

.. math::

  u_t  = - v u_x + \nu u_{xx}.

This right hand side of the equation will be split into stiff terms handled implicitly
(:math:`\nu u_{xx}`) and non-stiff terms handled explicitly (:math:`-v u_x`),
hence an IMEX SDC substepper will be used to evolve the equation through time.

The solution of  implicit equation is done using the FFT through the FFTW package.


to be continued . . .
