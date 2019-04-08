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

Using your favorite editor, open the file sdc.nml.  There are two namelists here, 'PF_PARAMS' and "PARAMS'.  The second of these is for local variables which in this simple example are just the values of lam1 and lam2, the simulation time and the number of steps.  The 'PF_PARAMS' variables are discussed below in the section `Parameters <parameters>`_
Any parameter in either list can be overwritten by adding it to the command line, as in

`$ ./main.exe sdc.nml lam1=3.0 niters=10`

The order of the command line parameters is not important except that they must come after the input file is specified.

To run an example that does the actualy PFASST algorithm

`$ mpirun -n 32 ./main.exe multilevel.nml`

The main program is in ``src/main.f90`` and can be used as a template for building applications using LibPFASST.
The most important part of the routine ``run_pfasst`` is loop over levels where the level, sweeper, and
data encapsulation are specified.  The abstract type defining a level must be extended with spatial restriction and
interpolation operators (in this example, these are the identity operators).  This is done in the file ``src/level.f90``.
Next, the encapsulation of the data type must be specified by assigning its factory.  In this example, a data type provided by the library corresponding to an N-dimensional array is used.  In the next Tutorial, a local encapsulation will be demonstrated.  Third, the type abstract sweeper type must be extended and assigned to the level. This is done in the file ``src/sweeper.f90`` discussed below.  Finally, the level type carries a one-dimensional integer array that to carry information about the size of the data.  This array must be set, and here the size is simply a one-dimensional array of length 1.

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
   
The local sweeper type needs to define
functions to evaluate each term in the IMEX splitting and
a routine to solve an implicit equation equivalent to a
backward-Euler step.  These routines are in ``src/sweeper.f90`` and are called
``f_eval`` and ``f_comp``.

There is also a file ``src/probin.f90`` to read in parameters for LibPFASST and the local application.


Example 2
---------
This example solves exactly the same equation as Example 1, but using more user generated code.
The main difference is that instead of using the LibPFASST data encapsulation ndarray, a local data
encapsulation called scalar_encap is defined in ``src/encap.f90``.  To define a data encapsulation, the user must
define 7 simple actions on the data set corresponding the the procedure

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

In this example, these are all trivial and should be self-explanatory from the code.  The last of these, eprint, is not typically needed for by the LibPFASST but is included for convenience in debugging.

Another difference in this example, is that a local hook is defined in the file ``src/hooks.f90`` to print the error to the screen.  The user can construct custom hooks following this template.
		
  

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


