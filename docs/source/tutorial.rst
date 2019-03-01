Tutorial
========

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


The main program is in ``src/main.f90``.  The IMEX sweeper needs 
functions to evaluate each term in the IMEX splitting and
a routine to solve an implicit equation equivalent to a
backward-Euler step.  These routines are in ``src/feval.f90`` and are called
``f_eval`` and ``f_comp``.

  
Please see the `adv_diff_fft`_ example included in LibPFASST for a simple PDE application of LibPFASST.

This example solves a 1d linear advection diffusion equation

.. math::

  u_t  = - v u_x + \nu u_{xx}.

This right hand side of the equation will be split into stiff terms handled implicitly
(:math:`\nu u_{xx}`) and non-stiff terms handled explicitly (:math:`-v u_x`),
hence an IMEX SDC substepper will be used to evolve the equation through time.

The solution of  implicit equation is done using the FFT through the FFTW package.

Three PFASST levels will be used.  The coarsest level will have 2 SDC
nodes and 32 spatial points; the intermediate level will have 3 SDC
nodes and 64 spatial points; and the fine level will have 5 SDC nodes
and 128 spatial points.  The routines to spatially interpolate and
restrict solutions are in ``src/transfer.f90``.  Two coarse sweeps
will be performed per iteration on the coarsest levels.  This helps
reduce the total number of PFASST iterations required.

The solution :math:`u` will be stored in a flat Fortran array, and
hence this application will use LibPFASSTs built in ``ndarray``
encapsulation.  Note that LibPFASST doesn't impose any particular
storage format on your solver -- instead, you must tell LibPFASST how
to perform a few basic operations on your solution (eg, how to create
solutions, perform ``y <- y + a x``, etc).  Various hooks are added to
echo the error (on the finest level) and residual (on all levels)
throughout the algorithm.  These hooks are in ``src/hooks.f90``.

Note that the ``feval_create_workspace`` routine is specific to the
problem being solved (ie, not part of LibPFASST, but part of the user
code).  It creates FFTW plans, creates a complex workspace for use
with FFTW (so that we can avoid allocating and deallocating these
workspaces during each call to the function evaluation routines), and
pre-computes various spectral operators.

LibPFASST allows you, the user, to attach an arbitrary C pointer to
each PFASST level.  This is called a context (typically called
``levelctx`` in the source) pointer (as in, "for the problem I'm
solving I have a specific context that I will be working in").  Your
context pointer gets passed to your function evaluation routines and
to your transfer routines.  In most of the examples this context
pointer is used to hold FFTW plans etc as described above.  Note that
each level gets its own context because each level has a different
number of degrees-of-freedom (``nvars``).

C pointers are used because they provide a lot of flexibility.  The
drawback to this is that we loose the ability for the compiler to do
type checking for us.


