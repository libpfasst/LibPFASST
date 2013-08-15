Tutorial
========

Please see the ``mpi-advection`` example included in LIBPFASST for a
simple application of LIBPFASST.

This example solves a 1d linear advection diffusion equation

.. math::

  u_t + v u_x = \nu u_{xx}.

This equation will be split into stiff (:math:`\nu u_{xx}`) and
non-stiff (:math:`-v u_x`) pieces, and hence an IMEX SDC substepper
will be used to evolve the equation through time.

The main program is in ``src/main.f90``.  The IMEX sweeper needs three
functions: one to evaluate the non-stiff piece, another to evaluate
the stiff piece, and one more to solve a backward-Euler step.  These
routines are in ``src/feval.f90``.

Three PFASST levels will be used.  The coarsest level will have 2 SDC
nodes and 32 spatial points; the intermediate level will have 3 SDC
nodes and 64 spatial points; and the fine level will have 5 SDC nodes
and 128 spatial points.  The routines to spatially interpolate and
restrict solutions are in ``src/transfer.f90``.  Two coarse sweeps
will be performed per iteration on the coarsest levels.  This helps
reduce the total number of PFASST iterations required.

The solution :math:`u` will be stored in a flat Fortran array, and
hence this application will use LIBPFASSTs built in ``ndarray``
encapsulation.  Various hooks are added to echo the error (on the
finest level) and residual (on all levels) throughout the algorithm.
These hooks are in ``src/hooks.f90``.

Note that the ``feval_create_workspace`` routine is specific to the
problem being solved (ie, not part of LIBPFASST, but part of the user
code).  It creates FFTW plans, creates a complex workspace for use
with FFTW (so that we can avoid allocating and deallocating these
workspaces during each call to the function evaluation routines), and
pre-computes various spectral operators.

LIBPFASST allows you, the user, to attach an arbitrary C pointer to
each PFASST level.  This is called a context (typically called ``ctx``
in the source) pointer (as in, "for the problem I'm solving I have a
specific context that I will be working in").  Your context pointer
gets passed to your function evaluation routines and to your transfer
routines.  In most of the examples this context pointer is used to
hold FFTW plans etc as described above.  Note that each level gets its
own context because each level has a different number of
degrees-of-freedom (``nvars``).

C pointers are used because they provide a lot of flexibility.  The
drawback to this is that we loose the ability for the compiler to do
type checking for us.
