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
