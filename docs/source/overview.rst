
Overview
========

The LibPFASST library evolves systems of ODEs in time:

.. math::

   u'(t) = f(u(t), t).

LibPFASST comes with builtin sub-steppers for fully explicit, fully
implicit, and implicit-explicit (IMEX) systems.

The user provides:

#. Function evaluation routines: :math:`F(U, t)`.  These functions may
   depend on the PFASST level, which allows for multi-order and/or
   multi-physics operators.

#. For implicit or IMEX schemes, a backward Euler solver: :math:`U -
   dt F(U) = RHS`.

#. Spatial interpolation and restriction routines: :math:`R(U)` and
   :math:`P(U)`.

#. Solution encapsulation routines to tell LibPFASST how to manipulate
   solutions.  This enables LibPFASST to interoperate nicely with
   other libraries that have particular ways of storing solutions (eg,
   BoxLib, PEPC, or PETSc).  In particular, this allows LibPFASST to
   work with spatially distributed solutions.


User interactions with LibPFASST are typically marshaled through the
``type(pf_pfasst_t)`` class.  This class acts as the overall
controller of the algorithm.  Implementing a ODE/PDE solver in
LibPFASST generally consists of the following steps:

#. Create an ``type(pf_pfasst_t)`` controller.

#. Create an ``type(pf_encap_t)`` encapsulation object.  This
   may be level dependent.

#. Create an ``type(pf_comm_t)`` communicator object.  This allows
   LibPFASST to use MPI parallelization (independent of any spatial
   parallelization).

#. Create an ``type(pf_sweeper_t)`` SDC sweeper.  Builtin
   sweepers include: explicit, implicit, or imex.  Custom sweepers can
   also be built.  This can also be level dependent.

#. Add levels to the ``type(pf_pfasst_t)`` controller.  This
   includes telling LibPFASST how many degrees-of-freedom are on each
   level, how many SDC nodes should be used on each level and their
   type, which interpolation and restriction routines to use, which
   encapsulation object to use etc.

#. Add any desired hooks to the controller and set any controller
   options.

#. Call the communicator's and controller's setup routines.

#. Run.

#. Tidy up.
