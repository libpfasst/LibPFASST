Tutorial
========

Please see the `mpi-advection`_ example for a simple example of how to
use LIBPFASST.

This example solves the linear advection diffusion equation

.. math::

  u_t + v u_x = \nu u_xx.

This equation will be split into stiff (:math:`\nu u_xx`) and
non-stiff (:math:`-v u_x`) pieces, and an IMEX SDC substepper will be
used.


.. _`mpi-advection`: https://bitbucket.org/memmett/libpfasst/src/ffa8de3a1feafa7e8f811da0a418b939aca3e2c0/examples/mpi-advection?at=master
