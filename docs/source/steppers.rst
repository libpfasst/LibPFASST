Steppers
========

Under construction!

LibPFASST contains two stepper types which execute Runge-Kutta time steps over an interval.
These are usually applied during a Parareal implementation, but can also be called to initialize the coarsest level in PFASST
These are contained in `src/` within module files that contain `sweeper` in the name.

* ``rkstepper``: IMEX Runge-Kutta integrator for equations with an implicit  and explicit splitting
* ``erkstepper``: Exponential Runge-Kutta integrator for equations that can be split into
  linear and nonlinear terms
