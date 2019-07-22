Sweepers
========

Under construction!

LibPFASST contains several sweeper types for solving different types of equations.  These are contained in ``src/sweepers/``.


* exp:  Exponential integrator
* imex:  Semi-implicit or IMEX integrator for equations with an implicit  and explicit splitting
* misdc: Multi-implicit integrator for equations with an two implicit  and 1 explicit splitting
* imk:  Implicit Munthe-Kass Runge-Kutta integrators 
* magpicard:   Magnus  integrator based on Picard iterations
* verlet:   Verlet type integrator for 2nd order equations
* imex_oc:  Optimal control version of imex
* misdc_oc:  Optimal control version of misdc

For historical reasons, there are two versions of some of these sweeper, with a `Q` denoting the more modern version that can utilize DIRK style sweepers.  
