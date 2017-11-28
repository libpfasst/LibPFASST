# LIBPFASST Tasks

1. Understand the logic of Matt's communication routine so that we can eventually remove all the nonworking stuff.
- Document the logic of communication so that everybody can remember what it does.
- Change the possibilities for status and pstatus to converged or not_converged
  or just a boolean.
- Have the sweepers not sweep if the residual is small enough.  Note that a
  sweeper may choose not to sweep on a given iteration but then receive a new
  initial condition that requires further sweeps.  The residual may change because
  the solution changes or the initial condition changes.  I still wonder if the
  residual computation should always be the first thing in the sweeper.  I also
  wish the code would not recompute the first function value if it hasn't changed,
  but this logic is also tricky.
- Separate the posting and sending logic from the "check_convergence" routine.
- Remove the logic from inside of communication routines.
- Revisit consolidating the traditional Euler based sweepers
- Make unit tests that do all relevant permutations on
  - Number of levels
  - Number of processors
  - Vcycle or not Vcycle
  - Type of sweeper
  - if tolerances are set or just zero
  - predictor types
  - node types?
- Eventually merge branches nwc with master.  Why are they still different?
