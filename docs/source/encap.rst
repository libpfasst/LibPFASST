Data encapsulation
==================

!  Under construction

The temporal integration routines in LibPFASST do not specify the form of the solution in terms of its data layout, but rather this must be described by an encapsulation module which tells LibPFASST how to operate on the data.  For convenience LibPFASST provides sample encapsulations corresponding to regular Fortran arrays. These are provided in
`src/` within module files whose names  include `encap`.

Curently, the provided encaps include

* ``ndarray``:  An N-dimensional array
* ``ndsysarray``: System of ndarrays
* ``zndarray``: An N-dimensional complex array
* ``zndsysarray``: System of complex ndarrays
* ``ndarray_oc``:  Optimal control version of ndarray

A discussion of how to construct an encapsulation is included in the 2nd tutorial.  
The advanced user might also consider the encapsulations

* ``petscVec``:  uses the standard Vec structure from PETSc
* ``hamsys``:    provides encapsulation of two Hamiltonian variables for Verlet sweepers
* ``amrex``:  based on the multifab structure from the AMReX code (under construction)
