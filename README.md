# LIBPFASST

LIBPFASST is a lightweight implementation of the Parallel Full
Approximation Scheme in Space and Time (PFASST) algorithm.  It is
written in Fortran (mostly F90, with some F03), but can be interfaced
with C and C++ fairly easily.

# NOTICE ON NEW RELEASE

You might notice that the code here has undergone substantial recent changes.
We are currently in the process of updating the documentation and examples to reflect the changes.
Please stay tuned.

You can see the old documentation in the 'docs' directory, and
an online version of the documentation is available at

  https://libpfasst.readthedocs.org/en/latest/

There is (new) Doxygen info at

https://pfasst.bitbucket.io/Doxygen/html/index.html



# COMPILERS

LIBPFASST has been successfully compiled with

  + GNU 4.5 and 4.8,
  + Intel 13.0.1, and
  + PGI 12.9 (almost, minor issues, but overall OK).
   

Before compiling, you will need an MPI library installed.  On Debian
or Ubuntu machines, we recommend MPICH2:

  `$ sudo apt-get install mpich2 libmpich2-dev`


# COMPILING

Before compiling, you should check the definitions of the Fortran (FC)  and C (CC) compiler in Makefile.defaults
and change them if necessary

+ FC =mpif90
+ CC = mpicc

Once you have configured your compilers, you can build LIBPFASST by
running 'make' in the LIBPFASST directory.  The resulting (static)
library should be in the 'lib' directory.


