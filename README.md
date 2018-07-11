
# LIBPFASST

LIBPFASST is a lightweight implementation of the Parallel Full
Approximation Scheme in Space and Time (PFASST) algorithm.  It is
written in Fortran (mostly F90, with some F03), but can be interfaced
with C and C++ fairly easily.


Libpfasst Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) and Sebastian Goetschel. All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. IBPFASST

# NOTICE ON NEW RELEASE

We have changed the License Agreement.  Please see LICENSE for details

There is (newish) Doxygen info at
https://pfasst.bitbucket.io/Doxygen/html/index.html

We are currently writing better documentation, examples, and a tutorial,  https://pfasst.bitbucket.io/docs/build/html/index.html


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


