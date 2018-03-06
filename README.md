# LIBPFASST

LIBPFASST is a lightweight implementation of the Parallel Full
Approximation Scheme in Space and Time (PFASST) algorithm.  It is
written in Fortran (mostly F90, with some F03), but can be interfaced
with C and C++ fairly easily.

# NOTICE ON NEW RELEASE

You might notice that the code here has undergone substantial recent changes.
We are currently in the process of updating the documentation and examples to reflect the changes.
Please stay tuned.

You can see the old documentation in the 'docs' directory and,
an online version of the documentation is available at

  https://libpfasst.readthedocs.org/en/latest/

There is (new) Doxygen info at

https://pfasst.bitbucket.io/Doxygen/html/index.html

<!---
and examples in 'examples' for
usage information.
--->



# COMPILERS

LIBPFASST has been successfully compiled with

  + GNU 4.5 and 4.8,
  + Intel 13.0.1, and
  + PGI 12.9 (almost, minor issues, but overall OK).
   

Before compiling, you will need an MPI library installed.  On Debian
or Ubuntu machines, we recommend MPICH2:

  `$ sudo apt-get install mpich2 libmpich2-dev`


# COMPILING

The Makefile included in the LIBPFASST distribution uses the following
defaults:

  + Fortran: 'mpif90'
  + C:       'gcc' 

You can change these defaults by creating (or symlinking, see
mk/*.defs) a .defs file (eg, Makefile.defs) in the LIBPFASST
directory.

Once you have configured your compilers, you can build LIBPFASST by
running 'make' in the LIBPFASST directory.  The resulting (static)
library should be in the 'build' directory.


# EXAMPLES

Once LIBPFASST is compiled, you can build the examples individually.
The 'mpi-advection' and 'pth-advection' examples require FFTW3.  If
you already have FFTW3 installed on your machine, please set the FFTW3
environment variable appropriately.  For example:

  `$ export FFTW3=$HOME/opt`

If you do not have FFTW3 installed, simply run the 'fftw3' make rule
from the LIBPFASST directory:

  `$ make fftw3`

This will download and compile a local (static) copy of FFTW3.  Once
you have FFTW installed, you can build the mpi-advection and
pth-advection examples by:

  `$ cd examples/mpi-advection`
  `$ make`

Similarly, some other examples may require HDF5.  If you do not
already have HDF5 installed, simply run the 'hdf5' make rule from the
LIBPFASST directory:

  `$ make hdf5`

This will download and compile a local (static) copy of HDF5.


# LICENSE

LIBPFASST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

LIBPFASST is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.


# TAU

Note that using TAU with MPICH 3.X doesn't work.  MPICH 1.5 seems to be fine.

Here is how memmett installed TAU using MPI and OpenMP:

~~~
  $ H=/home/memmett
  $ O=$H/opt
  $ export PATH=$H/bin:$O/bin:$H/gcc-4.7/bin
  $ export PATH=$PATH:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
  $ ./configure -mpi -mpiinc=$O/include -mpilib=$O/lib -openmp -useropt="-D__USE_XOPEN2K8"
~~~

Note the ordering of the paths above: we want TAU to pick up our GCC,
not the system GCC.  The ouput of the TAU configure scripts should
show that TAU detected GCC 4.7.3.

Next:

~~~
  $ export PATH=$H/src/tau-2.22.2/x86_64/bin:$PATH
  $ make install
~~~

To use TAU on libpfasst programs, set the following in your *.defs file:

~~~
  export TAU_MAKEFILE=/home/memmett/src/tau-2.22.2/x86_64/lib/Makefile.tau-mpi-openmp
  export MPICH_F90=/home/memmett/gcc-4.7/bin/gfortran
  FC = tau_f90.sh
  CC = tau_cc.sh
~~~

Some useful environment variables:

~~~
  $ export TAU_TRACE=1
  $ export TAU_CALLPATH=1
  $ export TAU_CALLPATH_DEPTH=30
~~~

To bring up Jumpshot:

~~~
  $ mpirun -n 4 ./main.exe
  $ tau_treemerge.pl
  $ tau2slog2 tau.trc tau.edf -o app.slog2
  $ jumpshot app.slog2
~~~