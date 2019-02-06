Compiling
=========

Libpfasst comes with its own Makefiles, and it is possible that simply entering

   `~/libpfasst $ make`

in the directory /libpfasst will compile the code succesfully.  A successful make will build the file

libpfasst/lib/libpfasst.a  

Before compiling, you will need an MPI library installed.  On Debian
or Ubuntu machines, we recommend MPICH2:

  `$ sudo apt-get install mpich2 libmpich2-dev`

The user will probably have to adjust some flags in 'Makefile.defaults'.

The first three flags correspond to the compilers and linkers.  These can be changed to correspond to those used on your system.

FC = mpif90

CC = mpicc

AR = ar rcs

Some optional flags are

DEBUG = FALSE

which if made true will add many diagnostic flags to the build, and

MKVERBOSE=FALSE

which if made true will show more info in the build steps.  These can also both be changed on the command line as in

   `~/libpfasst $ make MKVERBOSE=TRUE DEBUG=TRUE`

There is is a Fortran dependency file included called .depend.  If you want to remake this, it can be done using the makedpef90 package by uncommenting out the appropriate lines in Makefile.rules and entering

   `~/libpfasst $ make depend`

Finally, enter

   `~/libpfasst $ make clean`

to remove all intermediate files and start from scratch.

