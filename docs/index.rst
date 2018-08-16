.. LIBPFASST documentation master file, created by
   sphinx-quickstart on Tue Jul 10 15:22:48 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the  LIBPFASST documentation!
========================================

LIBPFASST is a Fortran implementation of the PFASST algorithm.

Libpfasst Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) and Sebastian Goetschel.  All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a  non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :download: 
**Main parts of the documentation**

* :doc:`Download <download>` - download and installation instructions.
* :doc:`Compiling <compiling>` - instructions for compiling
* :doc:`Tutorial <tutorial>` - getting started and basic usage.
* :doc:`Overview <overview>` - design and interface overview.
* :doc:`Reference <reference>` - information about the internals of LIBPFASST.

Download
========

You can obtain a copy of LIBPFASST through the `PFASST@lbl.gov project
page`_ on `bitbucket.org`_.  

.. _`LIBPFASST project page`: https://pfasst.lbl.gov
.. _`bitbucket.org`: https://bitbucket.org/berkeleylab/libpfasst/src/master/

Compiling
=========

libpfasst comes with its own Makefiles, and it is possible that simply typing 'make' in the directory /libpfasst will compile the code succesfully.  If not, the user will probably have to adjust some flags in 'Makefile.defaults' as discussed below.  A successful make will build the file libpfasst/lib/libpfasst.a

The first three flags correspond to the compilers and linkers.  These can be changed to correspond to those used on your system.

FC = mpif90
CC = mpicc
AR = ar rcs

Some optional flags are

DEBUG = FALSE

which if made true will add many diagnostic flags to the build, and

MKVERBOSE=FALSE

which if made true will show more info in the build steps.  These can both be changed on the command line as in

 make MKVERBOSE=TRUE DEBUG=TRUE

There is is a Fortran dependency file included called .depend.  If you want to remake this, it can be done using the makedpef90 package by uncommenting out the appropriate lines in Makefile.rules and typing

make depend

Finally, one can type

make clean

to remove all intermediate files and start from scratch.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :hidden:

   self
   maths
   tutorial
   overview
   reference
   communicators
   parameters
