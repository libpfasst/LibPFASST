
# LibPFASST

LibPFASST is a lightweight implementation of the Parallel Full
Approximation Scheme in Space and Time (PFASST) algorithm.  It is
written in Fortran (mostly F90, with some F03), but can be interfaced
with C and C++ fairly easily.

## References
- Matthew Emmett, Michael Minion: Toward an efficient parallel in time method for partial differential equations, Commun. Appl. Math. Comput. Sci. 7(1), 105-132, 2012. [http://dx.doi.org/10.2140/camcos.2012.7.105](http://dx.doi.org/10.2140/camcos.2012.7.105)
- Matthew Emmett, Michael Minion: Efficient Implementation of a Multi-Level Parallel in Time Algorithm. In: Erhel, J., Gander, M., Halpern, L., Pichot, G., Sassi, T., Widlund, O. (eds) Domain Decomposition Methods in Science and Engineering XXI. Lecture Notes in Computational Science and Engineering, vol 98. Springer, Cham, 2014. [https://doi.org/10.1007/978-3-319-05789-7_33](https://doi.org/10.1007/978-3-319-05789-7_33)

##  Documentation
We are currently writing better documentation, examples, and a tutorial, at [https://libpfasst.github.io/LibPFASST](https://libpfasst.github.io/LibPFASST).
The documentation can be edited on the gh-pages branch of this repo.

## Requirements
- a Fortran compiler, for some external libraries like FFTW C/C++ compiler
- an MPI library

## Quickstart
- compile the library, see [https://libpfasst.github.io/LibPFASST/docs/build/html/compiling.html](https://libpfasst.github.io/LibPFASST/docs/build/html/compiling.html)
- look through the tutorials, see [https://libpfasst.github.io/LibPFASST/docs/build/html/tutorial.html](https://libpfasst.github.io/LibPFASST/docs/build/html/tutorial.html)

## License

Libpfasst Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) and Sebastian Goetschel.  All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a  non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.










