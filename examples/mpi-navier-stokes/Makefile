#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3
HDF5      ?= $(LIBPFASST)/hdf5

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

VPATHS += src
FSRC   += src/main.f90 src/feval.f90 src/transfer.f90 src/hooks.f90 src/encap.f90 src/fftwpp.f90 src/initial.f90 src/probin.f90
CSRC   += src/numpy.c
CXXSRC += src/cfftw++.cpp src/fftw++.cpp src/convolution.cpp

FFLAGS  += -O3 -I$(FFTW3)/include -I$(HDF5)/include -fopenmp
CFLAGS  += -I$(FFTW3)/include -fopenmp
LDFLAGS += -fopenmp -lstdc++ -ldl -lz -L$(FFTW3)/lib -lfftw3 -lfftw3_omp -L$(HDF5)/lib -lhdf5_fortran -lhdf5

vpath %.cpp $(VPATHS)

build/%.o: %.cpp
	@mkdir -p build
	$(CXX) $(CFLAGS) -c $< $(OUTPUT_OPTION)

all: $(EXE)

OBJ = $(subst src,build,$(FSRC:.f90=.o) $(CSRC:.c=.o) $(CXXSRC:.cpp=.o))

include $(LIBPFASST)/Makefile.rules

test_operators: test_operators.f90 $(OBJ)
	mpif90 -Ibuild -Jbuild -o test_operators test_operators.f90 $(filter-out build/main.o,$(OBJ)) $(LDFLAGS)

#
# dependencies
#

build/feval.o:     build/pfasst.o build/encap.o build/fftwpp.o
build/hooks.o:     build/numpy.o build/initial.o build/transfer.o
build/transfer.o:  build/feval.o
build/main.o:      build/feval.o build/hooks.o build/transfer.o build/initial.o build/probin.o
