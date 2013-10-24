#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

VPATHS += src
FSRC   += src/main.f90 src/feval.f90 src/transfer.f90 src/hooks.f90 src/encap.f90 src/fftwpp.f90
CSRC   += src/numpy.c
CXXSRC += src/cfftw++.cpp src/fftw++.cpp src/convolution.cpp 

FFLAGS  += -I$(FFTW3)/include 
LDFLAGS += -L$(FFTW3)/lib -lfftw3

build/%.o: %.cpp
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< $(OUTPUT_OPTION)

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/solutions.o: build/pfasst.o 
build/feval.o:     build/pfasst.o build/encap.o build/fftwpp.o 
build/hooks.o:     build/numpy.o build/solutions.o
build/main.o:      build/feval.o build/hooks.o build/transfer.o
