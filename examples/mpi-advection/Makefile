#
# Makefile for mpi-advection example.  This Makefile builds LIBPFASST
# inplace (instead of using the static library form of LIBPFASST) as
# well as the application specific code in the src directory.
#

LIBPFASST ?= ../..
FFTW3     ?= ../../fftw3

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

VPATHS += src
FSRC   += src/main.f90 src/feval.f90 src/transfer.f90 src/hooks.f90

FFLAGS  += -I$(FFTW3)/include 
LDFLAGS += -L$(FFTW3)/lib -lfftw3

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/pfasst.o
build/hooks.o: build/pfasst.o
build/main.o: build/feval.o build/hooks.o build/transfer.o
