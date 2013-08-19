#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

VPATHS += src
FSRC   += src/main.f90 src/feval.f90 src/transfer.f90 src/hooks.f90 src/solutions.f90 src/probin.f90
CSRC   += src/numpy.c

FFLAGS  += -I$(FFTW3)/include 
LDFLAGS += -L$(FFTW3)/lib -lfftw3

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/solutions.o: build/pfasst.o build/probin.o
build/feval.o:     build/pfasst.o build/probin.o
build/hooks.o:     build/numpy.o build/solutions.o
build/main.o:      build/feval.o build/hooks.o build/transfer.o
