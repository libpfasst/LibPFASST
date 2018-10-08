#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../../../libpfasst
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = main_split.exe

include $(LIBPFASST)/Makefile.defaults

VPATHS += src
FSRC   =  src/feval.f90 src/hooks.f90 src/solutions.f90 src/probin.f90 src/pf_optimization.f90 src/main_split.f90
CSRC   = src/numpy.c
FFLAGS  += -I$(FFTW3)/include -I$(LIBPFASST)/include 
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst 

all: $(EXE)


include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/solutions.o:     build/probin.o
build/feval.o:         build/probin.o build/solutions.o
build/pf_optimization.o: build/solutions.o build/feval.o build/probin.o
build/hooks.o:         build/solutions.o
build/main_split.o:    build/feval.o build/pf_optimization.o build/hooks.o build/probin.o build/numpy.o $(LIBPFASST)/lib/libpfasst.a
