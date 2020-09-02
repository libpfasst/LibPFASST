##
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../../../LibPFASST
FFTW3     ?= $(LIBPFASST)/fftw3
LIBNPY ?= $(LIBPFASST)/libnpy

include $(LIBPFASST)/Makefile.defaults

BUILDDIR = build

USEFFTW=TRUE
NOFFTW=FALSE
EXE = main.exe

DEBUG=FALSE

VPATHS = src
FSRC   = src/probin.f90  src/sweeper.f90 src/hooks.f90  src/pf_optimization_1d.f90 src/level.f90 src/main.f90

OBJ = $(subst src, build, $(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include -I$(LIBNPY)/include  -g -fcheck=all,bounds -fbacktrace
#-I$(FFTW3)/include
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst -L$(LIBPFASST)/libnpy/libnpy_fortran_mod -lnpy
#-lblas -llapack
LDFLAGS += -L$(FFTW3)/lib -L$(FFTW3)/lib64 -lfftw3


all: $(EXE)

include $(LIBPFASST)/Makefile.rules

main.exe : $(LIBPFASST)/lib/libpfasst.a 
#
# dependencies
#

build/feval.o:         build/probin.o
build/pf_optimization_1d.o: build/sweeper.o build/probin.o
build/main.o:          build/feval.o build/level.o build/hooks.o build/probin.o $(LIBPFASST)/lib/libpfasst.a
