#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..

EXE = rossby_adjustment

include $(LIBPFASST)/Makefile.defs

FSRC = $(wildcard src/*.f90)
OBJ  = $(subst src,build,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/build -I$(CURDIR)/inc -I/home/daniel/Programs/hdf5/include
LDFLAGS += -L$(LIBPFASST)/build -lpfasst -L$(CURDIR)/lib -lspatialdiscretization -lserialio -L/home/daniel/Programs/hdf5/lib -lhdf5_fortran -lhdf5 -lz

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/encap.o
build/hooks.o: build/encap.o
build/main.o: build/feval.o build/hooks.o build/transfer.o
