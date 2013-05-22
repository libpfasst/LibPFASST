#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..

EXE = rossby_adjustment

include $(LIBPFASST)/Makefile.defaults

FSRC = $(wildcard src/*.f90)
OBJ  = $(subst src,build,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/build -I$(CURDIR)/inc -I/Users/ruprechd/Programs/MPICH_GCC46/include -I$(CURDIR)/libserialio -I$(CURDIR)/libspatialdiscretization
LDFLAGS += -L$(LIBPFASST)/build -lpfasst -L$(CURDIR)/libspatialdiscretization -lspatialdiscretization -L$(CURDIR)/libserialio -lserialio -L/Users/ruprechd/Programs/MPICH_GCC46/lib -lhdf5_fortran -lhdf5 -lz

VPATHS += src

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#
build/feval.o: build/encap.o
build/hooks.o: build/encap.o
build/main.o: build/feval.o build/hooks.o build/transfer.o
