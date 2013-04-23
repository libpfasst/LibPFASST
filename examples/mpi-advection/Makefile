#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(HOME)/opt

EXE = mpi-advection


include $(LIBPFASST)/Makefile.defs

FSRC = $(shell ls src/*.f90)
OBJ  = $(subst src,build,$(FSRC:.f90=.o))

FFLAGS  += -Ibuild -Jbuild -I$(FFTW3)/include -I$(LIBPFASST)/build 
LDFLAGS += -L$(LIBPFASST)/build -L$(FFTW3)/lib -lfftw3

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/encap.o
build/hooks.o: build/encap.o
build/main.o: build/feval.o build/hooks.o build/transfer.o
