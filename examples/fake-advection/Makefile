#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = fake-advection

include $(LIBPFASST)/Makefile.defaults

FSRC = $(wildcard src/*.f90)
OBJ  = $(subst src,build,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/build -I$(FFTW3)/include 
LDFLAGS += -L$(LIBPFASST)/build -lpfasst -L$(FFTW3)/lib -lfftw3

VPATHS += src

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/encap.o
build/hooks.o: build/encap.o
build/main.o: build/feval.o build/hooks.o build/transfer.o
