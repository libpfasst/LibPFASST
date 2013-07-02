#
# Makefile for mpi-advection example.  This Makefile builds LIBPFASST
# inplace (instead of using the static library form of LIBPFASST) as
# well as the application specific code in the src directory.
#

LIBPFASST ?= ../..
FFTW3     ?= ../../fftw3

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = $(wildcard src/*.f90 $(LIBPFASST)/src/*.f90)
CSRC = $(wildcard src/*.c $(LIBPFASST)/src/*.c)
OBJ  = $(addprefix build/, $(notdir $(FSRC:.f90=.o) $(CSRC:.c=.o)))

FFLAGS  += -I$(FFTW3)/include 
LDFLAGS += -L$(FFTW3)/lib -lfftw3

VPATHS = $(LIBPFASST)/src src

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/pfasst.o
build/hooks.o: build/pfasst.o
build/main.o: build/feval.o build/hooks.o build/transfer.o
