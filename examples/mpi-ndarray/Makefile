#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3

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

build/solutions.o: build/pfasst.o build/probin.o
build/feval.o:     build/pfasst.o build/probin.o
build/hooks.o:     build/numpy.o build/solutions.o
build/main.o:      build/feval.o build/hooks.o build/transfer.o
