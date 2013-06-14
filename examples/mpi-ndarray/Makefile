#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = $(wildcard src/*.f90)
CSRC = $(wildcard src/*.c)
OBJ  = $(subst src,build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

FFLAGS  += -I$(LIBPFASST)/build -I$(FFTW3)/include 
LDFLAGS += -L$(LIBPFASST)/build -lpfasst -L$(FFTW3)/lib -lfftw3

VPATHS += src

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/solutions.o: build/encap.o build/probin.o
build/feval.o:     build/encap.o build/probin.o
build/hooks.o:     build/numpy.o build/solutions.o
build/main.o:      build/feval.o build/hooks.o build/transfer.o
