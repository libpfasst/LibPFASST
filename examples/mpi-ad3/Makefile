#
# Makefile for mpi-advection example
#

LIBPFASST ?= ../..
include $(LIBPFASST)/Makefile.defaults

FFTW3 ?= /home/akreie/opt/fftw/3.3.4

FFLAGS += -g -I$(LIBPFASST)/include -I$(FFTW3)/include
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst -L$(FFTW3)/lib -lfftw3

FSRC = main.f90 feval.f90 transfer.f90 hooks.f90 solutions.f90 probin.f90 multigrid.f90
OBJ = $(addprefix build/,$(FSRC:.f90=.o))

EXE = main.exe
all: $(EXE)

VPATHS += src
include $(LIBPFASST)/Makefile.rules

#
# Dependencies
#

build/solutions.o: build/probin.o
build/feval.o: build/probin.o build/multigrid.o
build/hooks.o: build/numpy.o build/solutions.o
build/main.o: build/feval.o build/hooks.o build/transfer.o build/probin.o
