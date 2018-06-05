#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../../..
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = main_split.exe

include $(LIBPFASST)/Makefile.defaults

VPATHS += src
FSRC   =  src/feval.f90 src/hooks.f90 src/solutions.f90 src/probin.f90 src/pf_optimization_flex.f90 src/main_split.f90

FFLAGS  += -I$(FFTW3)/include -I$(LIBPFASST)/include 
#-g -fcheck=all,bounds
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst 
#-lblas -llapack
LDFLAGS += -L$(FFTW3)/lib -lfftw3

all: $(EXE)


include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/solutions.o:     build/probin.o
build/feval.o:         build/probin.o build/solutions.o
build/pf_optimization_flex.o: build/solutions.o build/feval.o build/probin.o
build/hooks.o:         build/numpy.o build/solutions.o
build/main.o:          build/feval.o build/hooks.o build/probin.o $(LIBPFASST)/lib/libpfasst.a
build/main_split.o:    build/feval.o build/pf_optimization_flex.o build/hooks.o build/probin.o $(LIBPFASST)/lib/libpfasst.a
build/main_heat.o:    build/feval.o build/pf_optimization_flex.o build/hooks.o build/probin.o $(LIBPFASST)/lib/libpfasst.a
