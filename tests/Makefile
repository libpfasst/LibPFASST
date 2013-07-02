#
# Makefile for mpi-advection example.  This Makefile builds LIBPFASST
# inplace (instead of using the static library form of LIBPFASST) as
# well as the application specific code in the src directory.
#

LIBPFASST ?= ..
FFTW3     ?= ../fftw3

all: t1.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = $(wildcard $(LIBPFASST)/src/*.f90)
CSRC = $(wildcard $(LIBPFASST)/src/*.c)
OBJ  = $(addprefix build/, $(notdir $(FSRC:.f90=.o) $(CSRC:.c=.o)))

FFLAGS  += -I$(FFTW3)/include 
LDFLAGS += -L$(FFTW3)/lib -lfftw3

VPATHS = $(LIBPFASST)/src src

t1.exe: src/spec_ad.f90 src/t1.f90 $(OBJ)
	$(FC) $(FFLAGS) $^ $(OUTPUT_OPTION) $(LDFLAGS)





include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/pfasst.o
build/hooks.o: build/pfasst.o
build/main.o: build/feval.o build/hooks.o build/transfer.o

