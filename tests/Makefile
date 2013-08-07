#
# Makefile for mpi-advection example.  This Makefile builds LIBPFASST
# inplace (instead of using the static library form of LIBPFASST) as
# well as the application specific code in the src directory.
#

LIBPFASST ?= ..
FFTW3     ?= ../fftw3

include $(LIBPFASST)/Makefile.defaults

FFLAGS  += -I$(FFTW3)/include 
LDFLAGS += -L$(FFTW3)/lib -lfftw3

VPATHS = $(LIBPFASST)/src src

all: t1.exe t2.exe t3.exe

t1.exe: src/spec_ad.f90 src/t1.f90 $(OBJ)
	$(FC) $(FFLAGS) $^ $(OUTPUT_OPTION) $(LDFLAGS)

t2.exe: src/spec_ad.f90 src/t2.f90 $(OBJ)
	$(FC) $(FFLAGS) $^ $(OUTPUT_OPTION) $(LDFLAGS)

t3.exe: src/spec_ad.f90 src/t3.f90 $(OBJ)
	$(FC) $(FFLAGS) $^ $(OUTPUT_OPTION) $(LDFLAGS)


include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/pfasst.o
build/hooks.o: build/pfasst.o
build/main.o: build/feval.o build/hooks.o build/transfer.o

