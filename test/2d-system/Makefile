#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = src/main.f90 \
       src/feval.f90 \
       src/hooks.f90 \
       src/probin.f90 

OBJ = $(subst src,build,$(FSRC:.f90=.o))

FFLAGS  +=  -I$(LIBPFASST)/include 
LDFLAGS +=  -L$(LIBPFASST)/lib -lpfasst

all: $(EXE)

VPATHS = src

include $(LIBPFASST)/Makefile.rules

build/%.o: %.f90
	@mkdir -p build
	@mkdir -p include
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)


#
# dependencies
#
main.exe : $(LIBPFASST)/lib/libpfasst.a
build/main.o : src/main.f90 ./build/probin.o ./build/hooks.o ./build/feval.o 
build/feval.o : src/feval.f90 ./build/probin.o 
build/hooks.o : src/hooks.f90 ./build/feval.o 
build/probin.o : src/probin.f90
