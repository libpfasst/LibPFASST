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

FFLAGS  += -g -I$(LIBPFASST)/include -I$(FFTW3)/include
LDFLAGS += -g -L$(LIBPFASST)/lib -lpfasst
LDFLAGS += -L$(FFTW3)/lib -lfftw3

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


