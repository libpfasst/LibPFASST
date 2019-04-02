#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../../..


FFTW=TRUE
EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = src/sweeper.f90 \
       src/hooks.f90 \
       src/main.f90 \
       src/probin.f90 \
       src/phi_mod.f90

OBJ = $(subst src,build,$(FSRC:.f90=.o))
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst
FFLAGS  += -I$(LIBPFASST)/include 


all: $(EXE)

VPATHS = src 


include $(LIBPFASST)/Makefile.rules

build/%.o: %.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< $(OUTPUT_OPTION)


build/%.o: %.f90
	@mkdir -p build
	@mkdir -p include
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)


#
# dependencies
#
./build/sweeper.o : src/sweeper.f90 ./build/probin.o ./build/phi_mod.o
./build/hooks.o : src/hooks.f90 ./build/sweeper.o
./build/main.o : src/main.f90 ./build/probin.o ./build/hooks.o ./build/sweeper.o ./build/phi_mod.o
./build/probin.o : src/probin.f90
./build/phi_mod.o : src/phi_mod.f90
main.exe : $(LIBPFASST)/lib/libpfasst.a 


