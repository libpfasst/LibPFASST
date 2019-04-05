#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = src/main.f90 \
       src/sweeper.f90 \
       src/level.f90 \
       src/probin.f90 


OBJ = $(subst src,build,$(FSRC:.f90=.o))
FFLAGS  += -g -I$(LIBPFASST)/include
LDFLAGS += -g -L$(LIBPFASST)/lib -lpfasst


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
main.exe : $(LIBPFASST)/lib/libpfasst.a 
build/sweeper.o : src/sweeper.f90 ./build/probin.o 
build/level.o : src/level.f90 ./build/probin.o 
build/main.o : src/main.f90 ./build/probin.o  ./build/sweeper.o  ./build/level.o 
build/probin.o : src/probin.f90 
