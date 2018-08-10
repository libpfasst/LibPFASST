#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = src/feval.f90 \
       src/hooks.f90 \
       src/main.f90 \
       src/encap.f90 \
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
build/feval.o : src/feval.f90 ./build/probin.o  ./build/encap.o 
build/hooks.o : src/hooks.f90 ./build/feval.o 
build/main.o : src/main.f90 ./build/probin.o ./build/hooks.o ./build/feval.o  ./build/encap.o 
build/probin.o : src/probin.f90 
build/encap.o : src/encap.f90 
