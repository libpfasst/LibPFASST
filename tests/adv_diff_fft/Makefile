#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = feval.f90 \
       hooks.f90 \
       main.f90 \
       probin.f90 \
       fftpack5.1d.f90 
#       transfer.f90

OBJ = $(addprefix build/,$(FSRC:.f90=.o))

FFLAGS  += -g -I$(LIBPFASST)/include -I$(FFTW3)/include -fbounds-check
LDFLAGS += -g -L$(LIBPFASST)/lib -L$(FFTW3)/lib -lpfasst -lfftw3

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

build/feval.o:
build/probin.o:
build/hooks.o: build/feval.o
build/main.o:  build/hooks.o build/feval.o build/probin.o
build/transfer.o:  build/feval.o
