#
# Makefile for mpi-advection example.  This Makefile builds LIBPFASST
# inplace (instead of using the static library form of LIBPFASST) as
# well as the application specific code in the src directory.
#

LIBPFASST ?= ../..

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = feval.f90 \
       hooks.f90 \
       main.f90 \
       transfer.f90  

OBJ  = $(addprefix build/,$(FSRC:.f90=.o))
#CSRC =  pf_numpy.c
#OBJ  += $(addprefix build/,$(CSRC:.c=.o))

FFLAGS  += -I$(LIBPFASST)/include  -I$(FFTW3)/include -fbounds-check 
LDFLAGS +=  -L$(FFTW3)/lib  -L$(LIBPFASST)/lib -lpfasst -lfftw3

all: $(EXE)

VPATHS = src 

include $(LIBPFASST)/Makefile.rules

FFLAGS += -Wno-surprising
LDFLAGS += -g -lfftw3

build/%.o: %.f90
	@mkdir -p build
	@mkdir -p include
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)


#
# dependencies
#

build/feval.o: 
build/hooks.o: build/feval.o 
build/main.o:  build/hooks.o build/feval.o 
build/transfer.o:  build/feval.o  
