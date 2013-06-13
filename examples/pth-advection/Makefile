#
# Makefile for pth-advection example.
#

LIBPFASST ?= ../..
FFTW3     ?= $(LIBPFASST)/fftw3

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = $(wildcard src/*.f90)
OBJ  = $(subst src,build,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/build -I$(FFTW3)/include 
LDFLAGS += -L$(LIBPFASST)/build -lpfasst -L$(FFTW3)/lib -lfftw3

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/encap.o
build/hooks.o: build/encap.o
build/main.o: build/feval.o build/hooks.o build/transfer.o build/pmain.o
