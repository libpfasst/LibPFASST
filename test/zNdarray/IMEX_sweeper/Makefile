#
# Makefile for mpi-advection example.
#
DIM=1
LIBPFASST ?= ../../..

EXE = main.$(DIM)d.exe
all: $(EXE)
BUILDDIR = build

include $(LIBPFASST)/Makefile.defaults

FSRC = src/sweeper_$(DIM)d.f90  
include ../Make.package
OBJ  = $(subst src, build,$(FSRC:.f90=.o) $(CSRC:.c=.o))


FFLAGS  += -I$(LIBPFASST)/include 


VPATHS = src:../src


include $(LIBPFASST)/Makefile.rules




# dependencies
#

main.exe : $(LIBPFASST)/lib/libpfasst.a $(OBJ)
build/sweeper_$(DIM)d.o : src/sweeper_$(DIM)d.f90 ./build/probin.o ./build/utils_$(DIM)d.o $(LIBPFASST)/lib/libpfasst.a src/sweeper_include.f90


