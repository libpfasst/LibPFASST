#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../../..

EXE = main.exe
BUILDDIR = build

include $(LIBPFASST)/Makefile.defaults
FSRC = main.f90 sweeper.f90 level.f90 hooks.f90 probin.f90 phi_mod.f90

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include 

all: $(EXE)

VPATHS = src 


include $(LIBPFASST)/Makefile.rules


# dependencies
#
include $(PWD)/.depend
main.exe : $(LIBPFASST)/lib/libpfasst.a 



