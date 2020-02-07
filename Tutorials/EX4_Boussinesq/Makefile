#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../../
LIBNPY ?= $(LIBPFASST)/libnpy

EXE = main.exe
BUILDDIR = build

include $(LIBPFASST)/Makefile.defaults
FSRC = main.f90 sweeper.f90 level.f90 hooks.f90 probin.f90 

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include  -I$(LIBNPY)/include 
LDFLAGS  += -L$(LIBNPY)/libnpy_fortran_mod -lnpy

all: $(EXE)

VPATHS = src 


include $(LIBPFASST)/Makefile.rules


# dependencies
#
include $(PWD)/.depend
main.exe : $(LIBPFASST)/lib/libpfasst.a 



