#
# Makefile for 2D Boussinesq Example
#

LIBPFASST ?= ../../
LIBNPY ?= $(LIBPFASST)/libnpy

EXE = main.exe
BUILDDIR = build

include $(LIBPFASST)/Makefile.defaults
FSRC = src/main.f90 src/sweeper.f90 src/level.f90 src/hooks.f90 src/probin.f90 

OBJ  = $(subst src, build,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include  -I$(LIBNPY)/include 
LDFLAGS  += -L$(LIBNPY)/libnpy_fortran_mod -lnpy

all: $(EXE)

VPATHS = src 


include $(LIBPFASST)/Makefile.rules


# dependencies
#
include $(PWD)/.depend
main.exe : $(LIBPFASST)/lib/libpfasst.a 



