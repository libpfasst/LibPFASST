#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..
BUILDDIR = build

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = main.f90  sweeper.f90 level.f90   probin.f90 

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include

all: $(EXE)

VPATHS = src 

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#
main.exe : $(LIBPFASST)/lib/libpfasst.a 
$(BUILDDIR)/main.o :  $(BUILDDIR)/probin.o  $(BUILDDIR)/sweeper.o  $(BUILDDIR)/level.o 
$(BUILDDIR)/sweeper.o :  $(BUILDDIR)/probin.o 
$(BUILDDIR)/level.o : 
$(BUILDDIR)/probin.o : 
