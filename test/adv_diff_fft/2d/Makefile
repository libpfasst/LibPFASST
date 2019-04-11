#
# Makefile for advection diffusion example.
#


LIBPFASST ?= ../../..
BUILDDIR = build

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = main.f90 feval.f90 hooks.f90 probin.f90 

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include 

all: $(EXE)

VPATHS = src 


include $(LIBPFASST)/Makefile.rules

#
# dependencies
#
include $(PWD)/.depend
main.exe : $(LIBPFASST)/lib/libpfasst.a 


