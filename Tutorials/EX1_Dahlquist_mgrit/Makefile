#
# Makefile for Dahlyquist example.
#

LIBPFASST ?= ../..
BUILDDIR = build

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = src/hooks.f90 src/stepper.f90 src/main.f90 src/level.f90 src/probin.f90 

OBJ  = $(subst src, build,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include

all: $(EXE)

VPATHS = src 

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#
include .depend
main.exe : $(LIBPFASST)/lib/libpfasst.a 
