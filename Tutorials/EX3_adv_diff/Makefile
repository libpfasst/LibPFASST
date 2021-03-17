#
# Makefile for advection diffusion example
#

LIBPFASST ?= ../..
BUILDDIR = build

EXE = main.exe
USE_FFT=TRUE
include $(LIBPFASST)/Makefile.defaults

FSRC = src/main.f90 src/level.f90 src/sweeper.f90 src/hooks.f90 src/probin.f90 

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
