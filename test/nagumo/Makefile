#
# Makefile for nagumo example.
#

LIBPFASST ?= ../..
BUILDDIR = build

EXE = main_split.exe

include $(LIBPFASST)/Makefile.defaults


FSRC   =  feval.f90 hooks.f90 solutions.f90 probin.f90 pf_optimization_flex.f90 main_split.f90
CSRC   = numpy.c

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o) $(CSRC:.c=.o))
FFLAGS  += -I$(LIBPFASST)/include
# -O -Wall -fcheck=all -g -fbacktrace 

all: $(EXE)

VPATHS = src

include $(LIBPFASST)/Makefile.rules
#
# dependencies
#
include $(PWD)/.depend
