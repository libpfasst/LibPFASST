#
# Makefile for Exponential RK examples
#

LIBPFASST ?= ../../..
DO_EXP=FALSE
DO_RK=TRUE

DIM=1
EXE = main.$(DIM)d.exe
all: $(EXE)
BUILDDIR = build

include $(LIBPFASST)/Makefile.defaults
FSRC = src/stepper_$(DIM)d.f90 

include ../Make.package
OBJ  = $(subst src, build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

FFLAGS  += -I$(LIBPFASST)/include 

all: $(EXE)

VPATHS = src:../src 


include $(LIBPFASST)/Makefile.rules


# dependencies
#

main.exe : $(LIBPFASST)/lib/libpfasst.a 
build/stepper_$(DIM)d.o : src/stepper_$(DIM)d.f90  build/probin.o   build/utils_$(DIM)d.o src/stepper_include.f90

