#
# Makefile for nagumo example.
#

LIBPFASST ?= ../..

EXE = main_split.exe

include $(LIBPFASST)/Makefile.defaults

VPATHS += src
FSRC   =  src/feval.f90 src/hooks.f90 src/solutions.f90 src/probin.f90 src/pf_optimization.f90 src/main_split.f90 src/level.f90
CSRC   = src/numpy.c
FFLAGS  += -I$(LIBPFASST)/include 
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst 

all: $(EXE)


include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

./build/feval.o : src/feval.f90 ./build/solutions.o ./build/probin.o 
./build/hooks.o : src/hooks.f90 ./build/solutions.o ./build/probin.o 
./build/solutions.o : src/solutions.f90 ./build/probin.o 
./build/probin.o : src/probin.f90 
./build/pf_optimization.o : src/pf_optimization.f90 ./build/solutions.o ./build/probin.o ./build/feval.o ./build/level.o 
./build/main_split.o : src/main_split.f90 ./build/solutions.o ./build/probin.o ./build/hooks.o ./build/feval.o ./build/level.o ./build/pf_optimization.o $(LIBPFASST)/lib/libpfasst.a
./build/level.o : src/level.f90 ./build/feval.o ./build/solutions.o ./build/probin.o 
