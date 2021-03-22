LIBPFASST ?= ../..
BUILDDIR = build

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = src/main.f90 src/sweeper.f90 src/level.f90 src/hooks.f90 src/probin.f90 src/encap.f90
CPPSRC = src/cencap.cpp src/cscalar.cpp

OBJ  = $(subst src, build,$(FSRC:.f90=.o))
OBJ += $(subst src, build,$(CPPSRC:.cpp=.o))

FFLAGS  += -I$(LIBPFASST)/include
LDFLAGS += -lstdc++

all: cpp_encap $(EXE)

VPATHS = src 

include $(LIBPFASST)/Makefile.rules

include .depend
main.exe : $(LIBPFASST)/lib/libpfasst.a 

cpp_encap:
	mkdir -p build
	g++ -c $(CPPSRC)
	mv cencap.o cscalar.o build/
