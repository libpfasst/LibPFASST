####################
# User-set variables
####################
# Hypre include directory
export HYPRE_DIR=./hypre/src/hypre
# C++ compiler (user must set this)
CPP_COMPILE = g++ -O3 -g
# LibPFASST directory
LIBPFASST ?= ../..

CPP_SRC_DIR=src/
export USE_HYPRE=TRUE

CPP_FILES = \
        $(CPP_SRC_DIR)hypre_struct.cpp \
        $(CPP_SRC_DIR)hypre_vector.cpp \
        $(CPP_SRC_DIR)hypre_solver.cpp \
        $(CPP_SRC_DIR)hypre_fortran.cpp

CPP_INCLUDE = \
        -I$(HYPRE_DIR)/include \
        -I./hypre/src/struct_ls

BUILDDIR = build

EXE = main.exe

# Must use Makefile.defaults included in this repo
include $(LIBPFASST)/Makefile.defaults

FSRC = src/comm.f90 src/level.f90 src/sweeper.f90 src/stepper.f90 src/hooks.f90 src/probin.f90 src/encap.f90 src/pfasst_hypre.f90 src/main.f90
CSRC = src/hypre_struct.cpp src/hypre_vector.cpp src/hypre_solver.cpp src/hypre_fortran.cpp

OBJ  = $(subst src, build,$(FSRC:.f90=.o))
OBJ += $(subst src, build,$(CSRC:.cpp=.o))

FFLAGS  += -I$(LIBPFASST)/include -g

all: hypre_cpp $(EXE)

hypre_cpp:
	mkdir -p build
	$(CPP_COMPILE) -c $(CPP_FILES) $(CPP_INCLUDE)
	mv $(subst src, .,$(CPP_FILES:.cpp=.o)) build/

VPATHS = src 

include $(LIBPFASST)/Makefile.rules

include .depend
main.exe : $(LIBPFASST)/lib/libpfasst.a
