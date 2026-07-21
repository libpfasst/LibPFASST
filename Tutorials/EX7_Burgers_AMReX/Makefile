#
# Makefile for advection-diffusion example in 1D/2D using LibPFASST + AMReX.
#
DIM=1
LIBPFASST ?= ../..7172/lib/LibPFASST
LIBNPY ?= $(LIBPFASST)/libnpy
BUILDDIR = BUILD

EXE = main_$(DIM)d.exe
USE_AMREX = TRUE

$(info "DIM is $(DIM)" )
$(info "use AMReX is $(USE_AMREX)")

include $(LIBPFASST)/Makefile.defaults

FSRC = src/sweeper_$(DIM)d.f90
#  Include the standard 
FSRC +=  src/probin.f90  src/utils_$(DIM)d.f90 src/dim_$(DIM)d.f90
#  Choose the main and level to include
FSRC += src/main.f90 src/hooks.f90 src/level.f90 src/mpi_comm.f90

OBJ  = $(subst src, build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

FFLAGS  += -I$(LIBPFASST)/include  -I$(LIBNPY)/include 
#LDFLAGS  += -L$(LIBNPY)/libnpy_fortran_mod -lnpy

AMREX_HOME=$(LIBPFASST)/../amrex
LDFLAGS += -lgfortran -lpthread -lstdc++

all: $(EXE)

VPATHS = src

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#
include .depend

