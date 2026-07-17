#
# Makefile for advection-diffusion example in 1D/2D using LibPFASST + AMReX.
#
DIM=1
LIBPFASST ?= /home/cfb7172/lib/LibPFASST
LIBNPY ?= $(LIBPFASST)/libnpy
BUILDDIR = BUILD

USE_AMREX = TRUE
USE_DIAG=FALSE
ifeq ($(USE_DIAG),TRUE)
EXE = main_$(DIM)d_diag.exe
else
EXE = main_$(DIM)d.exe
endif


$(info "DIM is $(DIM)" )
$(info "use AMReX is $(USE_AMREX)")
$(info "use Diagonal Sweeper is $(USE_DIAG)")

include $(LIBPFASST)/Makefile.defaults

ifeq ($(USE_DIAG),TRUE)
$(info Using diagonal sweeper)
FSRC = src/sweeper_diag_$(DIM)d.f90
else
$(info Using standard sweeper)
FSRC = src/sweeper_$(DIM)d.f90
endif

#  Include the standard 
FSRC +=  src/probin.f90  src/utils_$(DIM)d.f90 src/dim_$(DIM)d.f90
#  Choose the main and level to include
FSRC += src/main.f90 src/hooks.f90 src/level.f90 src/mpi_comm.f90

OBJ  = $(subst src, build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

FFLAGS  += -I$(LIBPFASST)/include  -I$(LIBNPY)/include 
#LDFLAGS  += -L$(LIBNPY)/libnpy_fortran_mod -lnpy

AMREX_HOME=$(LIBPFASST)/../amrex
LDFLAGS += -lgfortran -lpthread -lstdc++
# use flag below for debugging
#LDFLAGS += -fsanitize=address

all: $(EXE)

VPATHS = src

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#
include .depend

