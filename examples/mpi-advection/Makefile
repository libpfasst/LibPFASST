#
# Makefile for mpi-advection example.  This Makefile builds LIBPFASST
# inplace (instead of using the static library form of LIBPFASST) as
# well as the application specific code in the src directory.
#

LIBPFASST ?= ../..

EXE = main.exe

FC = mpif90
CC = mpicc
AR = ar rcs
PY = python
NOSE = nosetests

FFLAGS = -g -Ibuild -Jbuild -cpp -fbacktrace

include $(wildcard $(LIBPFASST)/*.defs)

VPATHS = $(LIBPFASST)/src

FSRC = src/pfasst.f90 \
       src/pf_dtype.f90 \
       src/pf_explicit.f90 \
       src/pf_hooks.f90 \
       src/pf_imex.f90 \
       src/pf_imexQ.f90 \
       src/pf_implicit.f90 \
       src/pf_interpolate.f90 \
       src/pf_logger.f90 \
       src/pf_ndarray.f90 \
       src/pf_options.f90 \
       src/pf_parallel.f90 \
       src/pf_pfasst.f90 \
       src/pf_quadrature.f90 \
       src/pf_restrict.f90 \
       src/pf_timer.f90 \
       src/pf_utils.f90 \
       src/sdc_poly.f90 \
       src/sdc_quadrature.f90

CSRC = src/pf_numpy.c

ifdef NOMPI
  FFLAGS += -DNOMPI
else
  FSRC += src/pf_mpi.f90
endif

OBJ = $(subst src,build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

VPATHS += src
FSRC   += src/main.f90 src/feval.f90 src/transfer.f90 src/hooks.f90

FFLAGS += -Wno-surprising
LDFLAGS += -g -lfftw3

all: $(EXE)

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#

build/feval.o: build/pfasst.o
build/hooks.o: build/pfasst.o
build/main.o: build/feval.o build/hooks.o build/transfer.o
