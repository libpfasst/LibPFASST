#
# Makefile for libpfasst.
#

#
# config
#

FC     = mpif90
FFLAGS = -fPIC -Wall -g -pg -Ibuild -Jbuild

STATIC = ar rcs
SHARED = gfortran -shared -Wl,-soname,pfasst

FSRC = $(shell ls src/*.f90)
OBJ  = $(subst src,build,$(FSRC:.f90=.o))

MKVERBOSE =

#
# rules
#

all: libpfasst.so libpfasst.a

libpfasst.so: $(OBJ)
	$(SHARED) -o libpfasst.so $(OBJ)

libpfasst.a: $(OBJ)
	$(STATIC) libpfasst.a $(OBJ)

build/%.o: src/%.f90
	@mkdir -p build
ifdef MKVERBOSE
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)
else
	@echo FC $(notdir $<)
	@$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)
endif

.PHONY: clean

clean:
	rm -rf build libpfasst.*

#
# dependencies (generated by: python tools/f90_mod_deps.py --o-prefix build/ src/*.f90)
#

build/pfasst.o: build/pf_dtype.o build/pf_hooks.o build/pf_parallel.o build/pf_pfasst.o build/pf_version.o build/pf_implicit.o build/pf_explicit.o build/pf_imex.o
build/pf_dtype.o:
build/pf_explicit.o: build/pf_dtype.o build/pf_timer.o
build/pf_hooks.o: build/pf_dtype.o build/pf_timer.o
build/pf_imex.o: build/pf_dtype.o build/pf_timer.o
build/pf_implicit.o: build/pf_dtype.o build/pf_timer.o
build/pf_interpolate.o: build/pf_dtype.o build/pf_restrict.o build/pf_timer.o
build/pf_mpi.o: build/pf_dtype.o build/pf_timer.o
build/pf_parallel.o: build/pf_dtype.o build/pf_interpolate.o build/pf_restrict.o build/pf_utils.o build/pf_timer.o build/pf_hooks.o
build/pf_pfasst.o: build/pf_dtype.o build/pf_utils.o build/pf_version.o build/pf_quadrature.o build/pf_mpi.o
build/pf_quadrature.o: build/pf_dtype.o build/sdc_quadrature.o
build/pf_restrict.o: build/pf_dtype.o build/pf_utils.o build/pf_timer.o
build/pf_timer.o: build/pf_dtype.o
build/pf_utils.o: build/pf_dtype.o
build/sdc_poly.o:
build/sdc_quadrature.o: build/sdc_poly.o build/pf_dtype.o
