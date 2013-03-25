#
# Makefile for libpfasst.
#

# generate src/pf_version.f90
$(shell python mk/version.py src/pf_version.f90)

# configure
include Makefile.defs

FFLAGS += -Ibuild -Jbuild

#
# rules
#

FSRC = $(shell ls src/*.f90)
OBJ  = $(subst src,build,$(FSRC:.f90=.o))

all: build/libpfasst.a

build/libpfasst.a: $(OBJ)
ifdef MKVERbOSE
	$(AR) build/libpfasst.a $(OBJ)
else
	@echo AR build/libpfasst.a
	@$(AR) build/libpfasst.a $(OBJ)
endif

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
# dependencies
#

build/pf_utils.o:       build/pf_dtype.o
build/pf_timer.o:       build/pf_dtype.o
build/pf_explicit.o:    build/pf_timer.o
build/pf_hooks.o:       build/pf_timer.o
build/pf_imex.o:        build/pf_timer.o
build/pf_implicit.o:    build/pf_timer.o
build/pf_mpi.o:         build/pf_timer.o
build/pf_restrict.o:    build/pf_utils.o build/pf_timer.o
build/pf_interpolate.o: build/pf_restrict.o
build/pf_parallel.o:    build/pf_interpolate.o build/pf_hooks.o
build/sdc_quadrature.o: build/pf_dtype.o build/sdc_poly.o
build/pf_quadrature.o:  build/sdc_quadrature.o
build/pf_pfasst.o:      build/pf_utils.o build/pf_quadrature.o

build/pfasst.o:         build/pf_parallel.o build/pf_pfasst.o build/pf_implicit.o build/pf_explicit.o build/pf_imex.o build/pf_mpi.o build/pf_version.o

