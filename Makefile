#
# Makefile for libpfasst.
#

include Makefile.defs

FFLAGS += -Ibuild -Jbuild

#
# rules
#

$(shell python mk/version.py src/pf_version.f90)

FSRC = $(shell ls src/*.f90)
CSRC = $(shell ls src/*.c)
OBJ  = $(subst src,build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

all: build/libpfasst.a

build/libpfasst.a: $(OBJ)
	$(AR) build/libpfasst.a $(OBJ)

build/%.o: src/%.f90
	@mkdir -p build
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< $(OUTPUT_OPTION)

advection: examples/advection/*.f90
	$(FC) $(FFLAGS) -I/home/memmett/opt/include -o $@ $^ -Lbuild -lpfasst -L/home/memmett/opt/lib -lfftw3


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
build/pf_pthreads.o:    build/pf_timer.o
build/pf_logger.o:      build/pf_hooks.o
build/pf_restrict.o:    build/pf_utils.o build/pf_timer.o
build/pf_interpolate.o: build/pf_restrict.o
build/pf_parallel.o:    build/pf_interpolate.o build/pf_hooks.o
build/sdc_quadrature.o: build/pf_dtype.o build/sdc_poly.o
build/pf_quadrature.o:  build/sdc_quadrature.o
build/pf_pfasst.o:      build/pf_utils.o build/pf_quadrature.o

build/pfasst.o:         build/pf_parallel.o build/pf_pfasst.o \
                        build/pf_implicit.o build/pf_explicit.o build/pf_imex.o \
                        build/pf_mpi.o build/pf_pthreads.o build/pf_cpthreads.o \
                        build/pf_version.o build/pf_logger.o
