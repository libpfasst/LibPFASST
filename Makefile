#
# Makefile for libpfasst.
#

LIBPFASST ?= $(PWD)

include $(LIBPFASST)/Makefile.defaults

#
# libpfasst
#

$(shell $(PY) mk/version.py src/pf_version.f90)

VPATHS = $(LIBPFASST)/src

FSRC = $(wildcard src/*.f90)
CSRC = $(wildcard ls src/*.c)
OBJ  = $(subst src,build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

build/libpfasst.a: $(OBJ)
	$(AR) build/libpfasst.a $(OBJ)

include $(LIBPFASST)/Makefile.rules

rossby_adjustment: examples/rossby_adjustment/lib/libspatialdiscretization.a examples/rossby_adjustment/*.f90
		$(FC) $(FFLAGS) -I$(CURDIR)/examples/rossby_adjustment/inc -o $@ $^ -Lbuild -lpfasst -L$(CURDIR)/examples/rossby_adjustment/lib -lspatialdiscretization -lserialio $(LDFLAGS)

boussinesq_2d: examples/boussinesq_2d/lib/libspatialdiscretization.a examples/boussinesq_2d/*.f90
		$(FC) $(FFLAGS) -I$(CURDIR)/examples/boussinesq_2d/inc -o $@ $^ -Lbuild -lpfasst -L$(CURDIR)/examples/boussinesq_2d/lib -lspatialdiscretization -lserialio $(LDFLAGS)

.PHONY: clean

#
# fftw3 (for the examples)
# 

FFTW3 = fftw-3.3.3
FFTW3_CONFIG = --enable-openmp --enable-sse2 --enable-threads # --enable-shared 

fftw3: $(FFTW3).tar.gz
	rm -rf $(FFTW3)
	tar xzf $(FFTW3).tar.gz
	cd $(FFTW3) && ./configure --prefix=$(PWD)/fftw3 $(FFTW3_CONFIG)
	cd $(FFTW3) && make install

$(FFTW3).tar.gz:
	wget http://www.fftw.org/$(FFTW3).tar.gz

#
# gcc-4.7
#

GCC    = gcc-4.7
GCCDL  = $(GCC).tar.xz xz.tar.gz gcc-infrastructure.tar.xz
GCCTMP = $(LIBPFASST)/gcc.tmp
XZCAT  = $(GCCTMP)/usr/bin/xzcat

gcc47: $(GCCDL)
	rm -rf $(GCCTMP) $(GCC) makefile.gcc47.defs
	mkdir $(GCCTMP)
	cd $(GCCTMP) && tar zxf $(LIBPFASST)/xz.tar.gz
	cd $(GCCTMP) && $(XZCAT) $(LIBPFASST)/$(GCC).tar.xz | tar x
	cd $(GCCTMP)/$(GCC) && $(XZCAT) $(LIBPFASST)/gcc-infrastructure.tar.xz | tar x
	mv $(GCCTMP)/$(GCC) $(LIBPFASST)
	rm -rf $(GCCTMP)
	ln -s $(LIBPFASST)/mk/makefile.gcc47.defs
	@echo
	@echo GCC 4.7 has been installed in 
	@echo 
	@echo   $(LIBPFASST)/$(GCC)
	@echo
	@echo and a symlink to makefile.gcc47.defs has been created.  These makefile
	@echo definitions should work either of the of OpenMPI or MPICH2 mpif90 
	@echo compiler wrappers.
	@echo
	@echo You may need to set your LD_LIBRARAY_PATH to include
	@echo
	@echo   $(LIBPFASST)/$(GCC)/lib:$(LIBPFASST)/$(GCC)/lib64
	@echo

$(GCC).tar.xz:
	wget http://gfortran.com/download/x86_64/snapshots/$(GCC).tar.xz

xz.tar.gz:
	wget http://gfortran.com/download/x86_64/xz.tar.gz

gcc-infrastructure.tar.xz:
	wget http://gfortran.com/download/x86_64/gcc-infrastructure.tar.xz


#
# dependencies
#

build/pf_utils.o:       build/pf_dtype.o
build/pf_timer.o:       build/pf_dtype.o
build/pf_cycle.o:       build/pf_dtype.o
build/pf_explicit.o:    build/pf_timer.o
build/pf_hooks.o:       build/pf_timer.o
build/pf_imex.o:        build/pf_timer.o
build/pf_implicit.o:    build/pf_timer.o
build/pf_mpi.o:         build/pf_timer.o
build/pf_pthreads.o:    build/pf_timer.o
build/pf_logger.o:      build/pf_hooks.o
build/pf_restrict.o:    build/pf_utils.o build/pf_timer.o build/pf_hooks.o
build/pf_interpolate.o: build/pf_restrict.o build/pf_hooks.o
build/pf_parallel.o:    build/pf_interpolate.o build/pf_hooks.o build/pf_cycle.o
build/sdc_quadrature.o: build/pf_dtype.o build/sdc_poly.o
build/pf_quadrature.o:  build/sdc_quadrature.o
build/pf_pfasst.o:      build/pf_utils.o build/pf_quadrature.o

build/pfasst.o:         build/pf_parallel.o build/pf_pfasst.o \
                        build/pf_implicit.o build/pf_explicit.o build/pf_imex.o \
                        build/pf_mpi.o build/pf_pthreads.o build/pf_cpthreads.o \
                        build/pf_version.o build/pf_logger.o build/pf_cycle.o
