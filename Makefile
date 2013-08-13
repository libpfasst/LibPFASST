#
# Makefile for libpfasst.  This builds the static LIBPFASST library.
#

LIBPFASST ?= $(PWD)

include $(LIBPFASST)/Makefile.defaults

#
# libpfasst
#

$(shell $(PY) mk/version.py src/pf_version.f90)

VPATHS = $(LIBPFASST)/src

FSRC = $(wildcard src/*.f90)
CSRC = $(wildcard src/*.c)
OBJ  = $(subst src,build,$(FSRC:.f90=.o) $(CSRC:.c=.o))

build/libpfasst.a: $(OBJ)
	$(AR) build/libpfasst.a $(OBJ)

include $(LIBPFASST)/Makefile.rules
include $(LIBPFASST)/Makefile.external


EXTESTS = examples/mpi-advection examples/mpi-ndarray

tests: $(EXTESTS)
	@echo
	@echo '===> running nosetests'
	@$(NOSE)

$(EXTESTS):
	@echo "===>" $@
	$(MAKE) -C $@

.PHONY: clean tests $(EXTESTS)
