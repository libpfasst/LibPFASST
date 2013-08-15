#
# Makefile for libpfasst.  This builds the static LIBPFASST library.
#

LIBPFASST ?= $(PWD)

include $(LIBPFASST)/Makefile.defaults

#
# libpfasst
#

$(shell $(PY) mk/version.py src/pf_version.f90)

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
