#
# Makefile for libpfasst.  This builds the static LIBPFASST library.
#

LIBPFASST ?= $(PWD)

include $(LIBPFASST)/Makefile.defaults

#
# libpfasst
#

build/libpfasst.a: $(OBJ)
	$(AR) build/libpfasst.a $(OBJ)

include $(LIBPFASST)/Makefile.rules
include $(LIBPFASST)/Makefile.external

EXTESTS = examples/mpi-advection examples/mpi-ndarray examples/fake-advection

tests: $(EXTESTS)
	@echo
	@echo '===> running nosetests'
	@$(NOSE)

$(EXTESTS):
	@echo "===>" $@
	$(MAKE) -C $@

.PHONY: clean tests $(EXTESTS)
