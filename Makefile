#
# Makefile for libpfasst.  This builds the static LIBPFASST library.
#

LIBPFASST ?= $(PWD)

include $(LIBPFASST)/Makefile.defaults

#
# libpfasst
#

build/libpfasst.a: $(OBJ)
	@mkdir -p lib
	$(AR) lib/libpfasst.a $(OBJ)

include $(LIBPFASST)/Makefile.rules
include $(LIBPFASST)/Makefile.external

EXTESTS = examples/mpi-advection

tests: $(EXTESTS)
	@echo
	@echo '===> running nosetests'
	@$(NOSE)

$(EXTESTS):
	@echo "===>" $@
	$(MAKE) -C $@

.PHONY: clean tests $(EXTESTS)
