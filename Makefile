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

echo:
	@echo $(OBJ)

# XXX: all of the examples should be buildable and runnable too
tests:
	@echo
	@echo '===> building tests'
	cd tests && make

	@echo
	@echo '===> running nosetests'
	@$(NOSE)


.PHONY: clean tests echo

