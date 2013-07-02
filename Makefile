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


# XXX: all of the examples should be buildable and runnable too
tests:
	@echo
	@echo '===> building tests'
	cd tests && make

	@echo
	@echo '===> running nosetests'
	@$(NOSE)


.PHONY: clean tests
