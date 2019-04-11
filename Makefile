#
# Makefile for libpfasst.  This builds the static LIBPFASST library.
#

LIBPFASST ?= $(PWD)
BUILDDIR = build
include $(LIBPFASST)/Makefile.defaults

#
# libpfasst
#

build/libpfasst.a: $(OBJ)
	@mkdir -p lib
	$(AR) lib/libpfasst.a $(OBJ)

include $(LIBPFASST)/Makefile.rules
include $(LIBPFASST)/Makefile.external

.PHONY: clean depend
