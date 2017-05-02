#
# Helper for this Makefile
#

HOSTNAME = $(shell hostname)

#
# Machines (see above PETSc `make` commands)
#

# Hostname
ifeq ($(findstring <hostname>,$(HOSTNAME)),<hostname>)
include host/<hostname>.mk
# Abort
else
$(error Your hostname $(HOSTNAME) is unknown)
endif

#
# Includes
#

FFLAGS += -I$(LIBPFASST_DIR)/include -cpp -fPIC -fopenmp -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fcheck=all -g3 -ggdb

#
# Linker
#

LDFLAGS += -L$(LIBPFASST_DIR)/lib -lpfasst -dynamic

#
# Sources
#

FSRC = mod_mpi.f90 mod_pf.f90 mod_deriv.f90 mod_interp.f90 mod_ctx.f90 mod_encap.f90 mod_transfer.f90 mod_prob.f90 mod_hook.f90 mod_petsc.f90 mod_feval.f90 main.f90
OBJ = $(addprefix build/,$(FSRC:.f90=.o))

#
# Executable
#

EXE = main
all: $(EXE)

#
# Source
#

VPATHS += src
include $(LIBPFASST_DIR)/Makefile.rules

#
# Dependencies
#

build/main.o: build/mod_mpi.o build/mod_pf.o build/mod_deriv.o build/mod_interp.o build/mod_ctx.o build/mod_encap.o build/mod_transfer.o build/mod_prob.o build/mod_hook.o build/mod_petsc.o build/mod_feval.o
