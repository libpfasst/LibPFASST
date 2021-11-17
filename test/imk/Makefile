# Makefile for pfasst-nwchem example.
LIBPFASST = ../..
BUILDDIR = build

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = probin.f90 factory.f90 sweeper.f90 level.f90 main.f90 hooks.f90 utils.f90

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include 
LDFLAGS +=  -lopenblas -llapack

all: $(EXE)

VPATHS = src

include $(LIBPFASST)/Makefile.rules


include $(PWD)/.depend
main.exe : $(LIBPFASST)/lib/libpfasst.a 
