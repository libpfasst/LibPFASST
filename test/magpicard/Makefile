# Makefile for pfasst-nwchem example.
LIBPFASST = ../..
BUILDDIR = build

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = probin.f90  sweeper.f90 main.f90 hooks.f90  utils.f90

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst -lopenblas -llapack

VPATHS = src
include $(LIBPFASST)/Makefile.rules
all: $(EXE)

$(BUILDDIR)/sweeper.o: $(BUILDDIR)/utils.o 
$(BUILDDIR)/utils.o: 
$(BUILDDIR)/hooks.o: $(BUILDDIR)/sweeper.o $(FSRC)
$(BUILDDIR)/main.o: $(LIBPFASST)/lib/libpfasst.a $(BUILDDIR)/sweeper.o $(BUILDDIR)/probin.o $(BUILDDIR)/hooks.o
