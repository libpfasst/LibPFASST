# Makefile for pfasst-nwchem example.
LIBPFASST = ../..
BUILDDIR = build

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

FSRC = probin.f90 factory.f90 sweeper.f90 level.f90 main.f90 hooks.f90 utils.f90

OBJ  = $(addprefix $(BUILDDIR)/,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include 
LDFLAGS += -L$(LIBPFASST)/lib -lpfasst -lopenblas -llapack

VPATHS = src
include $(LIBPFASST)/Makefile.rules
all: $(EXE)

$(BUILDDIR)/sweeper.o: src/sweeper.f90 $(BUILDDIR)/factory.o $(BUILDDIR)/utils.o
$(BUILDDIR)/level.o: src/level.f90 
$(BUILDDIR)/utils.o: src/utils.f90 $(BUILDDIR)/factory.o $(BUILDDIR)/probin.o
$(BUILDDIR)/hooks.o: src/hooks.f90 $(BUILDDIR)/factory.o $(BUILDDIR)/sweeper.o $(BUILDDIR)/utils.o $(FSRC)
$(BUILDDIR)/main.o:  src/main.f90 $(BUILDDIR)/factory.o  $(BUILDDIR)/sweeper.o $(BUILDDIR)/probin.o $(BUILDDIR)/hooks.o $(BUILDDIR)/level.o
main.exe : $(LIBPFASST)/lib/libpfasst.a 
