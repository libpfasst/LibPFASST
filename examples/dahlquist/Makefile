#
# Makefile for Dahlquist equation
#

LIBPFASST ?= ../..
LAPACKLIB = /usr/lib

EXE = main

include $(LIBPFASST)/Makefile.defaults

FSRC = probin.f90 \
       encap.f90 \
       feval.f90 \
       transfer.f90 \
       hooks.f90 \
       main.f90

OBJ  = $(addprefix build/,$(FSRC:.f90=.o))

FFLAGS += -I$(LIBPFASST)/build -fbounds-check -fbacktrace
LDFLAGS += -L$(LIBPFASST)/build -lpfasst 
LDFLAGS += -L$(LAPACKLIB) -llapack

all: $(EXE)

#
# f90 source files live in `src`
#

VPATHS += src

include $(LIBPFASST)/Makefile.rules

#
# Dependencies
#

build/feval.o: build/probin.o build/encap.o
built/transfer.o: build/probin.o build/encap.o build/feval.o
build/hooks.o: build/probin.o build/encap.o 
build/main.o: build/probin.o build/encap.o build/feval.o build/transfer.o build/hooks.o
