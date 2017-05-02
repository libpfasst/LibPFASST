#
# Applications
#

my_GCC_DIR = $(HOME)/opt/gcc/6.3.0
my_MPICH_DIR = $(my_GCC_DIR)/mpich/3.2
my_HDF5_DIR = $(my_MPICH_DIR)/hdf5/1.10.0-patch1
my_PETSC_DIR = $(my_MPICH_DIR)/petsc/3.7.4

#
# Compilers
#

CC = $(my_MPICH_DIR)/bin/mpicc
FC = $(my_MPICH_DIR)/bin/mpifort

#
# LIBPFASST
#

LIBPFASST_DIR ?= $(my_MPICH_DIR)/libpfasst
include $(LIBPFASST_DIR)/Makefile.defaults

#
# Run `make PETSC_DIR="$my_PETSC_DIR" getincludedirs` in PETSc build directory for "<petsc-fflags>" in the following
#

FFLAGS += <petsc-fflags> -I$(my_HDF5_DIR)/include

#
# Run `make PETSC_DIR="$my_PETSC_DIR" getlinklibs` in PETSc build directory for "<petsc-ldflags>" in the following
#

LDFLAGS += <petsc-ldflags> -L$(my_HDF5_DIR)/lib -lhdf5 -lhdf5_fortran
