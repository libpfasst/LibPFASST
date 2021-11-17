#  Flags for compilers and linkers on local machine
#  This example is from GNU compilers on Linux
FC = mpifort
CC = mpicc
LD = $(FC)
AR=ar rcs

FFLAGS = -Ibuild -Jinclude -cpp -ffree-line-length-none
#  Use this flag on newer compilers with stricter bounds checking
#FFLAGS = -fallow-argument-mismatch

ifeq ($(DEBUG),TRUE)
FFLAGS += -fcheck=all -fbacktrace -g -ffpe-trap=invalid,zero,overflow -fbounds-check -fimplicit-none -ffree-line-length-none
else
FFLAGS += -O3 
endif