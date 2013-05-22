FC = /Users/ruprechd/Programs/MPICH_GCC46/bin/mpif90
FFLAGS = -cpp -fopenmp -Wtabs -ffree-line-length-none 
#FFLAGS += -fbounds-check -frange-check -fbacktrace -fvar-tracking -C -Wconversion -fimplicit-none -Warray-temporaries -Wcharacter-truncation -Wline-truncation -fcheck=all -g

all: libserialio.a

libserialio.a: serialio.o
	ar rcs libserialio.a serialio.o

serialio.o: serialio.f90
	$(FC) $(FFLAGS) -I/home/daniel/Programs/hdf5/include -c serialio.f90

serialio.mod: serialio.o

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a