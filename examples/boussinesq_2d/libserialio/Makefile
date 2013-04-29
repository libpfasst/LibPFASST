include $(CURDIR)/../../../Makefile.defs

all: libserialio.a

libserialio.a: serialio.o
	ar rcs libserialio.a serialio.o

serialio.o: serialio.f90
	$(FC) $(FFLAGS) -I/home/daniel/Programs/hdf5/include -c serialio.f90

serialio.mod: serialio.o

install: libserialio.a serialio.mod
	cp libserialio.a $(CURDIR)/../lib
	cp *.mod $(CURDIR)/../inc

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
	rm -rf $(CURDIR)/../lib/libserialio.a
	rm -rf $(CURDIR)/../inc/serialio.mod
