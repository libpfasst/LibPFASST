include $(CURDIR)/../../../Makefile.defs

all: libinitmpi.a

libinitmpi.a: initmpi.o
	ar crs libinitmpi.a initmpi.o

initmpi.o: initmpi.f90
	$(FC) $(FFLAGS) -c initmpi.f90

initmpi.mod: initmpi.o

install: libinitmpi.a initmpi.mod
	cp libinitmpi.a $(CURDIR)/../lib
	cp *.mod $(CURDIR)/../inc

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
	rm -rf $(CURDIR)/../lib/libinitmpi.a
	rm -rf $(CURDIR)/../inc/initmpi.mod
