include $(CURDIR)/../../../Makefile.defs

FFLAGS += -I$(CURDIR)/../inc -I/home/daniel/Programs/hdf5/include -I$(CURDIR)

#LDFLAGS += -L$(CURDIR)/../lib -l lserial_io

all: libspatialdiscretization.a

# Rule for library file
libspatialdiscretization.a: timediscretization_parameter.o FVMParameters.o MPIParameter.o FiniteVolumes.o DistributedIO.o mpi_space_communication.o RHSFunctions.o spatialdiscretization.o
	ar crs libspatialdiscretization.a FVMParameters.o FiniteVolumes.o RHSFunctions.o MPIParameter.o mpi_space_communication.o DistributedIO.o spatialdiscretization.o timediscretization_parameter.o

# Rules for object files

timediscretization_parameter.o: $(CURDIR)/../lib/libserialio.a timediscretization_parameter.f90
	$(FC) $(FFLAGS) -c timediscretization_parameter.f90

FVMParameters.o: $(CURDIR)/../lib/libserialio.a FiniteVolumes.o FVMParameters.f90
	$(FC) $(FFLAGS) -c FVMParameters.f90 

MPIParameter.o: FVMParameters.o FiniteVolumes.o MPIParameter.f90
	$(FC) $(FFLAGS) -c MPIParameter.f90 

FiniteVolumes.o: WENO5_Boussinesq.f90
	$(FC) $(FFLAGS) -c WENO5_Boussinesq.f90 -o FiniteVolumes.o

DistributedIO.o: FVMParameters.o MPIParameter.o FiniteVolumes.o DistributedIO.f90
	$(FC) $(FFLAGS) -c DistributedIO.f90

mpi_space_communication.o: FVMParameters.o MPIParameter.o FiniteVolumes.o mpi_space_communication.f90
	$(FC) $(FFLAGS) -c mpi_space_communication.f90

RHSFunctions.o: FVMParameters.o FiniteVolumes.o MPIParameter.o DistributedIO.o mpi_space_communication.o RHSFunctions.f90
	$(FC) $(FFLAGS) -c RHSFunctions.f90 -o RHSFunctions.o

spatialdiscretization.o: FVMParameters.o MPIParameter.o FiniteVolumes.o DistributedIO.o mpi_space_communication.o RHSFunctions.o spatialdiscretization.f90
	$(FC) $(FFLAGS) -c spatialdiscretization.f90

spatialdiscretization.mod: spatialdiscretization.o

# Rule to update serial io library if need
$(CURDIR)/../lib/libserialio.a:
	cd $(CURDIR)/../libserialio; make
	cd $(CURDIR)/../libserialio; make install

install: libspatialdiscretization.a spatialdiscretization.mod
	cp libspatialdiscretization.a $(CURDIR)/../lib
	cp *.mod $(CURDIR)/../inc

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
	rm -rf $(CURDIR)/../lib/libspatialdiscretization.a
	rm -rf $(CURDIR)/../inc/spatialdiscretization.mod
