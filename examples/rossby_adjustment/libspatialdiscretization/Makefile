FC = /Users/ruprechd/Programs/MPICH_GCC46/bin/mpif90
FFLAGS = -fopenmp -Wtabs -ffree-line-length-none -cpp -I/Users/ruprechd/Programs/MPICH_GCC46/include -I$(CURDIR)/../libserialio -I$(CURDIR)
LDFLAGS = -L/Users/ruprechd/Programs/MPICH_GCC46/lib -lhdf5_fortran -lhdf5 -L$(CURDIR)/../libserialio -lserialio 

#FFLAGS += -frange-check -fbacktrace -fvar-tracking -C -Wconversion -fimplicit-none -Warray-temporaries -Wcharacter-truncation -Wline-truncation -fcheck=all -g -fbounds-check

all: libspatialdiscretization.a

# Rule for library file
libspatialdiscretization.a: timediscretization_parameter.o FVMParameters.o MPIParameter.o FiniteVolumes.o DistributedIO.o mpi_space_communication.o RHSFunctions.o spatialdiscretization.o
	ar crs libspatialdiscretization.a FVMParameters.o FiniteVolumes.o RHSFunctions.o MPIParameter.o mpi_space_communication.o DistributedIO.o spatialdiscretization.o timediscretization_parameter.o

# Rules for object files

timediscretization_parameter.o: $(CURDIR)/../libserialio/libserialio.a timediscretization_parameter.f90
	$(FC) $(FFLAGS) -c timediscretization_parameter.f90

FVMParameters.o: $(CURDIR)/../libserialio/libserialio.a FiniteVolumes.o FVMParameters.f90
	$(FC) $(FFLAGS) -c FVMParameters.f90 

MPIParameter.o: FVMParameters.o FiniteVolumes.o MPIParameter.f90
	$(FC) $(FFLAGS) -c MPIParameter.f90 

#FiniteVolumes.o: FiniteVolumes.f90
#	$(FC) $(FFLAGS) -c FiniteVolumes.f90 

FiniteVolumes.o: WENO5_1D.f90
	$(FC) $(FFLAGS) -c WENO5_1D.f90 -o FiniteVolumes.o

#FiniteVolumes.o: FD_Burgers.f90
#	$(FC) $(FFLAGS) -c FD_Burgers.f90 -o FiniteVolumes.o 

#FiniteVolumes.o: WENO5_Burgers_CCOverlap.f90
#	$(FC) $(FFLAGS) -c WENO5_Burgers_CCOverlap.f90 -o FiniteVolumes.o 

#FiniteVolumes.o: WENO5_Burgers.f90
#	$(FC) $(FFLAGS) -c WENO5_Burgers.f90 -o FiniteVolumes.o

DistributedIO.o: FVMParameters.o MPIParameter.o FiniteVolumes.o DistributedIO.f90
	$(FC) $(FFLAGS) -c DistributedIO.f90

mpi_space_communication.o: FVMParameters.o MPIParameter.o FiniteVolumes.o mpi_space_communication.f90
	$(FC) $(FFLAGS) -c mpi_space_communication.f90

#RHSFunctions.o: FVMParameters.o FiniteVolumes.o MPIParameter.o DistributedIO.o mpi_space_communication.o RHSFunctions_CCoverlap.f90
#	$(FC) $(FFLAGS) -c RHSFunctions_CCoverlap.f90 -o RHSFunctions.o

RHSFunctions.o: FVMParameters.o FiniteVolumes.o MPIParameter.o DistributedIO.o mpi_space_communication.o RHSFunctions.f90
	$(FC) $(FFLAGS) -c RHSFunctions.f90 -o RHSFunctions.o

spatialdiscretization.o: FVMParameters.o MPIParameter.o FiniteVolumes.o DistributedIO.o mpi_space_communication.o RHSFunctions.o spatialdiscretization.f90
	$(FC) $(FFLAGS) -c spatialdiscretization.f90

spatialdiscretization.mod: spatialdiscretization.o

# Rule to update serial io library if need
$(CURDIR)/../libserialio/libserialio.a:
	cd $(CURDIR)/../libserialio; make

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a