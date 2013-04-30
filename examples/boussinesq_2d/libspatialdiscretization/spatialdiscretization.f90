MODULE spatialdiscretization
! This is the encompassing module for the spatial parallelization, including modules for MPI communication between subdomains
! and distributed IO.
!
! The driver function should only use this module, through which all required functions and variable are then made available. Whatever is or should not be 
! needed is declated PRIVATE below.
!
! Daniel Ruprecht
! Institute of Computational Science, Lugano
! November 2, 2012
!
USE FVMParameters,           only : Nx, Ny, dim, ReadFVMParameter, c_s, stabFreq, grav
USE timediscretization_parameter, only : ReadIntegrationParameter, global_tend, maxit, Nsteps_fine_total
USE FiniteVolumes,           only : InitializeFiniteVolumes, nr_fields
USE RHSFunctions,            only : InitializeRHSFunctions, RHS, RHS_coarse, InitializeLinearAdvectionVelocity, Timer_BOne, Timer_BTwo, NrCalls_BOne, NrCalls_BTwo
USE MPIParameter,            only : InitializeMPI, FinalizeSpatialMPI, MPI_WTIME, myrank, nprocs, MPI_DOUBLE_PRECISION, cartesian_comm
USE mpi_space_communication, only : initialize_mpi_space_communication
USE DistributedIO,           only : InitializeDistributedIO, CloseDistributedIO, WriteData, WriteFinalData, WriteTimerData, ReadInitialValue, ReadLinearAdvectionVelocity

IMPLICIT NONE

PRIVATE
PUBLIC :: dim, Nx, Ny, nr_fields, ReadFVMParameter  ! from FVMParameters and FiniteVolumes
!PUBLIC :: FinalizeSpatialMPI, myrank, nprocs, MPI_DOUBLE_PRECISION, MPI_WTIME, cartesian_comm ! from MPIParameter
PUBLIC :: RHS, RHS_coarse !, InitializeLinearAdvectionVelocity, Timer_BOne, Timer_BTwo, NrCalls_BOne, NrCalls_BTwo ! from RHSFunctipns
PUBLIC :: ReadInitialValue, WriteData !,CloseDistributedIO, WriteFinalData, WriteTimerData, ReadLinearAdvectionVelocity ! from DistributedIO
PUBLIC :: InitializeSpatialDiscretization, CloseSpatialDiscretization ! From this module
PUBLIC :: ReadIntegrationParameter, global_tend, maxit, Nsteps_fine_total

CONTAINS

	SUBROUTINE InitializeSpatialDiscretization(maxit, Nparareal_restarts, mpi_init_thread_flag, mpi_communicator, Nthreads, echo_on, dim)
		LOGICAL, INTENT(IN) :: echo_on
		INTEGER, INTENT(IN) :: Nthreads, maxit, Nparareal_restarts, mpi_init_thread_flag, mpi_communicator
		INTEGER, INTENT(OUT) :: dim
			
		CALL InitializeFiniteVolumes(Nx, Ny,                                                                         Nthreads = Nthreads, mpi_init_thread = mpi_init_thread_flag, echo_on = echo_on, c_s = c_s, stabFreq = stabFreq, grav = grav)
		CALL InitializeMPI(mpi_communicator,                                                                         Nthreads = Nthreads, mpi_init_thread = mpi_init_thread_flag, echo_on = echo_on)
		CALL InitializeRHSFunctions(mpi_init_thread_flag,                                                            Nthreads = Nthreads, echo_on = echo_on)
		CALL initialize_mpi_space_communication(Nx, Ny, Nghost_max = 3, mpi_init_thread_flag = mpi_init_thread_flag, Nthreads = Nthreads, echo_on = echo_on)
		CALL InitializeDistributedIO(maxit, Nparareal_restarts,                                                      Nthreads = Nthreads, echo_on = echo_on)
		CALL InitializeLinearAdvectionVelocity()

	END SUBROUTINE InitializeSpatialDiscretization
	 
        SUBROUTINE CloseSpatialDiscretization()
          CALL CloseDistributedIO()
          CALL FinalizeSpatialMPI()
        END SUBROUTINE CloseSpatialDiscretization

END MODULE spatialdiscretization
