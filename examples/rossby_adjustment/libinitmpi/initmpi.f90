MODULE initmpi
! This module provides access to the routines MPI_INIT_THREAD and MPI_TERMINATE. It initializes
! multithread MPI depending on the flag set when calling 'StartMPI'. The function returns the MPI_COMM_WORLD
! communicator only, everything else has to be handled in the libraries for the spatial and temporal parallelization.
!
! Daniel Ruprecht
! Institute of Computational Science, Lugano
! October 30, 2012
IMPLICIT NONE

INCLUDE 'mpif.h'

! While here the parameter type contains only one single field, to match the design of the other libs, it is also
! encapsulated into a "param" data type. This also allows to easily add more parameter later if needed.
TYPE :: MPIInitialize_Parameter
	LOGICAL :: echo_on
END TYPE

TYPE(MPIInitialize_Parameter) :: param

PRIVATE
PUBLIC :: StartMPI, TerminateMPI

CONTAINS

	SUBROUTINE StartMPI(mpi_init_thread_flag, mpi_communicator, echo_on)
	
		INTEGER, INTENT(IN)  :: mpi_init_thread_flag
		INTEGER, INTENT(OUT) :: mpi_communicator
		LOGICAL, INTENT(IN)  :: echo_on
		
		INTEGER :: mpi_thread_provided, ierr
		
		param%echo_on = echo_on
		IF (param%echo_on) WRITE(*,*) 'Starting initialization of MPI'
		! Initialize MPI
		SELECT CASE (mpi_init_thread_flag)
		
			CASE (0) ! flag==0 means no MPI at all. There is only one, unsplit computational domain and ghost-cell values are computed depending on the prescribed BC
				
				CALL MPI_INIT(ierr)
! NOTE: In the end, this mode should work completely without ever calling any MPI routine >>>>>>				
				!IF (echo_on) WRITE(*,'(A)') ' Running without MPI.'
				!mpi_communicator = -1 ! In this mode, the value should never actually be used in any MPI function call.
		
			CASE (1)! flag==1 means either pure MPI without multithreading or funneled MPI
			
				CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)
			
				IF ( mpi_thread_provided < MPI_THREAD_FUNNELED) THEN
					WRITE(*,'(A)') ' MPI multithreaded initialization failed:'
					WRITE(*,'(A, I2)') ' Value for requested MPI_THREAD_FUNNELED : ', MPI_THREAD_FUNNELED
					WRITE(*,'(A, I2)') ' Returned value from MPI_INIT_THREAD     : ', mpi_thread_provided
					WRITE(*,*) ' Now exiting.'
					STOP
				END IF

				IF (echo_on) WRITE(*,'(A)') ' MPI multithreaded initialization succesful: MPI_THREAD_FUNNELED'

			CASE (2) ! flag==2 corresponds to serialized MPI multithreading
			
				CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, mpi_thread_provided, ierr)
				
				IF ( mpi_thread_provided < MPI_THREAD_SERIALIZED) THEN
					WRITE(*,'(A)') ' MPI multithreaded initialization failed:'
					WRITE(*,'(A, I2)') ' Value for requested MPI_THREAD_SERIALIZED : ', MPI_THREAD_SERIALIZED
					WRITE(*,'(A, I2)') ' Returned value from MPI_INIT_THREAD       : ', mpi_thread_provided
					WRITE(*,*) ' Now exiting.'
					STOP					
				END IF
				
				IF (echo_on) WRITE(*,'(A)') ' MPI multithreaded initialization succesful: MPI_THREAD_SERIALIZED'
				
			
			CASE (3) ! flag==3 means full threadsafe MPI
			
				CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, mpi_thread_provided, ierr)
				
				IF ( mpi_thread_provided < MPI_THREAD_MULTIPLE) THEN
					WRITE(*,'(A)') ' MPI multithreaded initialization failed:'
					WRITE(*,'(A, I2)') ' Value for requested MPI_THREAD_MULTIPLE : ', MPI_THREAD_MULTIPLE
					WRITE(*,'(A, I2)') ' Returned value from MPI_INIT_THREAD     : ', mpi_thread_provided
					WRITE(*,*) ' Now exiting.'
					STOP					
				END IF
				IF (echo_on) WRITE(*,'(A)') ' MPI multithreaded initialization succesful: MPI_THREAD_MULTIPLE'
				
			CASE DEFAULT
					
				WRITE(*,'(A, I2)') 'Encountered unknown value for mpi_init_thread_flag, now exiting. Value found: ', mpi_init_thread_flag
				STOP
				
		END SELECT
	
		mpi_communicator = MPI_COMM_WORLD
		
	END SUBROUTINE StartMPI
	
	SUBROUTINE TerminateMPI()

		INTEGER :: ierr
	
		CALL MPI_FINALIZE(ierr)
		IF ((param%echo_on) .and. (ierr==0)) WRITE(*, '(A)') ' Successful termination of MPI'
	
	END SUBROUTINE TerminateMPI

END MODULE initmpi
