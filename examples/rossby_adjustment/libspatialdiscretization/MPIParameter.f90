MODULE MPIParameter
! This module holds and provides all variables related to MPI communication. It also
! provides subroutines for initializing and terminating MPI.
!
! The module does not provide any communication routines, these are all included in the
! "MPICommunication" module.
!
! The reason for separating the variables and the communication routines is that the communication 
! subroutines require information about the
! size of the broadcasted buffer. This information is, however, read from a HDF5 data file
! and in order to properly set up HDF5 for distributed reading, a MPI communicator is needed.
 
USE FVMParameters,        only : Nx, Ny, nprocs_x, nprocs_y, BC, order_fine_advection, order_coarse_advection, order_fine_sound, order_coarse_sound
								  
USE FiniteVolumes,        only : GetMpiDatatypePar, nr_fields
								  
IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE mpiparam_parameter
	INTEGER :: Nthreads
	LOGICAL :: echo_on
END TYPE

TYPE(mpiparam_parameter) :: param

INTEGER :: ierr, nprocs, myrank, cart_coords(2), Nghost, Nghost_coarse

INTEGER, ALLOCATABLE, DIMENSION(:)     :: cartesian_comm
INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: mpi_request
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: mpi_status

! Ranks of the eight processors handling the eight subdomains around the subdomain
! belonging to the current process. Values are initialized with -1, so if later the
! values are not properly set, MPI will throw an error as no rank -1 should exist.
!
! Indices correspond to location in following way:
!
! 1   2   3
! 4   P   5
! 6   7   8
!
INTEGER, DIMENSION(8) :: ProcNeighbors = -1

! Define four MPI datatypes that correspond the the halo-cells of the neighboring processes.
! "UpperRowsType" are the Nghost upper rows which are send to the process above the current one
! etc.
INTEGER, DIMENSION(8,MPI_STATUS_SIZE) :: Mpi_Datatypes_AllThreads, Mpi_Datatypes_SingleThread
		   
! For the cartesian communicator: Periodic in both directions and no reordering		   
INTEGER :: periodic(2) = (/ 1 , 1 /) , reorder=0 		   
		   	   		   
CONTAINS

	SUBROUTINE InitializeMPI(mpi_communicator, Nthreads, mpi_init_thread, echo_on)
		INTEGER, INTENT(IN) :: mpi_communicator, Nthreads, mpi_init_thread
		LOGICAL, INTENT(IN) :: echo_on
		
		INTEGER               :: i
					
		param%echo_on = echo_on	
		param%Nthreads = Nthreads
		
!		IF (nr_fields > 1) THEN
!			IF (mpi_init_thread==0) THEN
!				ALLOCATE(cartesian_comm(0:0))
!			! If corner ghost-cells are communicated, there are 2x8 = 16 sends and receives
!			ELSE IF (mpi_init_thread==1) THEN
!				ALLOCATE(mpi_request(16,0:0))
!				ALLOCATE(mpi_status(MPI_STATUS_SIZE, 16, 0:0))
!				ALLOCATE(cartesian_comm(0:0))
!			ELSE IF ((mpi_init_thread==2) .or. (mpi_init_thread==3)) THEN
!				ALLOCATE(mpi_request(16,0:Nthreads-1))
!				ALLOCATE(mpi_status(MPI_STATUS_SIZE, 16, 0:Nthreads-1))	
!				! to make sure that each thread communicates with its corresponding counterpart in other MPI processes, for the THREAD_SERIALIZED
!				! paradigm there are as much communicators as there are threads, each communicator providing the context for threads of the same thread number
!				ALLOCATE(cartesian_comm(0:Nthreads-1))		
!			END IF
!		ELSE
			! Without corners, there are 2x4 = 8 sends and receives
			IF (mpi_init_thread==0) THEN
				ALLOCATE(cartesian_comm(0:0))
		
			ELSE IF (mpi_init_thread==1) THEN
			
				ALLOCATE(mpi_request(8,0:0))
				ALLOCATE(mpi_status(MPI_STATUS_SIZE, 8, 0:0))
				ALLOCATE(cartesian_comm(0:0)) 
			
			ELSE IF ((mpi_init_thread==2) .or. (mpi_init_thread==3)) THEN
		
				ALLOCATE(mpi_request(8,0:Nthreads-1))
				ALLOCATE(mpi_status(MPI_STATUS_SIZE, 8, 0:Nthreads-1))
				ALLOCATE(cartesian_comm(0:Nthreads-1))
			
			END IF
		
!		END IF
				
		! Halo size depends on order and thus support of stencil		
		SELECT CASE (nr_fields)
		
			CASE (1)
				
				!Nghost_coarse = CEILING( DBLE(order_coarse_advection)/2.0 )
				!Nghost        = CEILING( DBLE(order_fine_advection)/2.0 )

				! WENO-5 needs 3 cell-wide HALO
				Nghost        = 3
				Nghost_coarse = 3
! >>> NOTE: These values should, as nr_fields, be provided by the FiniteVolume module. 	(change later)			
								
			CASE (2)
			
				Nghost        = 3
				Nghost_coarse = 3
				
				
			CASE (3)
				
				! For acoustic-advection and Boussinesq, at least a two cell halo has to be communicated in order for the 
				! divergence damping to be applicable.
				Nghost_coarse = MAX( CEILING(DBLE( MAX(order_coarse_advection, order_coarse_sound) )/2.0) , 2 )
				Nghost        = MAX( CEILING(DBLE( MAX(order_fine_advection, order_fine_sound) )/2.0) , 2 )

			CASE (4)
		
				Nghost_coarse = MAX( CEILING(DBLE( MAX(order_coarse_advection, order_coarse_sound) )/2.0) , 2 )
				Nghost        = MAX( CEILING(DBLE( MAX(order_fine_advection, order_fine_sound) )/2.0) , 2 )		
				
			CASE DEFAULT	
		
				WRITE(*,*) 'No values for Nghost set for chosen number of nr_fieds. Now exiting!'
				STOP
				
		END SELECT
		
		! First, set topology of MPI according to BC (cased 2 and 3 not yet implemented)
		SELECT CASE (BC)
		
			CASE (1)
				! All periodic boundary condition
				periodic = (/1,1/)
			CASE (2)
				! Outflow in all directions, hence no periodicity
				periodic = (/0,0/)
			CASE (3)
				! Solid walls in vertical, outflow in horizontal
				periodic = (/0,1/)

		END SELECT
		
		! Retrieve number of processes
		CALL MPI_COMM_SIZE(mpi_communicator, nprocs, ierr)
		
		IF ((mpi_init_thread == 0) .and. (nprocs > 1)) THEN
			WRITE(*,'(A)') 'Compiled with mpi_mode=0 but found number of processes > 1. Now exiting.' 
			STOP
		END IF

		! Verify that the number of processors available matches the number prescribed for the topology
		IF (nprocs .ne. nprocs_x*nprocs_y) THEN
			WRITE(*,'(A)') ' Number of processes does not match number of processes prescribed for MPI topology! Now exiting.'
			CALL MPI_ABORT(mpi_communicator, ierr)
			STOP
		END IF
			
		! Retrieve rank of current process
		CALL MPI_COMM_RANK(mpi_communicator, myrank, ierr)
			
		! Create cartesian topology
		CALL MPI_CART_CREATE(mpi_communicator, 2, (/ nprocs_y, nprocs_x /), periodic, reorder, cartesian_comm(0), ierr)
			
		! Determine ranks of neighboring processes.
		CALL MPI_CART_SHIFT(cartesian_comm(0), 0, 1, ProcNeighbors(2), ProcNeighbors(7),  ierr)
		CALL MPI_CART_SHIFT(cartesian_comm(0), 1, 1, ProcNeighbors(4), ProcNeighbors(5), ierr)
			
!		IF (nr_fields > 1) THEN
!			CALL DetermineCornerProcessRanks()
!		END IF
		
		! Determine coordinates of current process in cartesian topology
		CALL MPI_CART_COORDS(cartesian_comm(0), myrank, 2, cart_coords, ierr)
		
		! If the SERIALIZED paradigm is used, create now a number of copies of the cartesian communicator, equal to the number of threads spawned by each MPI process
		IF ((mpi_init_thread==2) .or. (mpi_init_thread==3)) THEN
			DO i=1,Nthreads-1
				CALL MPI_COMM_DUP(cartesian_comm(0), cartesian_comm(i), ierr)
			END DO
		END IF
		
		! Generate MPI datatypes with values provided by the FiniteVolume module to be able to select the exact values out of a solution buffer
		! that represent ghost-cell values for a neighboring sub-domain
		! CALL GenerateMpiDatatypes()
		
		IF (param%echo_on) THEN
			CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
			CALL SLEEP(myrank) ! wait so that processes report with increasing rank
			WRITE(*,*)
			WRITE(*, '(A, I2)') ' Initialization of MPI successful. Used processes : ', nprocs
			WRITE(*,'(A)') ' Cartesian topology for MPI created: '
			WRITE(*,'(A,I2)') ' Processors in x direction : ', nprocs_x
			WRITE(*,'(A,I2)') ' Processors in y direction : ', nprocs_y
			WRITE(*,*)
			WRITE(*,'(A,I2,A)') ' Rank ', myrank, ' has neighbors : '
			WRITE(*,'(I2)',advance="no") ProcNeighbors(1)
			WRITE(*,'(I2)',advance="no") ProcNeighbors(2)
			WRITE(*,'(I2)') ProcNeighbors(3)
			WRITE(*,'(I2)',advance="no") ProcNeighbors(4)
			WRITE(*,'(I2)',advance="no") myrank
			WRITE(*,'(I2)') ProcNeighbors(5)
			WRITE(*,'(I2)',advance="no") ProcNeighbors(6)
			WRITE(*,'(I2)',advance="no") ProcNeighbors(7)
			WRITE(*,'(I2)') ProcNeighbors(8)
			CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
		END IF				

		CONTAINS 
		
			SUBROUTINE DetermineCornerProcessRanks()
			
				INTEGER, DIMENSION(2) :: send, recv

				! Ranks for processes at corners are not readily available from MPI_CART_SHIFT, so every process broadcasts its upper and lower
				! process rank to its left and right neighbor, so that every process knows the ranks of the processes at the corners.		
				send = (/ ProcNeighbors(2) , ProcNeighbors(7) /)
				CALL MPI_ISEND(send, 2, MPI_INTEGER, ProcNeighbors(4), 1, cartesian_comm, mpi_request(1,0),   ierr)
				CALL MPI_RECV( recv, 2, MPI_INTEGER, ProcNeighbors(5), 1, cartesian_comm, mpi_status(:,1,0),  ierr)
				CALL MPI_WAIT(mpi_request(1,0), mpi_status(:,2,0), ierr)
				
				ProcNeighbors(3) = recv(1)
				ProcNeighbors(8) = recv(2)
				
				CALL MPI_ISEND(send, 2, MPI_INTEGER, ProcNeighbors(5), 2, cartesian_comm, mpi_request(1,0),   ierr)
				CALL MPI_RECV( recv, 2, MPI_INTEGER, ProcNeighbors(4), 2, cartesian_comm, mpi_status(:,1,0),  ierr)
				CALL MPI_WAIT(mpi_request(1,0), mpi_status(:,2,0), ierr)
				
				ProcNeighbors(1) = recv(1)
				ProcNeighbors(6) = recv(2)
						
			END SUBROUTINE DetermineCornerProcessRanks
		
			SUBROUTINE GenerateMpiDatatypes()
			
				INTEGER, ALLOCATABLE, DIMENSION(:,:) :: blocklengths, indices
				INTEGER, DIMENSION(8)                :: length, length_singleThread
				
				! Generate and submit required MPI datatypes

				ALLOCATE(blocklengths(Nx*nr_fields*param%Nthreads,8))
				ALLOCATE(indices(     Nx*nr_fields*param%Nthreads,8))
				
				! Indices are according to following layout:
				!   1  2  3
				!   4  #  5 
				!   6  7  8
								
				CALL GetMpiDatatypePar(Nx, Ny, Nghost, blocklengths, indices, length, length_singleThread)
																																																																																																				
				! Now create and commit all dataypes		
				DO i=1,8
					CALL CreateCommitDatatype(blocklengths(1:length(i),i),              indices(1:length(i),i),              i, 0)
					CALL CreateCommitDatatype(blocklengths(1:length_singleThread(i),i), indices(1:length_singleThread(i),i), i, 1)				
				END DO				
																					
				DEALLOCATE(blocklengths, indices)
				
			END SUBROUTINE GenerateMpiDatatypes
			
			SUBROUTINE CreateCommitDatatype(blocklengths, indices, datatype_index, all_or_one_thread)
					INTEGER, DIMENSION(:), INTENT(IN) :: blocklengths, indices
					INTEGER,               INTENT(IN) :: datatype_index
					INTEGER,               INTENT(IN) :: all_or_one_thread
					
					INTEGER :: length, ierr
					
					length = SIZE(blocklengths)
					
					IF (all_or_one_thread==0) THEN
						CALL MPI_TYPE_INDEXED(length, blocklengths, indices, MPI_DOUBLE_PRECISION, Mpi_Datatypes_AllThreads(datatype_index,:), ierr)
						CALL MPI_TYPE_COMMIT(Mpi_Datatypes_AllThreads(datatype_index,:), ierr)					
					ELSE IF (all_or_one_thread==1) THEN
						CALL MPI_TYPE_INDEXED(length, blocklengths, indices, MPI_DOUBLE_PRECISION, Mpi_Datatypes_SingleThread(datatype_index,:), ierr)
						CALL MPI_TYPE_COMMIT(Mpi_Datatypes_SingleThread(datatype_index,:), ierr)					
					ELSE
						WRITE(*,*) 'Variable all_or_one_thread has to be either 0 or 1. Encountered other value, now exiting'
						STOP
					END IF
					
			END SUBROUTINE CreateCommitDatatype
	
	END SUBROUTINE InitializeMPI
	
	SUBROUTINE FinalizeSpatialMPI
		INTEGER :: ierr
		
		CALL MPI_COMM_FREE(cartesian_comm, ierr)

		IF ((param%echo_on) .and. (ierr==0)) THEN
			WRITE(*, '(A)') ' Successful termination of spatial communicator'
		END IF

	END SUBROUTINE
	
END MODULE MPIParameter
