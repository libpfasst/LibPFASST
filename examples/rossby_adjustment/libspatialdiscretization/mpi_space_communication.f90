MODULE mpi_space_communication
! This module provides a collection of routines and variables
! that are required for the spatial decomposition of the solution
! vector based on decomposing the domain into subdomains and then
! assigning each subdomain to one MPI process.
!
! As the code features a hybrid parallelization, distributed memory and
! MPI in space, shared memory and OpenMP in time, each MPI process is not
! assigned to a single core but to N_p processors, where N_p is the number of
! threads used by the time parallelization.
!
! Ideally, this module should hide all MPI variables and MPI subroutine calls
! from the rest of the code, by providing updated ghost cell values in a
! "black box fashion".
!

USE FVMParameters,         only : nr_fields
USE MPIParameter,          only : cartesian_comm, &
							      ProcNeighbors, Mpi_Datatypes_AllThreads, Mpi_Datatypes_SingleThread, mpi_request, mpi_status
USE FiniteVolumes,         only : FillinGhostcells
USE omp_lib,               only : omp_get_thread_num

IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE mpicomm_parameter
	LOGICAL :: echo_on
	INTEGER :: mpi_init_thread_flag
END TYPE

TYPE(mpicomm_parameter) :: param

! Two versions are available: One takes the full 4-D buffer Q containing values for all threads and sends it, the
! other receives a 3-D buffer corresponding to the solution for a single thread and sends it.
INTERFACE PostSends
	MODULE PROCEDURE PostSends_AllThreads, PostSends_SingleThread
END INTERFACE PostSends
							
! Because in fine and coarse propagator, messages of different lengths are send, linear receive buffers are used from which the
! data is then sorted into the corresponding ghost-cell buffers		
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Qup_recv, Qdown_recv, Qleft_recv, Qright_recv, &
												 Qupleft_recv, Qupright_recv, Qdownleft_recv, Qdownright_recv
											   				   	   
CONTAINS
		
	SUBROUTINE initialize_mpi_space_communication(Nx, Ny, Nghost_max, mpi_init_thread_flag, Nthreads, echo_on)
		
		INTEGER, INTENT(IN) :: Nx, Ny, Nghost_max, mpi_init_thread_flag, Nthreads
		LOGICAL, INTENT(IN) :: echo_on
		
		INTEGER             :: thread_nr, i
		
		param%echo_on = echo_on
		param%mpi_init_thread_flag = mpi_init_thread_flag
		
! If only the master thread is handling communication, the send and receive buffers have to contain the ghost cell information
! for all threads. Accordingly, the buffers are "first touched" only by the master thread, so that memory is allocated close to
! the CPU running the master thread. This will, however, require broadcasting data through the interconnect in the socket, as the
! data from the receive buffers have to be copied into the buffers for the ghost cells, which are located close to the CPUs running
! the corresponding thread. 

	   IF ( (mpi_init_thread_flag==0) .or. (mpi_init_thread_flag==1) ) THEN
	   
			ALLOCATE(Qup_recv(   Nghost_max*Nx*nr_fields*Nthreads, 0:0))
			ALLOCATE(Qdown_recv( Nghost_max*Nx*nr_fields*Nthreads, 0:0))
			ALLOCATE(Qleft_recv( Ny*Nghost_max*nr_fields*Nthreads, 0:0))
			ALLOCATE(Qright_recv(Ny*Nghost_max*nr_fields*Nthreads, 0:0))
			
			ALLOCATE(Qupleft_recv(   nr_fields*Nthreads,0:0))
			ALLOCATE(Qupright_recv(  nr_fields*Nthreads,0:0))
			ALLOCATE(Qdownleft_recv( nr_fields*Nthreads,0:0))
			ALLOCATE(Qdownright_recv(nr_fields*Nthreads,0:0))

			Qup_recv(       :,0) = DBLE(0.0)
			Qleft_recv(     :,0) = DBLE(0.0)
			Qupleft_recv(   :,0) = DBLE(0.0)
			Qupright_recv(  :,0) = DBLE(0.0)
			Qdownleft_recv( :,0) = DBLE(0.0)
			Qdownright_recv(:,0) = DBLE(0.0)
		
		ELSE IF ( (mpi_init_thread_flag==2) .or. (mpi_init_thread_flag==3) ) THEN
		
			ALLOCATE(Qup_recv(   Nghost_max*Nx*nr_fields, 0:Nthreads-1))
			ALLOCATE(Qdown_recv( Nghost_max*Nx*nr_fields, 0:Nthreads-1))
			ALLOCATE(Qleft_recv( Ny*Nghost_max*nr_fields, 0:Nthreads-1))
			ALLOCATE(Qright_recv(Ny*Nghost_max*nr_fields, 0:Nthreads-1))
			
			ALLOCATE(Qupleft_recv(   nr_fields,0:Nthreads-1))
			ALLOCATE(Qupright_recv(  nr_fields,0:Nthreads-1))
			ALLOCATE(Qdownleft_recv( nr_fields,0:Nthreads-1))
			ALLOCATE(Qdownright_recv(nr_fields,0:Nthreads-1))

			!$OMP PARALLEL private(thread_nr)
			!$OMP DO schedule(static)
			DO i=0,Nthreads-1
				thread_nr = omp_get_thread_num()
				Qup_recv(       :,thread_nr) = DBLE(0.0)
				Qleft_recv(     :,thread_nr) = DBLE(0.0)
				Qupleft_recv(   :,thread_nr) = DBLE(0.0)
				Qupright_recv(  :,thread_nr) = DBLE(0.0)
				Qdownleft_recv( :,thread_nr) = DBLE(0.0)
				Qdownright_recv(:,thread_nr) = DBLE(0.0)
			END DO
			!$OMP END DO
			!$OMP END PARALLEL
			
		ELSE	
			WRITE(*,*) ' Encountered unknown value for flag for mpi_init_thread_flag in InitializeMPICommunication. Now exiting.'
			STOP
		END IF

	END SUBROUTINE initialize_mpi_space_communication
	
	SUBROUTINE PostSends_AllThreads(Q, Nghost)
	
		DOUBLE PRECISION, DIMENSION(:,:,:,0:), INTENT(IN)  :: Q
		INTEGER,				               INTENT(IN)  :: Nghost
			
		INTEGER :: thread_nr, ierr

		! As this subroutine is invoked only by the master thread, passing the full buffer Q from RHSFunctions,
		! the number of threads is equal to the fourth dimension of Q.
		thread_nr  = omp_get_thread_num()
	
		CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(4,:), ProcNeighbors(4), 1, cartesian_comm(thread_nr), mpi_request(1,thread_nr), ierr )
		CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(5,:), ProcNeighbors(5), 2, cartesian_comm(thread_nr), mpi_request(2,thread_nr), ierr )		
		CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(2,:), ProcNeighbors(2), 3, cartesian_comm(thread_nr), mpi_request(3,thread_nr), ierr )
		CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(7,:), ProcNeighbors(7), 4, cartesian_comm(thread_nr), mpi_request(4,thread_nr), ierr )

		IF (nr_fields>1) THEN		
			CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(1,:), ProcNeighbors(1), 5, cartesian_comm(thread_nr), mpi_request(9,thread_nr),  ierr )
			CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(3,:), ProcNeighbors(3), 6, cartesian_comm(thread_nr), mpi_request(10,thread_nr), ierr )
			CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(6,:), ProcNeighbors(6), 7, cartesian_comm(thread_nr), mpi_request(11,thread_nr), ierr )
			CALL MPI_ISEND( Q(1,1,1,0), 1, Mpi_Datatypes_AllThreads(8,:), ProcNeighbors(8), 8, cartesian_comm(thread_nr), mpi_request(12,thread_nr), ierr )
		END IF			
		
		IF (param%echo_on) WRITE(*,'(A)') ' AllThreads sends prepared, posted MPI_ISENDs... '
																						
	END SUBROUTINE PostSends_AllThreads
		
	SUBROUTINE PostSends_SingleThread(Q, Nghost)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		INTEGER,				            INTENT(IN)  :: Nghost
			
		INTEGER :: thread_nr, ierr

		thread_nr  = omp_get_thread_num()

		! As this subroutine is invokded only by the master thread, passing the full buffer Q from RHSFunctions,
		! the number of threads is equal to the fourth dimension of Q.	
		CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(4,:), ProcNeighbors(4), 1, cartesian_comm(thread_nr), mpi_request(1,thread_nr), ierr )
		CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(5,:), ProcNeighbors(5), 2, cartesian_comm(thread_nr), mpi_request(2,thread_nr), ierr )		
		CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(2,:), ProcNeighbors(2), 3, cartesian_comm(thread_nr), mpi_request(3,thread_nr), ierr )
		CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(7,:), ProcNeighbors(7), 4, cartesian_comm(thread_nr), mpi_request(4,thread_nr), ierr )

		IF (nr_fields>1) THEN		
			CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(1,:), ProcNeighbors(1), 5, cartesian_comm(thread_nr), mpi_request(9,thread_nr),  ierr )
			CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(3,:), ProcNeighbors(3), 6, cartesian_comm(thread_nr), mpi_request(10,thread_nr), ierr )
			CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(6,:), ProcNeighbors(6), 7, cartesian_comm(thread_nr), mpi_request(11,thread_nr), ierr )
			CALL MPI_ISEND( Q(1,1,1), 1, Mpi_Datatypes_SingleThread(8,:), ProcNeighbors(8), 8, cartesian_comm(thread_nr), mpi_request(12,thread_nr), ierr )
		END IF			
		
		IF (param%echo_on) WRITE(*,'(A)') ' SingleThread sends prepared, posted MPI_ISENDs...  '
																						
	END SUBROUTINE PostSends_SingleThread		

	SUBROUTINE PostReceives(total_length_x, total_length_y, nr_threads)
		
		INTEGER, INTENT(IN) :: total_length_x, total_length_y, nr_threads
	
		INTEGER :: ierr, thread_nr
		
		ierr=-1
		
		thread_nr      = omp_get_thread_num()
		
		! The "send to right" message corresponds to the "receive from left" message, hence the tags are exchanged. Same for the "send upwards", "receive from downwards"						
		CALL MPI_IRECV( Qleft_recv( 1:total_length_x, thread_nr), total_length_x, MPI_DOUBLE_PRECISION, ProcNeighbors(4), 2, cartesian_comm(thread_nr), mpi_request(5,thread_nr), ierr)
		CALL MPI_IRECV( Qright_recv(1:total_length_x, thread_nr), total_length_x, MPI_DOUBLE_PRECISION, ProcNeighbors(5), 1, cartesian_comm(thread_nr), mpi_request(6,thread_nr), ierr)		
		CALL MPI_IRECV( Qup_recv(   1:total_length_y, thread_nr), total_length_y, MPI_DOUBLE_PRECISION, ProcNeighbors(2), 4, cartesian_comm(thread_nr), mpi_request(7,thread_nr), ierr)
		CALL MPI_IRECV( Qdown_recv( 1:total_length_y, thread_nr), total_length_y, MPI_DOUBLE_PRECISION, ProcNeighbors(7), 3, cartesian_comm(thread_nr), mpi_request(8,thread_nr), ierr)
				
		IF (nr_fields > 1) THEN
			CALL MPI_IRECV( Qupleft_recv(   1:nr_fields*nr_threads,0), nr_fields*nr_threads, MPI_DOUBLE_PRECISION, ProcNeighbors(1), 8, cartesian_comm(thread_nr), mpi_request(13,thread_nr), ierr)
			CALL MPI_IRECV( Qupright_recv(  1:nr_fields*nr_threads,0), nr_fields*nr_threads, MPI_DOUBLE_PRECISION, ProcNeighbors(3), 7, cartesian_comm(thread_nr), mpi_request(14,thread_nr), ierr)
			CALL MPI_IRECV( Qdownleft_recv( 1:nr_fields*nr_threads,0), nr_fields*nr_threads, MPI_DOUBLE_PRECISION, ProcNeighbors(6), 6, cartesian_comm(thread_nr), mpi_request(15,thread_nr), ierr)
			CALL MPI_IRECV( Qdownright_recv(1:nr_fields*nr_threads,0), nr_fields*nr_threads, MPI_DOUBLE_PRECISION, ProcNeighbors(8), 5, cartesian_comm(thread_nr), mpi_request(16,thread_nr), ierr)
		END IF

		IF (param%echo_on) WRITE(*,'(A)') 'Placed non-blocking receives for boundary values...'

	END SUBROUTINE PostReceives
	
	!
	! To guarantee that non-blocking communication has finished, call the corresponding MPI_WAIT routines
	!
	SUBROUTINE WaitForCommunication(Nx, Ny, Nghost, nr_threads)
	
		INTEGER, INTENT(IN) :: Nx, Ny, Nghost, nr_threads

		INTEGER :: ierr, i, j, k, counter, thread_nr
		
		INTEGER :: status(MPI_STATUS_SIZE,8)
		ierr=-1
		
		thread_nr = omp_get_thread_num()

		! Now post a wait for the non-blocking sends posted in the send routin
		IF (nr_fields > 1) THEN
			CALL MPI_WAITALL(16, mpi_request(:,thread_nr), mpi_status(:,:,thread_nr), ierr)
		ELSE
!			CALL MPI_WAITALL(8, mpi_request(:,thread_nr), mpi_status(:,:,thread_nr), ierr)		
						CALL MPI_WAITALL(8, mpi_request(:,thread_nr), mpi_status(:,:,thread_nr), ierr)		

		END IF

    	IF (param%echo_on) WRITE(*,'(A, I3)') ' Error flag of MPI_WAITALL : ', ierr

		IF (param%mpi_init_thread_flag == 1) THEN
		
				CALL FillinGhostcells(Qleft_recv(:,0), Qright_recv(:,0), Qup_recv(:,0), Qdown_recv(:,0), Qupleft_recv(:,0), Qupright_recv(:,0), Qdownleft_recv(:,0), Qdownright_recv(:,0), Nx, Ny, Nghost)
		
		ELSE IF ((param%mpi_init_thread_flag == 2) .or. (param%mpi_init_thread_flag==3)) THEN
		
			CALL FillinGhostcells(Qleft_recv(:,thread_nr), Qright_recv(:,thread_nr), Qup_recv(:,thread_nr), Qdown_recv(:,thread_nr), Qupleft_recv(:,thread_nr), &
						Qupright_recv(:,thread_nr), Qdownleft_recv(:,thread_nr), Qdownright_recv(:,thread_nr), Nx, Ny, Nghost)
		END IF
				
	END SUBROUTINE
	
	
END MODULE mpi_space_communication