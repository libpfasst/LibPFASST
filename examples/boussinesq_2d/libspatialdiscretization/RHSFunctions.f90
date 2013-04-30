MODULE RHSFunctions

#include <preprocessor.f90>

!
! This module provides the routines 'RHS' and 'RHS_coarse' which implement the interface 'RHS1'
! for the right hand side functions defined in the module 'TimeIntegrators'.
! 
! It completely hides the details of the underlying spatial discretization, provided by the module
! 'FiniteVolumes', from the employed integration scheme.
!
! In a hybrid parallelization, the module also performs all necessary communication between different processes:
! Before each RHS evaluation, a number of RECEIVES is posted to obtain the up-to-date ghost cell values from
! the neighbouring subdomains. After each RHS evaluation, the corresponding non-blocking SENDs are issued.
!
!

	USE FVMParameters,           only : dim, dim_coarse, Nx, Ny, Nx_coarse, Ny_coarse, order_coarse_advection, order_coarse_sound, &
							            order_fine_advection, order_fine_sound, dx, dy, dx_coarse, dy_coarse, BC, nu, nu_coarse
	USE FiniteVolumes,           only : GetRHS, PackSolution, UnpackSolution, buffer_layout, nr_fields, FillGhostCells
	USE omp_lib,                 only : omp_get_thread_num
	USE MPIParameter,            only : Nghost, Nghost_coarse, MPI_WTIME
	USE mpi_space_communication, only : PostReceives, WaitForCommunication, PostSends
	USE DistributedIO,           only : ReadLinearAdvectionVelocity
	
	IMPLICIT NONE

	TYPE rhs_parameter
		LOGICAL :: echo_on
		INTEGER :: mpi_init_thread_flag, Nthreads
	END TYPE
	
	TYPE(rhs_parameter) :: param
	
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: Uadv, Vadv, Uadv_coarse, Vadv_coarse
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ
	
	DOUBLE PRECISION, ALLOCATABLE, SAVE, DIMENSION(:) :: Timer_BOne, Timer_BTwo
	DOUBLE PRECISION, ALLOCATABLE,       DIMENSION(:) :: Timer_Temp
	INTEGER,          ALLOCATABLE, SAVE, DIMENSION(:) :: NrCalls_BOne, NrCalls_BTwo

	PRIVATE
	PUBLIC :: InitializeRHSFunctions, InitializeLinearAdvectionVelocity, RHS, RHS_coarse, Timer_BOne, Timer_BTwo, NrCalls_BOne, NrCalls_BTwo

	CONTAINS
		
		SUBROUTINE InitializeRHSFunctions(mpi_init_thread_flag, Nthreads, echo_on)
		
			INTEGER, INTENT(IN) :: mpi_init_thread_flag, Nthreads
			LOGICAL, INTENT(IN) :: echo_on
			
			INTEGER :: thread_nr, i
			
			param%echo_on = echo_on
			param%Nthreads = Nthreads
			param%mpi_init_thread_flag = mpi_init_thread_flag
			
			! buffer_layout determines which dimension counts through the different solution components
												
			IF (buffer_layout==0) THEN									
				ALLOCATE(Q( Ny, Nx, nr_fields, 0:Nthreads-1))
				ALLOCATE(RQ(Ny, Nx, nr_fields, 0:Nthreads-1))
			ELSE IF (buffer_layout==1) THEN
				ALLOCATE(Q(nr_fields, Ny, Nx, 0:Nthreads-1))
				ALLOCATE(RQ(nr_fields, Ny, Nx, 0:Nthreads-1))
			ELSE
				WRITE(*,*) 'Unknown flag for buffer_layout. Now exiting.'
				STOP
			END IF
			
			ALLOCATE(Uadv(Ny, Nx+1, 0:Nthreads-1))
			ALLOCATE(Vadv(Ny+1, Nx, 0:Nthreads-1))
			ALLOCATE(Uadv_coarse(Ny_coarse, Nx_coarse+1, 0:Nthreads-1))
			ALLOCATE(Vadv_coarse(Ny_coarse+1, Nx_coarse, 0:Nthreads-1))
			
			ALLOCATE(Timer_BOne(0:Nthreads-1))
			ALLOCATE(Timer_BTwo(0:Nthreads-1))
			ALLOCATE(Timer_Temp(0:Nthreads-1))
			ALLOCATE(NrCalls_BOne(0:Nthreads-1))
			ALLOCATE(NrCalls_BTwo(0:Nthreads-1))
			
			!$OMP PARALLEL private(thread_nr)
			!$OMP DO
			DO i=0,Nthreads-1
				thread_nr = omp_get_thread_num()
				Q(:,:,:,thread_nr)     = 0.0
				RQ(:,:,:,thread_nr)    = 0.0
				Uadv(:,:,thread_nr)    = 0.0
				Vadv(:,:,thread_nr)    = 0.0
				Uadv_coarse(:,:,thread_nr) = 0.0
				Vadv_coarse(:,:,thread_nr) = 0.0
				Timer_BOne(thread_nr)  = 0.0
				Timer_BTwo(thread_nr)  = 0.0
				Timer_Temp(thread_nr)  = 0.0
				NrCalls_BOne(thread_nr)= 0
				NrCalls_BTwo(thread_nr)= 0
			END DO
			!$OMP END DO
			!$OMP END PARALLEL
			
			IF (param%echo_on) WRITE(*, '(A,I2)') ' Module RHSFunctions successfully initialized for Nthreads =  ', Nthreads
			
		END SUBROUTINE InitializeRHSFunctions
		
		SUBROUTINE InitializeLinearAdvectionVelocity()
		
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Uadv_in, Vadv_in
			
			INTEGER :: thread_nr, i
			
			! to avoid confusion of the TAU instrumentor, we need one dummy line before the first OpenMP directive
			thread_nr = -1
			
			ALLOCATE(Uadv_in(Ny, Nx+1))
			ALLOCATE(Vadv_in(Ny+1, Nx))
			CALL ReadLinearAdvectionVelocity(Uadv_in, Vadv_in)
			
			!$OMP PARALLEL private(thread_nr) 
			!$OMP DO
			DO i=0,param%Nthreads-1
				thread_nr = omp_get_thread_num()
				Uadv(:,:,thread_nr)        = Uadv_in
				Vadv(:,:,thread_nr)        = Vadv_in
				Uadv_coarse(:,:,thread_nr) = Uadv_in
				Vadv_coarse(:,:,thread_nr) = Vadv_in
			END DO
			!$OMP END DO
			!$OMP END PARALLEL
			
			DEALLOCATE(Uadv_in, Vadv_in)
			
		END SUBROUTINE InitializeLinearAdvectionVelocity
	!	
	! NOTE: The in and out buffers Y and YDOT have to be of size dim = Ny*Nx*nr_fields, corresponding to
	! Ny x Nx cell values for every solution component. THERE IS NO CONSISTENCY CHECK HERE! If the sizes of the buffers
	! do not match Nx, Ny, nr_fields, the result will be either wrong or a segmentation fault.
	!	
		
		SUBROUTINE RHS(Y, YDOT, dt, par1, par2)
		
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN)  :: Y
			DOUBLE PRECISION,               INTENT(IN)  :: dt
			CHARACTER(len=1),               INTENT(IN)  :: par1, par2
			DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: YDOT
			
			INTEGER :: thread_nr, total_length_x, total_length_y
			
			! dummy line to avoid problems with TAU
			thread_nr = -1
			
			IF (param%echo_on) WRITE(*,'(A)') ' Starting RHS evaluation ... '
			
			thread_nr  = omp_get_thread_num()

			CALL UnpackSolution( Q(:,:,:,thread_nr), Y, nr_fields, Ny, Nx)
#if(mpi_mode==0)
			! No MPI, simply fill ghost cell values from own solution buffer according to BC
			CALL FillGhostCells(Q(:,:,:,thread_nr), Nghost, BC, thread_nr)
				
#elif(mpi_mode==1)	
			! All MPI communication is funneled through the master thread
			CALL BARRIER_ONE_FINE_PROFILING()
	
			!$OMP MASTER
				! Compute length of send messages
				total_length_x = Ny*Nghost*nr_fields*param%Nthreads
				total_length_y = Nghost*Nx*nr_fields*param%Nthreads
				CALL PostReceives(total_length_x,  total_length_y, param%Nthreads)
				CALL PostSends(Q, Nghost)
				CALL WaitForCommunication(Nx, Ny, Nghost, param%Nthreads)
			!$OMP END MASTER

			CALL BARRIER_TWO_FINE_PROFILING()
				
#elif(mpi_mode==2)
			
			total_length_x = Ny*Nghost*nr_fields
			total_length_y = Nghost*Nx*nr_fields
				
			!$OMP CRITICAL
			CALL PostSends(Q(:,:,:,thread_nr), Nghost)
			CALL PostReceives(total_length_x,  total_length_y, 1)								
			CALL WaitForCommunication(Nx, Ny, Nghost, 1)
			!$OMP END CRITICAL
							
#elif(mpi_mode==3)			
			total_length_x = Ny*Nghost*nr_fields
			total_length_y = Nghost*Nx*nr_fields
			
			CALL PostSends(Q(:,:,:,thread_nr), Nghost)
			CALL PostReceives(total_length_x,  total_length_y, 1)								
			CALL WaitForCommunication(Nx, Ny, Nghost, 1)							
#endif


			! Now that the ghost cell value are up-to-date, evaluate RHS function			
			CALL GetRHS(Q(:,:,:,thread_nr), order_fine_advection, order_fine_sound, RQ(:,:,:,thread_nr), dx, dy, dt, nu)

			CALL PackSolution( YDOT, RQ(:,:,:,thread_nr), nr_fields, Ny, Nx)		

			IF (param%echo_on) WRITE(*,*) ' Finished RHS evaluation. '
			
			CONTAINS								
		
				SUBROUTINE BARRIER_ONE_FINE_PROFILING()
				
					INTEGER :: thread_nr 

					thread_nr = omp_get_thread_num()
					Timer_Temp(thread_nr) = MPI_WTIME()
					
					!$OMP BARRIER

					Timer_BOne(thread_nr)   = Timer_BOne(thread_nr) + MPI_WTIME() - Timer_Temp(thread_nr)
					NrCalls_BOne(thread_nr) = NrCalls_BOne(thread_nr) + 1

					
				END SUBROUTINE BARRIER_ONE_FINE_PROFILING
				
				SUBROUTINE BARRIER_TWO_FINE_PROFILING()
				
					INTEGER :: thread_nr

					thread_nr = omp_get_thread_num()
					Timer_Temp(thread_nr) = MPI_WTIME()
					
					!$OMP BARRIER
					
					Timer_BTwo(thread_nr)   = Timer_BTwo(thread_nr) + MPI_WTIME() - Timer_Temp(thread_nr)
					NrCalls_BTwo(thread_nr) = NrCalls_BTwo(thread_nr) + 1
					
				END SUBROUTINE BARRIER_TWO_FINE_PROFILING
																							
		END SUBROUTINE RHS
		
		!
		! Provide an implementation of the interface RHS defined in the time integrator methods library for the coarse propagator.
		!					
		SUBROUTINE RHS_coarse(Y, YDOT, dt, par1, par2)
		
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN)  :: Y
			DOUBLE PRECISION,               INTENT(IN)  :: dt
			CHARACTER(len=1),               INTENT(IN)  :: par1, par2
			DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: YDOT
			
			INTEGER :: thread_nr, total_length_x, total_length_y
			
			thread_nr = -1 ! Dummy line to avoid problems with TAU profiler
			
			IF (param%echo_on) WRITE(*,'(A)') ' Starting RHS_coarse evaluation ... '
			
			thread_nr = omp_get_thread_num()
			
			CALL UnpackSolution(Q(:,:,:,thread_nr), Y, nr_fields, Ny_coarse, Nx_coarse)
				
#if(mpi_mode==0)
			! No MPI, simply fill ghost cell values from own solution buffer according to BC
			CALL FillGhostCells(Q(:,:,:,thread_nr), Nghost, BC, thread_nr)
								
#elif(mpi_mode==1)

			! In non-pipelined Parareal, the coarse propagator is run only
			! outside of the multi-threaded region, so no OMP directives are required here.
			total_length_x = Ny*Nghost*nr_fields*1
			total_length_y = Nx*Nghost*nr_fields*1					
			CALL PostReceives(total_length_x,  total_length_y, 1)
			CALL PostSends(Q(:,:,:,0), Nghost)
			CALL WaitForCommunication(Nx_coarse, Ny_coarse, Nghost, 1)
					
					
#elif (mpi_mode==2)	

			! Compute length of send messages: SIZE(Q,1)*Nghost*SIZE(Q,3)*SIZE(Q,4) and
			! Nghost*SIZE(Q,2)*SIZE(Q,3)*SIZE(Q,4)
			total_length_x = Ny*Nghost*nr_fields
			total_length_y = Nghost*Nx*nr_fields
	
			!$OMP CRITICAL
			CALL PostSends(Q(:,:,:,thread_nr), Nghost)
			CALL PostReceives(total_length_x,  total_length_y, 1)								
			CALL WaitForCommunication(Nx, Ny, Nghost, 1)
			!$OMP END CRITICAL
				
#elif (mpi_mode==3)			
			total_length_x = Ny*Nghost*nr_fields
			total_length_y = Nghost*Nx*nr_fields

			CALL PostSends(Q(:,:,:,thread_nr), Nghost)
			CALL PostReceives(total_length_x,  total_length_y, 1)								
			CALL WaitForCommunication(Nx, Ny, Nghost, 1)
#endif							
			
			! For nonlinear advection, simply put the index of the velocity field into the buffer otherwise containing the constant velocity field.
			CALL GetRHS( Q(:,:,:,thread_nr), order_coarse_advection, order_coarse_sound, RQ(:,:,:,thread_nr), dx_coarse, dy_coarse, dt, nu_coarse)
										
			CALL PackSolution( YDOT, RQ(:,:,:,thread_nr), nr_fields, Ny_coarse, Nx_coarse )

			IF (param%echo_on) WRITE(*,*) ' Finished RHS_coarse evaluation. '
			
			CONTAINS
			
				SUBROUTINE BARRIER_ONE_COARSE_PROFILING()
				
					INTEGER :: thread_nr 

					thread_nr = omp_get_thread_num()
					Timer_Temp(thread_nr) = MPI_WTIME()
					
					!$OMP BARRIER

					Timer_BOne(thread_nr)   = Timer_BOne(thread_nr) + MPI_WTIME() - Timer_Temp(thread_nr)
					NrCalls_BOne(thread_nr) = NrCalls_BOne(thread_nr) + 1

					
				END SUBROUTINE BARRIER_ONE_COARSE_PROFILING
				
				SUBROUTINE BARRIER_TWO_COARSE_PROFILING()
				
					INTEGER :: thread_nr

					thread_nr = omp_get_thread_num()
					Timer_Temp(thread_nr) = MPI_WTIME()
					
					!$OMP BARRIER
					
					Timer_BTwo(thread_nr)   = Timer_BTwo(thread_nr) + MPI_WTIME() - Timer_Temp(thread_nr)
					NrCalls_BTwo(thread_nr) = NrCalls_BTwo(thread_nr) + 1
					
				END SUBROUTINE BARRIER_TWO_COARSE_PROFILING
						
! Two important comments here:
! (i) The RESHAPE intrinsic routine is REALLY slow on the Cray, at least when compiled with PGI. The manual implementation with DO loops is magnitudes
!     faster
! (ii) This has NOT BEEN TESTED for the case where the coarse mesh is actually coarser than the fine mesh. For that case, leave the commented RESHAPE routines
!      in so that one can see what was working for this case.
					
		END SUBROUTINE RHS_coarse

END MODULE RHSFunctions
