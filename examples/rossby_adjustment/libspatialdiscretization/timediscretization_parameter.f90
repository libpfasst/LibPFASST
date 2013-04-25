MODULE timediscretization_parameter
	! This module reads and stores all parameter relevant for the time discretization. The nomenclature is as follows:
	!
	! Nsteps_fine_total       : the total number of fine steps from Tstart to Tend
	! Nsteps_coarse_total     : the total number of coarse steps from Tstart to Tend
	! Nsteps_coarse_per_slice : the number of coarse steps treated by one thread in Parareal
	! Nthreads                : the number of threads used by Parareal
	! global_tstart           : initial point in time
	! global_tend             : final point in time until which the problem is to be integrated
	!
	! Note that the following depencies must be fullfilled:
	!
	! Nparareal_restarts     = Nsteps_coarse/(Nsteps_coarse_per_slice*Nthreads) must be an integer 
	! Nsteps_fine_per_coarse = Nsteps_fine_total/Nsteps_fine_coarse must be an integer
	!
	! Based on these values, the interval [Tstart,Tend] is first divided into Nsteps_coarse many intervals of length dt_coarse = (Tend - Tstart)/Nsteps_coarse
	! 
	! Below is a sketch with
	!  Nsteps_fine_total       = 24 and
	!  Nsteps_coarse_total     = 8
	!  Nthreads                = 2 and
    !  Nsteps_coarse_per_slice = 2
	!
	! leading to
	!
	! Nparareal_restarts     = Nsteps_coarse_total/(Nsteps_coarse_per_slice*Nthreads) = 8/(2*2) = 2 (note that the first start of Parareal is also counted as a "restart")
	! Nsteps_fine_per_slice  = Nsteps_fine_total/(Nsteps_coarse_per_slice*Nthreads) = 24/(2*2) = 6  
	!
	!    <------ thread 0 --------><------ thread 1 --------> *restart* <------ thread 0 --------><------ thread 1 -------->
	!
	!    <===========> = dt_coarse
	! G: |-----------||-----------||-----------||-----------|           |-----------||-----------||-----------||-----------|
	!
	!    <===> = dt_fine
	! F: |---|---|---||---|---|---||---|---|---||---|---|---|           |---|---|---||---|---|---||---|---|---||---|---|---|
	!    *                                                  *           *                                                  *  
	!  Tstart                                              (Tend-Tstart)/2                                                Tend
	! 
	! Daniel Ruprecht
	! Institute of Computational Science
	! October 31, 2012

	USE serialio,     only : ReadIntegerDataset, ReadDoubleDataset, ActivateGroup, DeactivateGroup, WriteIntegerDataset, WriteDoubleDataset, InitializeInput, CloseInput
	USE omp_lib

	IMPLICIT NONE

	! Define solution variables
	DOUBLE PRECISION :: global_tstart, &
						global_tend,   &
						nu_coarse,     &
						nu,            &
						dt_coarse,     &
						dt_fine
													 
												   
	INTEGER, SAVE :: Nthreads, Nthreads_provided

	INTEGER, SAVE :: Nsteps_coarse_per_slice, &
			   Nparareal_restarts,      &
			   Nsteps_fine_per_coarse,  &
			   Ncores,                  &
			   Nsteps_sound_coarse,     &
			   Nsteps_sound_fine

	DOUBLE PRECISION :: tolit
	INTEGER			 :: maxit
	
	TYPE :: timedisc_parameter
		LOGICAL :: time_serial, echo_on, time_parallel_mpi
		INTEGER :: mpi_time_communicator, myrank_time, dim	
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: Tchunks
	END TYPE
	
	TYPE(timedisc_parameter) :: param
	
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: T_mesh
	
		CONTAINS
		
			SUBROUTINE ReadIntegrationParameter(time_serial, echo_on)
			
				LOGICAL,          INTENT(IN) :: time_serial, echo_on
							
				DOUBLE PRECISION :: dt_c_temp, dt_temp
				INTEGER          :: i, Nsteps_coarse_per_parareal, Nsteps_fine_total, Nsteps_coarse_total				
				
				param%time_serial = time_serial
				param%echo_on     = echo_on

				CALL InitializeInput(echo_on = echo_on, readonly = .true.)
				
				! Read all values from HDF5 group 'integration'
				CALL ActivateGroup('integration')
				CALL ReadDoubleDataset('global_tstart',            global_tstart)
				CALL ReadDoubleDataset('global_tend',              global_tend)
				CALL ReadDoubleDataset('nu_coarse',                nu_coarse)
				CALL ReadDoubleDataset('nu',                       nu)
				CALL ReadIntegerDataset('Nsteps_fine_total',       Nsteps_fine_total)
				CALL ReadIntegerDataset('Nsteps_coarse_total',     Nsteps_coarse_total)
				CALL ReadIntegerDataset('Nsteps_coarse_per_slice', Nsteps_coarse_per_slice)
				CALL ReadIntegerDataset('Nsteps_sound_coarse',     Nsteps_sound_coarse)
				CALL ReadIntegerDataset('Nsteps_sound_fine',       Nsteps_sound_fine)
				CALL DeactivateGroup()
				! Close HDF5 group integration
				
				
				! Read all values from HDF5 group 'parareal'
				CALL ActivateGroup('parareal')
				CALL ReadDoubleDataset('tolit',     tolit)
				CALL ReadIntegerDataset('maxit',    maxit)
				CALL ReadIntegerDataset('nthreads', Nthreads)			
				CALL DeactivateGroup()
				! Close HDF5 group 'parareal'
				CALL CloseInput()

				Ncores = omp_get_num_procs()
							
				Nsteps_coarse_per_parareal = Nthreads*Nsteps_coarse_per_slice
				
				! Verify that Nsteps_coarse_total is a multiple of Nsteps_coarse_per_parareal, so that division will yield an integer
				IF ( mod(Nsteps_coarse_total, Nsteps_coarse_per_parareal) .ne. 0 ) THEN
					WRITE(*,*) 'The total number of coarse steps must be a multiple of the coarse steps per each restart of parareal. Now exiting.'
					STOP
				END IF
				Nparareal_restarts = Nsteps_coarse_total/(Nsteps_coarse_per_parareal) 
								
				! Verify that Nsteps_fine_total is a multiple of Nsteps_coarse_total, so that division will yield an integer
				IF ( mod(Nsteps_fine_total, Nsteps_coarse_total) .ne. 0) THEN
					WRITE(*,*) 'The total number of fine steps must be a multiple of the total number of coarse steps. Now exiting'
					STOP
				END IF				
				Nsteps_fine_per_coarse = Nsteps_fine_total/Nsteps_coarse_total	
				
				! Compute timestep lengths
				dt_coarse = ( global_tend - global_tstart)/DBLE(Nsteps_coarse_total)
				dt_fine = ( global_tend - global_tstart )/DBLE(Nsteps_fine_total)
				
				ALLOCATE(T_mesh(Nparareal_restarts+1))
				
				DO i=1,Nparareal_restarts+1
					T_mesh(i) = DBLE(i-1)*DBLE(Nsteps_coarse_per_parareal)*dt_coarse
				END DO
				
				Nthreads_provided = omp_get_max_threads()
				! Finally, save the derived integration parameters into the group /input/derived/ so that they can
				! conviently be used during postprocessing
							
				IF (echo_on) WRITE(*, '(A, I2)' ) ' Integration parameters read. Number of provided threads : ', Nthreads_provided				
				
			END SUBROUTINE ReadIntegrationParameter
			
			SUBROUTINE WriteDerivedParameter(echo_on)
			
				LOGICAL, INTENT(IN) :: echo_on
			
				CALL InitializeInput(echo_on = echo_on, readonly = .false.)
			
				CALL ActivateGroup('derived')
				CALL WriteIntegerDataset('Nparareal_restarts',     Nparareal_restarts)
				CALL WriteIntegerDataset('Nsteps_fine_per_coarse', Nsteps_fine_per_coarse)
				CALL WriteIntegerDataset('Nthreads_provided',      Nthreads_provided)
				CALL WriteIntegerDataset('Ncores',                 Ncores)
				CALL WriteDoubleDataset('dt_coarse',               dt_coarse)
				CALL WriteDoubleDataset('dt_fine',                 dt_fine)
				CALL WriteDoubleDataset('dt_parallel',             T_mesh(2)-T_mesh(1))
				CALL DeactivateGroup()
				CALL CloseInput()
				IF (echo_on) WRITE(*,*) 'Successful write of derived parameter'
				
			END SUBROUTINE WriteDerivedParameter
			
			SUBROUTINE EchoIntegrationParameter(nprocs, u_max, dx, dx_coarse)
			
					INTEGER,                    INTENT(IN) :: nprocs
					DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: u_max, dx, dx_coarse
			
					WRITE (*,*)
					WRITE (*, '(A)') '********** Information from Temporal Mesh Generation **********'
					WRITE (*,*)
					WRITE (*, '(A, F8.5)') ' Starting time                            : ', T_mesh(1)
					WRITE (*, '(A, F8.5)') ' End time                                 : ', T_mesh(Nparareal_restarts+1)
					WRITE (*,*)
					WRITE (*, '(A, I7)') ' Number of parareal restarts              : ', Nparareal_restarts
					WRITE (*, '(A, I7)') ' Number of coarse steps per parareal run  : ', Nthreads*Nsteps_coarse_per_slice
					WRITE (*, '(A, I7)') ' Number of fine steps per parareal run    : ', Nthreads*Nsteps_coarse_per_slice*Nsteps_fine_per_coarse
					WRITE (*, '(A, I7)') ' Total number of coarse steps             : ', Nparareal_restarts*Nthreads*Nsteps_coarse_per_slice
					WRITE (*, '(A, I7)') ' Total number of fine steps               : ', Nparareal_restarts*Nthreads*Nsteps_coarse_per_slice*Nsteps_fine_per_coarse
					WRITE (*, '(A, F5.3)') ' Length of each parareal interval         : ', T_mesh(2) - T_mesh(1)
					WRITE (*, '(A, I7)') ' Number of threads per parallel step      : ', Nthreads
					WRITE (*, '(A, I7)') ' Number of coarse steps per thread        : ', Nsteps_coarse_per_slice
					WRITE (*, '(A, I7)') ' Number of fine steps per thread          : ', Nsteps_coarse_per_slice*Nsteps_fine_per_coarse
					WRITE (*,*)
					WRITE (*, '(A, F10.7)') ' Coarse time step length                  : ', dt_coarse
					WRITE (*, '(A, F10.7)') ' Fine time step length                    : ', dt_fine
					IF (PRESENT(dx_coarse) .and. PRESENT(u_max)) WRITE (*, '(A, F8.5)') ' Coarse CFL number                        : ', dt_coarse*u_max/dx_coarse
					IF (PRESENT(dx) .and. PRESENT(u_max)) WRITE (*, '(A, F8.5)') ' Fine CFL number                          : ', dt_fine*u_max/dx
					WRITE (*,*)
					WRITE (*, '(A, I2)') ' Number of iterations by Parareal         : ', maxit
					WRITE (*, '(A)') '***************************************************************'
					WRITE (*,*)
			
					IF (param%time_serial) THEN
						WRITE(*,*)
						WRITE(*, '(A)') ' ******************************'
						WRITE(*, '(A, I7)') ' Used MPI processes : ', nprocs
						WRITE(*, '(A)') ' --> No multithreading, running sequentially in time '
						WRITE(*, '(A)') ' ******************************'
						WRITE(*,*)					
					ELSE
						WRITE(*,*)
						WRITE(*, '(A)') ' ******************************'
						WRITE(*, '(A, I7)') ' Used MPI processes : ', nprocs
						WRITE(*, '(A, I7)') ' Requested threads  : ', Nthreads
						WRITE(*, '(A, I7)') ' Provided threads   : ', Nthreads_provided
						WRITE(*, '(A, I7)') ' Available cores    : ', Ncores
						WRITE(*, '(A)') ' ******************************'
						WRITE(*,*)
					END IF
					
			END SUBROUTINE EchoIntegrationParameter
					
END MODULE timediscretization_parameter