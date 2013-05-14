MODULE DistributedIO

	USE HDF5
	USE FVMParameters,         only : Nx, Ny, nr_fields
	USE MPIParameter,          only : cartesian_comm, MPI_INFO_NULL, nprocs, nprocs_x, nprocs_y, cart_coords, myrank
        USE FiniteVolumes,         only : PackSolution, buffer_layout

	IMPLICIT NONE
		
	CHARACTER(len=32) :: filename

	INTEGER :: hdf5_error
	
	TYPE :: distributedio_parameter
		INTEGER :: maxit, Nsteps_parallel, Nthreads
		LOGICAL :: echo_on
	END TYPE
	
	TYPE(distributedio_parameter) :: param
	
	! Define dimensionality of different objects to be used during I/O
	! keywords used in names : 
	!
	! "global"   : contains values from all processes
	! "local"    : contains values from only the current process
	! "constant" : a non-timedependent quantity, no temporal dimension
	! "common"   : a quantity that is identical for all processes, so that only one value is stored (e.g. timestep size)
	!
	INTEGER(HSIZE_T) :: dims_global_solution(4),        &! Full solution over time, target I/O buffer           4 dimensions: x, z, solution component, time
						dims_local_solution(3),         &! Full solution on one processor at one point in time, 3 dimsions:   x, z, solution component
						dims_global_scalar(2),          &! Process dependent scalar at one point in time,       2 dimensions: processor rank, time
						dims_local_scalar(1),           &! Scalar on one processor at one point in tme,         single value
						dims_global_iterate_scalar(3),  &! Scalar value produced by every processor for every Parareal iteration, 3 dimensions: processor rank, iteration, time
						dims_local_iterate_scalar(1),   &! Scalar value on one processor, one value for every Parareal iteration, 1 dimensions: iteration
						dims_global_constant_scalar(1), &! Process dependent, time-independent scalar,        1 dimension: processor rank
						dims_local_constant_scalar(1),  &! Process dependent, time-independent scalar on one processor,   single value
						dims_global_threads(2),         &! Process dependent, nr_thread many entries per process, 2 dimensions: processor rank, thread number
						dims_local_threads(1)            ! Single process variable with one entry per thread, 1 dimension: thread number 
	
	! Define names for property lists
	INTEGER(HID_T) :: proplist_file,     &!
					  proplist_transfer, &!
					  file_id,           &!
					  group_output_id,   &!
					  group_timer_id
					  	
	! Define names of data spaces for objects
	INTEGER(HID_T) :: dataspace_global_solution, &!
					  dataspace_global_scalar,   &!
					  dataspace_global_iterate_scalar,     &!
					  dataspace_global_constant_scalar,    &!
					  dataspace_global_threads,            &!
					  memspace_local_solution,             &!
					  memspace_local_scalar,               &!
					  memspace_local_iterate_scalar,       &!
					  memspace_local_threads
					  										  	  
	! Define the number of variables passed to the subroutines "WriteProblemData", "WriteData" and "WriteFinalData". Required
	! to allocate the correct number of IDs in the buffers defined below
	INTEGER, PARAMETER :: nr_vars_output = 11
											  
	! Define buffers for IDs of data sets. As the number of variables to be written can be specified
	INTEGER(HID_T) :: dataset_ids_data(nr_vars_output)
					  
	! NOTE: The mappings used for writing data are 
	!   (i) memspace_local_solution          ---> (Ny x Nx x nr_fields) hyperslab in dataspace_timeslice_global_solution at specified time point index
	!  (ii) memspace_local_scalar            --->                   (1) hyperslab in dataspace_global_scalar
	! (iii) memspace_local_iterate_scalar    --->               (maxit) hyperslab in dataspace_global_iterate_scalar
	!  (iv) dataspace_memspace_common_scalar ---> dataspace_memspace_common_scalar
	
	PRIVATE
	PUBLIC :: InitializeDistributedIO, ReadInitialValue, WriteData, WriteTimerData, WriteFinalData, CloseDistributedIO, ReadLinearAdvectionVelocity
	
	CONTAINS
	

		! This routine creates a group "output" in the HDF5 data file used for the simulation. It also creates the required data spaces and datasets for
		! writing data.
		SUBROUTINE InitializeDistributedIO(maxit, Nsteps_parallel, Nthreads, echo_on)
		
			INTEGER, INTENT(IN) :: maxit, Nsteps_parallel, Nthreads
			LOGICAL, INTENT(IN) :: echo_on
		
			param%maxit = maxit
			param%Nsteps_parallel = Nsteps_parallel
			param%Nthreads = Nthreads
			param%echo_on = echo_on
		
			! ------------------------------------------------------------------------
			! STEP 1 : Open the HDF5 file again and prepare for collective I/O, generate output file
		
			CALL getarg(1, filename)

			! Start the HDF5 library
			CALL H5OPEN_F(hdf5_error)
					
			! Generate property list for file access
			CALL H5PCREATE_F(H5P_FILE_ACCESS_F, proplist_file, hdf5_error)			

			! Provide MPI communicator to file access property list. Note: IO should only be performed by the master thread
			CALL H5PSET_FAPL_MPIO_F(proplist_file, cartesian_comm(0), MPI_INFO_NULL, hdf5_error)
			
			! Create property list for data transfer
			CALL H5PCREATE_F(H5P_DATASET_XFER_F, proplist_transfer, hdf5_error)
			
			! Start up file transfer for collective I/O
			CALL H5PSET_DXPL_MPIO_F(proplist_transfer, H5FD_MPIO_COLLECTIVE_F, hdf5_error)

			CALL H5FOPEN_F(filename, H5F_ACC_RDWR_F, file_id, hdf5_error, access_prp = proplist_file)

			! Create the output group
			CALL H5GCREATE_F(file_id, 'output', group_output_id, hdf5_error)
			
			! Create the timer group
			CALL H5GCREATE_F(file_id, 'timer', group_timer_id, hdf5_error)
			
			! ---------------------------------------------------------------------------------------
			! STEP 2 : Define all data space (that is, dimensionalities) that will be used during I/O
			
			! a) Define global data spaces, containing data from all processors. These are used as targets for output
			!dims_global_solution        = (/ Ny*nprocs_y, Nx*nprocs_x, nr_fields, Nsteps_parallel /)
                        dims_global_solution        = (/ nr_fields, Ny*nprocs_y, Nx*nprocs_x, Nsteps_parallel /)
			dims_global_iterate_scalar  = (/ nprocs, maxit, Nsteps_parallel /)
			dims_global_scalar          = (/ nprocs, Nsteps_parallel /)
			dims_global_constant_scalar = (/ nprocs /)
			dims_global_threads         = (/ nprocs, Nthreads /)

			CALL H5SCREATE_SIMPLE_F( INT(4), dims_global_solution,        dataspace_global_solution,        hdf5_error)			
			CALL H5SCREATE_SIMPLE_F( INT(3), dims_global_iterate_scalar,  dataspace_global_iterate_scalar,  hdf5_error)			
			CALL H5SCREATE_SIMPLE_F( INT(2), dims_global_scalar,          dataspace_global_scalar,          hdf5_error)						
			CALL H5SCREATE_SIMPLE_F( INT(1), dims_global_constant_scalar, dataspace_global_constant_scalar, hdf5_error)
			CALL H5SCREATE_SIMPLE_F( INT(2), dims_global_threads,         dataspace_global_threads,         hdf5_error)
			
			! b) Define local data spaces, containing data from one processor. These are used as sources for output
			!dims_local_solution       = (/ Ny, Nx, nr_fields /)
			dims_local_solution       = (/ nr_fields, Ny, Nx /)
                        dims_local_iterate_scalar = (/ maxit /)
			dims_local_scalar         = (/ 1 /)
			dims_local_threads        = (/ Nthreads /)
			
			CALL H5SCREATE_SIMPLE_F( INT(3), dims_local_solution,       memspace_local_solution,       hdf5_error)			
			CALL H5SCREATE_SIMPLE_F( INT(1), dims_local_iterate_scalar, memspace_local_iterate_scalar, hdf5_error)
			CALL H5SCREATE_SIMPLE_F( INT(1), dims_local_scalar,         memspace_local_scalar,         hdf5_error)			
			CALL H5SCREATE_SIMPLE_F( INT(1), dims_local_threads,        memspace_local_threads,        hdf5_error)
			
			! ------------------------------------------------------------------------
			! STEP 3: Create data sets contained in the group "output"
			CALL H5DCREATE_F(group_output_id, 'solution',       H5T_NATIVE_DOUBLE,  dataspace_global_solution,        dataset_ids_data(1), hdf5_error)
			CALL H5DCREATE_F(group_output_id, 'residuals',      H5T_NATIVE_DOUBLE,  dataspace_global_iterate_scalar,  dataset_ids_data(2), hdf5_error)
			CALL H5DCREATE_F(group_output_id, 'runtime_coarse', H5T_NATIVE_DOUBLE,  dataspace_global_scalar,          dataset_ids_data(3), hdf5_error)
			CALL H5DCREATE_F(group_output_id, 'runtime_fine',   H5T_NATIVE_DOUBLE,  dataspace_global_scalar,          dataset_ids_data(4), hdf5_error)
			CALL H5DCREATE_F(group_output_id, 'runtime_qr',     H5T_NATIVE_DOUBLE,  dataspace_global_scalar,          dataset_ids_data(5), hdf5_error)
			CALL H5DCREATE_F(group_output_id, 'iterations',     H5T_NATIVE_INTEGER, dataspace_global_scalar,          dataset_ids_data(6), hdf5_error)
			CALL H5DCREATE_F(group_output_id, 'total_runtime',  H5T_NATIVE_DOUBLE,  dataspace_global_constant_scalar, dataset_ids_data(7), hdf5_error)			
			
			! ------------------------------------------------------------------------
			! STEP 4: Create data sets contained in the group "timer"
			CALL H5DCREATE_F(group_timer_id, 'Timer_BarrierOne', H5T_NATIVE_DOUBLE, dataspace_global_threads,         dataset_ids_data(8), hdf5_error)			
			CALL H5DCREATE_F(group_timer_id, 'Timer_BarrierTwo', H5T_NATIVE_DOUBLE, dataspace_global_threads,         dataset_ids_data(9), hdf5_error)
			CALL H5DCREATE_F(group_timer_id, 'NrCalls_BarrierOne', H5T_NATIVE_INTEGER, dataspace_global_threads,      dataset_ids_data(10), hdf5_error)
			CALL H5DCREATE_F(group_timer_id, 'NrCalls_BarrierTwo', H5T_NATIVE_INTEGER, dataspace_global_threads,      dataset_ids_data(11), hdf5_error)
			
			IF (echo_on) WRITE(*,'(A, I4)') ' HDF5 succesfully initalized for I/O. Number of dumps : ', Nsteps_parallel

		END SUBROUTINE InitializeDistributedIO
		
		SUBROUTINE ReadLinearAdvectionVelocity(Uadv, Vadv)

                        DOUBLE PRECISION, DIMENSION(:,:),   INTENT(OUT) :: Uadv, Vadv

			INTEGER(HSIZE_T)  :: dims_local(3), array_adv_count(2), dims_adv(2), dims_adv_local(2)
			INTEGER(HSSIZE_T) :: array_adv_offset(2)
			INTEGER(HID_T)    :: group_input_id, group_problem_id, dataset, dataspace, dataspace_local

			CALL H5GOPEN_F(file_id,        'input',             group_input_id,   hdf5_error)
			CALL H5GOPEN_F(group_input_id, 'problemdefinition', group_problem_id, hdf5_error)
			
		! ---- Load Uadv -----
			array_adv_offset = (/ cart_coords(1)*Ny , cart_coords(2)*Nx /)
			array_adv_count  = (/ Ny , Nx+1  /)
			
			CALL H5DOPEN_F(group_problem_id, 'Uadv', dataset, hdf5_error)
			
			dims_adv = (/ nprocs_y*Ny, nprocs_x*Nx+1 /)
			CALL H5SCREATE_SIMPLE_F( INT(2), dims_adv, dataspace, hdf5_error)
			CALL H5SSELECT_HYPERSLAB_F(dataspace, H5S_SELECT_SET_F, array_adv_offset, array_adv_count, hdf5_error)
			
			dims_adv_local = (/ Ny, Nx+1 /)
			CALL H5SCREATE_SIMPLE_F( INT(2), dims_adv_local, dataspace_local, hdf5_error)
			
			CALL H5DREAD_F(dataset, H5T_NATIVE_DOUBLE, Uadv, dims_local, hdf5_error, &
				mem_space_id  = dataspace_local, &
				file_space_id = dataspace,       &
				xfer_prp      = proplist_transfer)
			
			CALL H5DCLOSE_F(dataset,         hdf5_error)
			CALL H5SCLOSE_F(dataspace,       hdf5_error)
			CALL H5SCLOSE_F(dataspace_local, hdf5_error)
			
			! ---- Load Vadv ----
			array_adv_offset = (/ cart_coords(1)*Ny , cart_coords(2)*Nx /)
			array_adv_count	 = (/ Ny+1 , Nx  /)
			
			CALL H5DOPEN_F(group_problem_id, 'Vadv', dataset, hdf5_error)
			
			dims_adv = (/ nprocs_y*Ny+1, nprocs_x*Nx /)
			CALL H5SCREATE_SIMPLE_F( INT(2), dims_adv, dataspace, hdf5_error)
			CALL H5SSELECT_HYPERSLAB_F(dataspace, H5S_SELECT_SET_F, array_adv_offset, array_adv_count, hdf5_error)
			
			dims_adv_local = (/ Ny+1, Nx /)
			CALL H5SCREATE_SIMPLE_F( INT(2), dims_adv_local, dataspace_local, hdf5_error)
			
			CALL H5DREAD_F(dataset, H5T_NATIVE_DOUBLE, Vadv, dims_adv_local, hdf5_error, &
				mem_space_id  = dataspace_local, &
				file_space_id = dataspace, &
				xfer_prp      = proplist_transfer)
			
			CALL H5DCLOSE_F(dataset,         hdf5_error)
			CALL H5SCLOSE_F(dataspace,       hdf5_error)
			CALL H5SCLOSE_F(dataspace_local, hdf5_error)
						
			! -- close groups --
			CALL H5GCLOSE_F(group_problem_id, hdf5_error)
			CALL H5GCLOSE_F(group_input_id,   hdf5_error)

			
		END SUBROUTINE ReadLinearAdvectionVelocity
		
		SUBROUTINE ReadInitialValue(Y0)
		
			DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Y0
			
			DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q_initial
			INTEGER(HSIZE_T)  :: array_count(3), dims(3), dims_local(3)
			INTEGER(HSSIZE_T) :: array_offset(3)
			INTEGER(HID_T)    :: group_input_id, group_problem_id, dataset, dataspace, dataspace_local
			
			!ALLOCATE(Q_initial(Ny, Nx, nr_fields))
			ALLOCATE(Q_initial(nr_fields, Ny, Nx))

			CALL H5GOPEN_F(file_id,        'input',             group_input_id,   hdf5_error)
			CALL H5GOPEN_F(group_input_id, 'problemdefinition', group_problem_id, hdf5_error)
												
			! ---- Load initial value ----	
                        ! array_offset = (/ cart_coords(1)*Ny , cart_coords(2)*Nx, 0 /)
                        ! array_count  = (/ Ny , Nx, nr_fields /) 
                        array_offset = (/ 0, cart_coords(1)*Ny, cart_coords(2)*Nx /)
                        array_count  = (/ nr_fields, Ny, Nx /)  

			CALL H5DOPEN_F(group_problem_id, 'q_initial', dataset, hdf5_error)
		
			!dims = (/ nprocs_y*Ny , nprocs_x*Nx, nr_fields /)
			dims = (/ nr_fields, nprocs_y*Ny, nprocs_x*Nx /)
                        CALL H5SCREATE_SIMPLE_F( INT(3), dims, dataspace, hdf5_error)			
			CALL H5SSELECT_HYPERSLAB_F(dataspace, H5S_SELECT_SET_F, array_offset, array_count, hdf5_error)
			
			!dims_local = (/ Ny, Nx, nr_fields /)
			dims_local = (/ nr_fields, Ny, Nx /)
                        CALL H5DREAD_F(dataset, H5T_NATIVE_DOUBLE, Q_initial, dims_local, hdf5_error, &
				mem_space_id  = memspace_local_solution,  &
				file_space_id = dataspace,                &
				xfer_prp      = proplist_transfer)
			
			CALL H5DCLOSE_F(dataset,          hdf5_error)
			CALL H5SCLOSE_F(dataspace,        hdf5_error)

			! -- close groups --
			CALL H5GCLOSE_F(group_problem_id, hdf5_error)
			CALL H5GCLOSE_F(group_input_id,   hdf5_error)
			
			CALL PackSolution(Y0, Q_initial, nr_fields, Ny, Nx)

			DEALLOCATE(Q_initial)
			
		END SUBROUTINE ReadInitialValue
		
		SUBROUTINE WriteData(residuals, runtime_coarse, runtime_fine, runtime_qr, iterations, time_index, Q)
			
			DOUBLE PRECISION, DIMENSION(:),               INTENT(IN) :: residuals
			DOUBLE PRECISION,                             INTENT(IN) :: runtime_coarse, runtime_fine, runtime_qr
			INTEGER,                                      INTENT(IN) :: iterations, time_index
			DOUBLE PRECISION, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: Q
			
			INTEGER(HSIZE_T)  :: array_count(4),  scalar_count(2),  scalar_iterate_count(3)
			INTEGER(HSSIZE_T) :: array_offset(4), scalar_offset(2), scalar_iterate_offset(3)
			
			! Write array solution

			IF (PRESENT(Q)) THEN
				!array_offset = (/ cart_coords(1)*Ny , cart_coords(2)*Nx, 0, time_index-1 /)
				!array_count  = (/ Ny, Nx, nr_fields, 1 /)
				
				array_offset = (/ 0, cart_coords(1)*Ny , cart_coords(2)*Nx, time_index-1 /)
				array_count  = (/ nr_fields, Ny, Nx, 1 /)

                                CALL H5SSELECT_HYPERSLAB_F(dataspace_global_solution, H5S_SELECT_SET_F, array_offset, array_count, hdf5_error)

				CALL H5DWRITE_F(dataset_ids_data(1), H5T_NATIVE_DOUBLE, Q, dims_local_solution, hdf5_error, &
					file_space_id = dataspace_global_solution, &
					mem_space_id  = memspace_local_solution,   &
					xfer_prp      = proplist_transfer)

                                IF (hdf5_error .ne. 0) WRITE(*,*) 'Error in writing "Q"'

			END IF
			

			! Write array over iterations	
			scalar_iterate_offset = (/ myrank, 0, time_index-1 /)
			scalar_iterate_count  = (/ 1, iterations, 1 /)
			CALL H5SSELECT_HYPERSLAB_F(dataspace_global_iterate_scalar, H5S_SELECT_SET_F, scalar_iterate_offset, scalar_iterate_count, hdf5_error)
			
			CALL H5DWRITE_F(dataset_ids_data(2), H5T_NATIVE_DOUBLE, residuals, dims_local_iterate_scalar, hdf5_error, &
				file_space_id = dataspace_global_iterate_scalar, &
				mem_space_id  = memspace_local_iterate_scalar,   &
				xfer_prp      = proplist_transfer)

                        IF (hdf5_error .ne. 0) WRITE(*,*) 'Error in writing "residuals"'

			! Write scalar values	
			scalar_offset = (/ myrank, time_index-1 /)
			scalar_count  = (/ 1 , 1 /) 
			CALL H5SSELECT_HYPERSLAB_F(dataspace_global_scalar, H5S_SELECT_SET_F, scalar_offset, scalar_count, hdf5_error)

			CALL H5DWRITE_F(dataset_ids_data(3), H5T_NATIVE_DOUBLE, runtime_coarse, dims_local_scalar, hdf5_error, &
				file_space_id = dataspace_global_scalar, &
				mem_space_id  = memspace_local_scalar,   &
				xfer_prp      = proplist_transfer)
				
                        IF (hdf5_error .ne. 0) WRITE(*,*) 'Error in writing "runtime_coarse"'

			CALL H5DWRITE_F(dataset_ids_data(4), H5T_NATIVE_DOUBLE, runtime_fine, dims_local_scalar, hdf5_error, &
				file_space_id = dataspace_global_scalar, &
				mem_space_id  = memspace_local_scalar,   &
				xfer_prp      = proplist_transfer)

                        IF (hdf5_error .ne. 0) WRITE(*,*) 'Error in writing "runtime_fine"'
				
			CALL H5DWRITE_F(dataset_ids_data(5), H5T_NATIVE_DOUBLE, runtime_qr, dims_local_scalar, hdf5_error, &
				file_space_id = dataspace_global_scalar, &
				mem_space_id  = memspace_local_scalar,   &
				xfer_prp      = proplist_transfer)			

                        IF (hdf5_error .ne. 0) WRITE(*,*) 'Error in writing "runtime_qr"'																										
			CALL H5DWRITE_F(dataset_ids_data(6), H5T_NATIVE_INTEGER, iterations, dims_local_scalar, hdf5_error, &
				file_space_id = dataspace_global_scalar, &
				mem_space_id  = memspace_local_scalar,   &
				xfer_prp      = proplist_transfer)																			
                        IF (hdf5_error .ne. 0) WRITE(*,*) 'Error in writing "iterations"'																																				
		END SUBROUTINE WriteData
		
		SUBROUTINE WriteTimerData(Timer_B1, Timer_B2, NrCalls_B1, NrCalls_B2)
		
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Timer_B1, Timer_B2
			INTEGER,          DIMENSION(:), INTENT(IN) :: NrCalls_B1, NrCalls_B2
			
			INTEGER(HSIZE_T)  :: count(2)
			INTEGER(HSSIZE_T) :: offset(2)
			
			offset = (/ myrank, 0 /)
			count  = (/ 1, param%Nthreads /)
			
			CALL H5SSELECT_HYPERSLAB_F(dataspace_global_threads, H5S_SELECT_SET_F, offset, count, hdf5_error)
			
			CALL H5DWRITE_F(dataset_ids_data(8), H5T_NATIVE_DOUBLE, Timer_B1, dims_local_threads, hdf5_error, &
				file_space_id = dataspace_global_threads, &
				mem_space_id  = memspace_local_threads,   &
				xfer_prp      = proplist_transfer)
				
			CALL H5DWRITE_F(dataset_ids_data(9), H5T_NATIVE_DOUBLE, Timer_B2, dims_local_threads, hdf5_error, &
				file_space_id = dataspace_global_threads, &
				mem_space_id  = memspace_local_threads,   &
				xfer_prp      = proplist_transfer)				
			
			CALL H5DWRITE_F(dataset_ids_data(10), H5T_NATIVE_INTEGER, NrCalls_B1, dims_local_threads, hdf5_error, &
				file_space_id = dataspace_global_threads, &
				mem_space_id  = memspace_local_threads,   &
				xfer_prp      = proplist_transfer)	
				
			CALL H5DWRITE_F(dataset_ids_data(11), H5T_NATIVE_INTEGER, NrCalls_B2, dims_local_threads, hdf5_error, &
				file_space_id = dataspace_global_threads, &
				mem_space_id  = memspace_local_threads,   &
				xfer_prp      = proplist_transfer)	
														
		END SUBROUTINE WriteTimerData
		
		! 
		SUBROUTINE WriteFinalData(total_runtime)
			
			DOUBLE PRECISION, INTENT(IN) :: total_runtime
			
			INTEGER(HSIZE_T)  :: scalar_constant_count(1)
			INTEGER(HSSIZE_T) :: scalar_constant_offset(1)
			
			scalar_constant_offset = (/ myrank /)
			scalar_constant_count  = (/ 1 /)
			
			CALL H5SSELECT_HYPERSLAB_F(dataspace_global_constant_scalar, H5S_SELECT_SET_F, scalar_constant_offset, scalar_constant_count, hdf5_error)
			
			CALL H5DWRITE_F(dataset_ids_data(7), H5T_NATIVE_DOUBLE, total_runtime, dims_local_scalar, hdf5_error, &
				file_space_id = dataspace_global_constant_scalar, &
				mem_space_id  = memspace_local_scalar, &
				xfer_prp      = proplist_transfer)

		END SUBROUTINE WriteFinalData
		
		! Close output file
		SUBROUTINE CloseDistributedIO()
				
			INTEGER :: i

			! Close data sets used multiple times for I/O in "WriteData"
			DO i=1,nr_vars_output
                                CALL H5DCLOSE_F(dataset_ids_data(i), hdf5_error)
			END DO
		
			! NOTE: All other data sets should be closed by now, as closing of data sets is done by the "write" routines
		
			! Close generated data spaces
			CALL H5SCLOSE_F(dataspace_global_solution, hdf5_error)
			CALL H5SCLOSE_F(dataspace_global_scalar, hdf5_error)
			CALL H5SCLOSE_F(dataspace_global_constant_scalar, hdf5_error)
			CALL H5SCLOSE_F(dataspace_global_iterate_scalar, hdf5_error)
			CALL H5SCLOSE_F(dataspace_global_threads, hdf5_error)
					
			CALL H5SCLOSE_F(memspace_local_solution, hdf5_error)
			CALL H5SCLOSE_F(memspace_local_scalar, hdf5_error)
			CALL H5SCLOSE_F(memspace_local_iterate_scalar, hdf5_error)
			CALL H5SCLOSE_F(memspace_local_threads, hdf5_error)
													
			! Close property lists
			CALL H5PCLOSE_F(proplist_transfer, hdf5_error)
			CALL H5PCLOSE_F(proplist_file, hdf5_error)
			
			! Close groups
			CALL H5GCLOSE_F(group_output_id, hdf5_error)
			CALL H5GCLOSE_F(group_timer_id, hdf5_error)
			
			! Close file
			CALL H5FCLOSE_F(file_id, hdf5_error)
			
			! Close HDF5
			CALL H5CLOSE_F(hdf5_error)
						
			IF (hdf5_error >= 0) THEN
				IF (param%echo_on) WRITE(*,*) 'SUCCESS generating output'
			ELSE
				WRITE(*,*) 'Negative error flag in H5_CLOSE, something is wrong and HDF5 could not be closed properly'
			END IF
			
		END SUBROUTINE CloseDistributedIO
	
END MODULE DistributedIO
