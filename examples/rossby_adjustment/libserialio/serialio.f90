MODULE serialio
	! Simple module that provides routines that allow to read integer and double precision data sets out of a HDF5 file.
	!
	! Daniel Ruprecht
	! Institute of Computational Science, Lugano
	! October 18, 2012 
	
	USE HDF5

	IMPLICIT NONE
	
	CHARACTER(len=32) :: filename	
	
	TYPE serial_io_parameter
		INTEGER(HID_T) :: group_id_input, file_id
		INTEGER(HID_T) :: group_id_active
		LOGICAL        :: echo_on
	END TYPE
	
	TYPE(serial_io_parameter) :: param
	
	PRIVATE
	PUBLIC :: InitializeInput, CloseInput, ReadIntegerDataset, ReadDoubleDataset, WriteIntegerDataset, WriteDoubleDataset, ActivateGroup, DeactivateGroup
												
	CONTAINS
	
		! This routine open the HDF5 data file that contains all required parameters for the simulation. It genererates a file ID, that can then be used
		! by the parameter modules to read the values from the corresponding HDF5 groups.
		SUBROUTINE InitializeInput(echo_on, readonly)
		
			LOGICAL,        INTENT(IN)  :: echo_on, &         ! if true, successful access of the file is reported
										   readonly
										   	
			INTEGER(HID_T) :: proplist_file
			INTEGER        :: hdf5_error
		
			! Start the HDF5 library
			CALL H5OPEN_F(hdf5_error)

			param%echo_on = echo_on
					
			! Generate property list for file access
			CALL H5PCREATE_F(H5P_FILE_ACCESS_F, proplist_file, hdf5_error)
			
			! File to be used is specified as first command line argument
			CALL getarg(1, filename)

			! Open file and retrieve a file ID
			IF (readonly) THEN
				CALL H5FOPEN_F(filename, H5F_ACC_RDONLY_F, param%file_id, hdf5_error, access_prp = proplist_file)
			ELSE
				CALL H5FOPEN_F(filename, H5F_ACC_RDWR_F, param%file_id, hdf5_error, access_prp = proplist_file)
			END IF
			
			IF (hdf5_error >= 0) THEN
				IF (param%echo_on) WRITE(*,*) ' File successfully opened, now commencing to read parameters .... '
			ELSE
					WRITE(*,*) ' ERROR OPENING FILE ! '
			END IF
		
			! Close property lists
			CALL H5PCLOSE_F(proplist_file, hdf5_error)
			
			CALL H5GOPEN_F(param%file_id, 'input', param%group_id_input, hdf5_error)
								
		END SUBROUTINE InitializeInput
			
		! close file for input
		SUBROUTINE CloseInput()
					
			INTEGER :: hdf5_error

			! Close groups
			CALL H5GCLOSE_F(param%group_id_input, hdf5_error)

			! Close file
			CALL H5FCLOSE_F(param%file_id, hdf5_error)

			! Close HDF5
			CALL H5CLOSE_F(hdf5_error)
						
			IF (hdf5_error >= 0) THEN
				IF (param%echo_on) WRITE(*,*) 'SUCCESS closing input module'
			ELSE
				WRITE(*,*) 'Negative error flag in H5_CLOSE, something is wrong and HDF5 could not be closed properly'
			END IF
			
		END SUBROUTINE CloseInput
	
		SUBROUTINE ActivateGroup(group_name)
			CHARACTER(len=*), INTENT(IN) :: group_name
			INTEGER :: hdf5_error
			LOGICAL :: group_exists
			CALL H5LEXISTS_F(param%group_id_input, group_name, group_exists, hdf5_error)
			
			IF (.not.(group_exists)) THEN
				IF (param%echo_on) WRITE(*,*) ' Creating group : ', group_name
				CALL H5GCREATE_F(param%group_id_input, group_name, param%group_id_active, hdf5_error)			
			END IF
			CALL H5GOPEN_F(param%group_id_input, group_name, param%group_id_active, hdf5_error)
			
		END SUBROUTINE ActivateGroup
		
		SUBROUTINE DeactivateGroup()
			INTEGER :: hdf5_error
					CALL H5GCLOSE_F(param%group_id_active, hdf5_error)
		END SUBROUTINE DeactivateGroup
		
		SUBROUTINE ReadIntegerDataset(dataset_name, value)
		
			CHARACTER(len=*), INTENT(IN)  :: dataset_name ! name of the HDF5 dataseto to be read
			INTEGER,          INTENT(OUT) :: value        ! the data read from the file
			
			INTEGER(HID_T)   :: dataset_id
			INTEGER(HSIZE_T) :: dims(1) = (/ 1 /)
			INTEGER          :: hdf5_error
			
			CALL H5DOPEN_F(param%group_id_active, dataset_name, dataset_id, hdf5_error)
			CALL H5DREAD_F(dataset_id, H5T_NATIVE_INTEGER, value, dims, hdf5_error)
			CALL H5DCLOSE_F(dataset_id, hdf5_error)					 
									 
		END SUBROUTINE ReadIntegerDataset
		
		SUBROUTINE ReadDoubleDataset(dataset_name, value)
		
			CHARACTER(len=*), INTENT(IN) :: dataset_name ! name of the HDF5 dataseto to be read
			DOUBLE PRECISION, INTENT(OUT) :: value       ! the data read from the file
			
			INTEGER(HID_T)   :: dataset_id
			INTEGER(HSIZE_T) :: dims(1) = (/ 1 /)
			INTEGER          :: hdf5_error
			
			CALL H5DOPEN_F(param%group_id_active, dataset_name, dataset_id, hdf5_error)
			CALL H5DREAD_F(dataset_id, H5T_NATIVE_DOUBLE, value, dims, hdf5_error)
			CALL H5DCLOSE_F(dataset_id, hdf5_error)					 
									 
		END SUBROUTINE ReadDoubleDataset
		
		SUBROUTINE WriteIntegerDataset(dataset_name, value)
		
			CHARACTER(len=*), INTENT(IN) :: dataset_name
			INTEGER,          INTENT(IN) :: value		
		
			INTEGER(HID_T) :: dataspace, dataset, hdf5_error
			INTEGER(HSIZE_T) :: dims(1)
			
			dims = (/ 1 /)
			CALL H5SCREATE_SIMPLE_F( INT(1), dims, dataspace, hdf5_error)
					
			CALL H5DCREATE_F(param%group_id_active, dataset_name, H5T_NATIVE_INTEGER, dataspace, dataset, hdf5_error)
			CALL H5DWRITE_F(dataset, H5T_NATIVE_INTEGER, value, dims, hdf5_error)
			
			CALL H5DCLOSE_F(dataset, hdf5_error)
			CALL H5SCLOSE_F(dataspace, hdf5_error)
			
		END SUBROUTINE WriteIntegerDataset

		SUBROUTINE WriteDoubleDataset(dataset_name,  value)
		
			CHARACTER(len=*), INTENT(IN) :: dataset_name
			DOUBLE PRECISION, INTENT(IN) :: value		
		
			INTEGER(HID_T) :: dataspace, dataset, hdf5_error
			INTEGER(HSIZE_T) :: dims(1)
			
			dims = (/ 1 /)
			CALL H5SCREATE_SIMPLE_F( INT(1), dims, dataspace, hdf5_error)
					
			CALL H5DCREATE_F(param%group_id_active, dataset_name, H5T_NATIVE_DOUBLE, dataspace, dataset, hdf5_error)
			CALL H5DWRITE_F(dataset, H5T_NATIVE_DOUBLE, value, dims, hdf5_error)
			
			CALL H5DCLOSE_F(dataset, hdf5_error)
			CALL H5SCLOSE_F(dataspace, hdf5_error)
			
		END SUBROUTINE WriteDoubleDataset		
		
END MODULE serialio