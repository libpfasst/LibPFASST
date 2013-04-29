MODULE FVMParameters

USE serialio,      only : ReadIntegerDataset, ReadDoubleDataset, ActivateGroup, DeactivateGroup, InitializeInput, CloseInput
USE FiniteVolumes, only : nr_fields

IMPLICIT NONE

INTEGER :: order_coarse_advection, &
		   order_coarse_sound    , &
		   order_fine_advection  , &
		   order_fine_sound      , &
		   Nx					 , &
		   Ny					 , &
		   Nx_coarse             , &
		   Ny_coarse             , &
		   dim                   , &
		   dim_coarse            , &
		   nprocs_x              , &
		   nprocs_y              , &
		   BC                    
		   
DOUBLE PRECISION :: dx,            &
				    dy,            &
					dx_coarse,     &
					dy_coarse,     &
					c_s,           &
					stabFreq,      &
                                        coriolisPar,   &
                                        grav,          &
					domain_height, &
					domain_width,  &
					nu,            &
					nu_coarse
					
CONTAINS
		
	! IMPORTANT: Before invoking this subroutine, the "InitializeInput" subroutine in the input module
	! has to be called, in order to provided meaningful values for	
	SUBROUTINE ReadFVMParameter(echo_on, dim)

			LOGICAL, INTENT(IN) :: echo_on
			INTEGER, INTENT(OUT) :: dim
			
			CALL InitializeInput(echo_on = echo_on, readonly = .true.)

			
			! STEP ONE: Open group "discretization" in group "input" and read the parameters contained there
			
			CALL ActivateGroup('discretization')
			CALL ReadIntegerDataset('Nx',                     Nx)
			CALL ReadIntegerDataset('Ny',                     Ny)
			CALL ReadIntegerDataset('Nx_coarse',              Nx_coarse)
			CALL ReadIntegerDataset('Ny_coarse',              Ny_coarse)
			CALL ReadIntegerDataset('order_coarse_advection', order_coarse_advection)
			CALL ReadIntegerDataset('order_coarse_sound',     order_coarse_sound)
			CALL ReadIntegerDataset('order_fine_advection',   order_fine_advection)
			CALL ReadIntegerDataset('order_fine_sound',       order_fine_sound)
			CALL ReadDoubleDataset('dx', dx)
			CALL ReadDoubleDataset('dy', dy)
			CALL ReadDoubleDataset('dx_coarse', dx_coarse)
			CALL ReadDoubleDataset('dy_coarse', dy_coarse)
			CALL DeactivateGroup()
						
			! STEP TWO: Open group "problemdefinition" and read the parameters contained there		
			CALL ActivateGroup('problemdefinition')
			CALL ReadDoubleDataset('domain_height', domain_height)
			CALL ReadDoubleDataset('domain_width',  domain_width)
			CALL ReadDoubleDataset('sound_speed',   c_s)
			CALL ReadDoubleDataset('stabFreq',      stabFreq)
                        CALL ReadDoubleDataset('coriolisPar',   coriolisPar)
                        CALL ReadDoubleDataset('grav',          grav)
			CALL ReadIntegerDataset('BC',           BC)
			CALL DeactivateGroup()
						
			! STEP THREE: Open group "topology" and read values contained there
			CALL ActivateGroup('topology')
			CALL ReadIntegerDataset('nprocs_x', nprocs_x)
			CALL ReadIntegerDataset('nprocs_y', nprocs_y)
			CALL DeactivateGroup()
			
			!CALL ActivateGroup('parareal')
			!CALL ReadIntegerDataset('nthreads', Nthreads)
			!CALL ReadIntegerDataset('maxit', maxit)
			!CALL ReadIntegerDataset('Nparareal_restarts', Nparareal_restarts)
			!CALL DeactivateGroup()
			
			CALL ActivateGroup('integration')
			CALL ReadDoubleDataset('nu', nu)
			CALL ReadDoubleDataset('nu_coarse', nu_coarse)
			CALL DeactivateGroup()
			
			IF (echo_on) WRITE(*,'(A)') ' FVM parameters read succesfully '

			CALL CloseInput()

			dim = Nx*Ny*nr_fields
			dim_coarse = Nx_coarse*Ny_coarse*nr_fields

	END SUBROUTINE ReadFVMParameter
									 							 
END MODULE FVMParameters
