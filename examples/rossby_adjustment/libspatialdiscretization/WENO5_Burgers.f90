MODULE FiniteVolumes
! This module provides the implementation of a two-dimensional, finite difference WENO-5 scheme.
!
! Daniel Ruprecht, 8.2.2012
! ICS Lugano

USE omp_lib,       only : omp_get_thread_num

IMPLICIT NONE

INTEGER, PARAMETER :: nr_fields = 1, buffer_layout = 1

TYPE fdm_parameter
	INTEGER :: Nthreads, mpi_init_thread_flag
	LOGICAL :: echo_on
END TYPE

TYPE(fdm_parameter) :: param

PRIVATE
PUBLIC :: GetRHS, InitializeFiniteVolumes, FillinGhostcells, GhostLeft, GhostRight, GhostUp, GhostDown, GhostCorners, &
	GetMpiDatatypePar, PackSolution, UnpackSolution, nr_fields, buffer_layout, FillGhostCells

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: GhostCorners

! Define buffers storing ghost-cell values. NOTE: In the 1-D module, GhostUp and GhostDown are only listed
! to ensure compatibility, they are neither used nor allocated.
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: GhostLeft, GhostRight, GhostUp, GhostDown, GhostFluxLeft, GhostFluxRight, GhostFluxUp, GhostFluxDown

! Define buffers storing the horizontal cell and interface flux values
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FluxInt_hor, FluxCell_hor, FluxInt_ver, FluxCell_ver

! Define fixed parameters used by the WENO-5 method (see Shu .... e.g.)
DOUBLE PRECISION, PARAMETER, DIMENSION(3) :: weights_plus = (/ 0.3, 0.6, 0.1 /)
DOUBLE PRECISION, PARAMETER, DIMENSION(5) :: stencil_weights = (/ 2.0, 5.0, -1.0, -7.0, 11.0 /)*(1.0/6.0)
DOUBLE PRECISION, PARAMETER               :: coeff_1 = 13.0/12.0, coeff_2 = 1.0/4.0
DOUBLE PRECISION, PARAMETER               :: weno_tol = 1.0e-6
INTEGER,          PARAMETER               :: weno_n   = 2

CONTAINS

	SUBROUTINE GetRHS(Q, order_advection, order_diffusion, RQ, dy, dx, nu)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: dy, dx, nu
		INTEGER,                            INTENT(IN)  :: order_advection, order_diffusion

		INTEGER :: thread_nr, Nx, Ny, i, j
		DOUBLE PRECISION :: coeff_xx, coeff_yy

		thread_nr = omp_get_thread_num()
		
		! Nonlinear advection
		FluxCell_hor(  1,:,:,thread_nr) = 0.5*Q(1,:,:)*Q(1,:,:)
		GhostFluxLeft( 1,:,:,thread_nr) = 0.5*GhostLeft( 1,:,:,thread_nr)*GhostLeft( 1,:,:,thread_nr)
		GhostFluxRight(1,:,:,thread_nr) = 0.5*GhostRight(1,:,:,thread_nr)*GhostRight(1,:,:,thread_nr)
				
		FluxCell_ver( 1,:,:,thread_nr)  = 0.5*Q(1,:,:)*Q(1,:,:)
		GhostFluxUp(  1,:,:,thread_nr)  = 0.5*GhostUp(  1,:,:,thread_nr)*GhostUp(  1,:,:,thread_nr)
		GhostFluxDown(1,:,:,thread_nr)  = 0.5*GhostDown(1,:,:,thread_nr)*GhostDown(1,:,:,thread_nr)
			
		! Linear advection
		!FluxCell_hor(  1,:,:,thread_nr) = Q(1,:,:)
		!GhostFluxLeft( 1,:,:,thread_nr) = GhostLeft( 1,:,:,thread_nr)
		!GhostFluxRight(1,:,:,thread_nr) = GhostRight(1,:,:,thread_nr)
				
		!FluxCell_ver( 1,:,:,thread_nr)  = Q(1,:,:)
		!GhostFluxUp(  1,:,:,thread_nr)  = GhostUp(  1,:,:,thread_nr)
		!GhostFluxDown(1,:,:,thread_nr)  = GhostDown(1,:,:,thread_nr)
							
		! Now update interface values of horizontal flux
		CALL UpdateHorizontalFlux(Q, dble(1.0) )
		CALL UpdateVerticalFlux(Q,   dble(1.0) )
		
		!CALL UpdateHorizontalFlux_Upwind(Q)
		!CALL UpdateVerticalFlux_Upwind(Q)
		!FluxInt_ver(:,:,:,thread_nr) = 0.0
		
		! Compute flux divergence
		CALL GetFluxDivergence(RQ, dy, dx)
		
		! Add contribution from diffusion
		CALL AddDiffusion(Q, RQ, dx, dy, nu, order_diffusion)
		
		CONTAINS 
		
			SUBROUTINE UpdateHorizontalFlux_Upwind(Q)
				DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: Q
				
				INTEGER          :: Nx, Ny, i, j, thread_nr
				DOUBLE PRECISION :: Uvel
				DOUBLE PRECISION, PARAMETER :: weight = 0.5
								
				thread_nr = omp_get_thread_num()
				
				Ny = SIZE(Q,2)
				Nx = SIZE(Q,3)
							
				i=1
				DO j=1,Ny
					Uvel = 0.5
					FluxInt_hor(1,j,i,thread_nr) = Uvel*GhostLeft(1,j,1,thread_nr)
				END DO
				
				DO i=2,Nx
					DO j=1,Ny
						Uvel = 0.5
						FluxInt_hor(1,j,i,thread_nr) = Uvel*Q(1,j,i-1)
					END DO
				END DO
				
				i=Nx+1
				DO j=1,Ny
					Uvel = 0.5
					FluxInt_hor(1,j,i,thread_nr) = Uvel*Q(1,j,i-1)		
				END DO
										
			END SUBROUTINE UpdateHorizontalFlux_Upwind
			
			SUBROUTINE UpdateVerticalFlux_Upwind(Q)
				DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: Q
				
				INTEGER          :: Nx, Ny, i, j, thread_nr
				DOUBLE PRECISION :: Uvel, weight = 1.0/2.0
								
				thread_nr = omp_get_thread_num()
				
				Ny = SIZE(Q,2)
				Nx = SIZE(Q,3)			
			
			END SUBROUTINE UpdateVerticalFlux_Upwind
			
			SUBROUTINE AddDiffusion(Q, RQ, dx, dy, nu, order)
			
				DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)    :: Q
				DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: RQ
				DOUBLE PRECISION,                   INTENT(IN)    :: dx, dy, nu
				INTEGER,                            INTENT(IN)    :: order
			
				INTEGER :: Nx, Ny, i, j
				
				Ny = SIZE(Q,2)
				Nx = SIZE(Q,3)
							
				SELECT CASE (order)
				
					CASE (2)
			
					coeff_xx = nu/(dx*dx)
					coeff_yy = nu/(dy*dy)

				
					! -- Stencils using left ghost cells --
					i=1
					j=1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( Q(1,j,i+1) - 2.0*Q(1,j,i) + GhostLeft(1,j,1,thread_nr) ) + coeff_yy*( GhostUp(1,1,i,thread_nr) - 2.0*Q(1,j,i) + Q(1,j+1,i) )
					
					DO j=2,Ny-1
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( Q(1,j,i+1) - 2.0*Q(1,j,i) + GhostLeft(1,j,1,thread_nr) ) + coeff_yy*( Q(1,j-1,i) - 2.0*Q(1,j,i) + Q(1,j+1,i) )
					END DO
					
					j=Ny
					RQ(1,j,i) = RQ(1,j,i) +coeff_xx*( Q(1,j,i+1) - 2.0*Q(1,j,i) + GhostLeft(1,j,1,thread_nr) ) + coeff_yy*( Q(1,j-1,i) - 2.0*Q(1,j,i) + GhostDown(1,1,i,thread_nr) )
					
					! -- Inner stencils --
					DO i=2,Nx-1
						
						j=1
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( Q(1,j,i+1) - 2.0*Q(1,j,i) + Q(1,j,i-1) ) + coeff_yy*( GhostUp(1,1,i,thread_nr) - 2.0*Q(1,j,i) + Q(1,j+1,i) )
						
						DO j=2,Ny-1
							RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( Q(1,j,i+1) - 2.0*Q(1,j,i) + Q(1,j,i-1) ) + coeff_yy*( Q(1,j-1,i) - 2.0*Q(1,j,i) + Q(1,j+1,i) )
						END DO
						
						j=Ny
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( Q(1,j,i+1) - 2.0*Q(1,j,i) + Q(1,j,i-1) ) + coeff_yy*( Q(1,j-1,i) - 2.0*Q(1,j,i) + GhostDown(1,1,i,thread_nr) )
						
					END DO
					
					! -- Stencils using right ghost cell --
					
					i=Nx
					j=1
					RQ(1,j,i) = RQ(1,j,i) +coeff_xx*( GhostRight(1,j,1,thread_nr) - 2.0*Q(1,j,i) + Q(1,j,i-1) ) +coeff_yy*( GhostUp(1,1,i,thread_nr) - 2.0*Q(1,j,i) + Q(1,j+1,i) )
					
					DO j=2,Ny-1
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( GhostRight(1,j,1,thread_nr) - 2.0*Q(1,j,i) + Q(1,j,i-1) ) + coeff_yy*( Q(1,j-1,i) - 2.0*Q(1,j,i) + Q(1,j+1,i) )
					END DO
					
					j=Ny
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( GhostRight(1,j,1,thread_nr) - 2.0*Q(1,j,i) + Q(1,j,i-1) ) + coeff_yy*( Q(1,j-1,i) - 2.0*Q(1,j,i) + GhostDown(1,1,i,thread_nr) )
					
							
					CASE (4)
					
					coeff_xx = nu/(12.0*dx*dx)
					coeff_yy = nu/(12.0*dy*dy) 
					
					! -- Stencils using left ghost cell --
					
					i=1 ! ----- !
					j=1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -Q(1,j,i+2)               + 16.0*Q(1,j,i+1)               - 30.0*Q(1,j,i) + 16.0*GhostLeft(1,j,1,thread_nr) - GhostLeft(1,j,2,thread_nr) ) &
										  + coeff_yy*( -GhostUp(1,2,i,thread_nr) + 16.0*GhostUp(1,1,i,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i)                 - Q(1,j+2,i) )
				
					j=2
					RQ(1,j,i) = RQ(1,j,i)+ coeff_xx*( -Q(1,j,i+2)               + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Ghostleft(1,j,1,thread_nr) - GhostLeft(1,j,2,thread_nr) ) &
										  +coeff_yy*( -GhostUp(1,1,i,thread_nr) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i)                 - Q(1,j+2,i) )
					
					DO j=3,Ny-2
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*GhostLeft(1,j,1,thread_nr) - GhostLeft(1,j,2,thread_nr) ) &
											  + coeff_yy*( -Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i)                 - Q(1,j+2,i) )
					END DO
					
					j=Ny-1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*GhostLeft(1,j,1,thread_nr) - GhostLeft(1,j,2,thread_nr) ) &
										  + coeff_yy*(-Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i)                 - GhostDown(1,1,i,thread_nr) ) 
					
					j=Ny
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*GhostLeft(1,j,1,thread_nr) - GhostLeft(1,j,2,thread_nr) ) &
										  + coeff_yy*(-Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*GhostDown(1,1,i,thread_nr) - GhostDown(1,2,i,thread_nr) )
					 
					i=2 ! ----- !
					j=1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2)               + 16.0*Q(1,j,i+1)               - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - GhostLeft(1,j,1,thread_nr) )&
										  + coeff_yy*(-GhostUp(1,2,i,thread_nr) + 16.0*GhostUp(1,1,i,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) )
					
					j=2
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2)               + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - GhostLeft(1,j,1,thread_nr) ) &
										  + coeff_yy*(-GhostUp(1,1,i,thread_nr) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) )
					
					DO j=3,Ny-2
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - GhostLeft(1,j,1,thread_nr) ) &
											  + coeff_yy*(-Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) )
					END DO
					
					j=Ny-1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - GhostLeft(1,j,1,thread_nr) ) &
										  + coeff_yy*(-Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - GhostDown(1,1,i,thread_nr) )
					
					j=Ny
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1)                 - GhostLeft(1,j,1,thread_nr) ) &
										  + coeff_yy*(-Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*GhostDown(1,1,i,thread_nr) - GhostDown(1,2,i,thread_nr) )
					
					! -- Inner stencils --
					DO i=3,Nx-2
					
						j=1
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2)               + 16.0*Q(1,j,i+1)               - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
											  + coeff_yy*(-GhostUp(1,2,i,thread_nr) + 16.0*GhostUp(1,1,i,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) ) 
						
						j=2
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2)               + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
											  + coeff_yy*(-GhostUp(1,1,i,thread_nr) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) )	
						
						DO j=3,Ny-2
							RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
												  + coeff_yy*( -Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) )
						END DO
						
						j=Ny-1
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
											  + coeff_yy*(-Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - GhostDown(1,1,i,thread_nr) )
						
						j=Ny
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -Q(1,j,i+2) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1)                 - Q(1,j,i-2) ) &
											  + coeff_yy*( -Q(1,j-2,i) + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*GhostDown(1,1,i,thread_nr) - GhostDown(1,2,i,thread_nr) )
									
					END DO
					
					! -- Stencils using right ghost cell --
					
					i=Nx-1 ! ---- !
					j=1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,1,thread_nr) + 16.0*Q(1,j,i+1)               - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
										  + coeff_yy*( -GhostUp(1,2,i,thread_nr)    + 16.0*GhostUp(1,1,i,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) )
					
					j=2
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,1,thread_nr) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
										  + coeff_yy*( -GhostUp(1,1,i,thread_nr)    + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) ) 
					
					DO j=3,Ny-2
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,1,thread_nr) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
											  + coeff_yy*( -Q(1,j-2,i)                  + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) ) 
					END DO
					
					j=Ny-1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,1,thread_nr) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
										  + coeff_yy*( -Q(1,j-2,i)                  + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - GhostDown(1,1,i,thread_nr) )
					 
					j=Ny
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,1,thread_nr) + 16.0*Q(1,j,i+1) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1)                 - Q(1,j,i-2) ) &
										  + coeff_yy*( -Q(1,j-2,i)                  + 16.0*Q(1,j-1,i) - 30.0*Q(1,j,i) + 16.0*GhostDown(1,1,i,thread_nr) - GhostDown(1,2,i,thread_nr) )
					
					i=Nx ! ---- !
					j=1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,2,thread_nr) + 16.0*GhostRight(1,j,1,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
										  + coeff_yy*( -GhostUp(1,2,i,thread_nr)    + 16.0*GhostUp(1,1,i,thread_nr)    - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) ) 
								
					j=2
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*(-GhostRight(1,j,2,thread_nr) + 16.0*GhostRight(1,j,1,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
										  + coeff_yy*(-GhostUp(1,1,i,thread_nr)    + 16.0*Q(1,j-1,i)                  - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) )
					
					DO j=3,Ny-2
						RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,2,thread_nr) + 16.0*GhostRight(1,j,1,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
											  + coeff_yy*( -Q(1,j-2,i)                  + 16.0*Q(1,j-1,i)                  - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - Q(1,j+2,i) ) 
					END DO
					
					j=Ny-1
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,2,thread_nr) + 16.0*GhostRight(1,j,1,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1) - Q(1,j,i-2) ) &
										  + coeff_yy*( -Q(1,j-2,i)                  + 16.0*Q(1,j-1,i)                  - 30.0*Q(1,j,i) + 16.0*Q(1,j+1,i) - GhostDown(1,1,i,thread_nr) )
					
					j=Ny
					RQ(1,j,i) = RQ(1,j,i) + coeff_xx*( -GhostRight(1,j,2,thread_nr) + 16.0*GhostRight(1,j,1,thread_nr) - 30.0*Q(1,j,i) + 16.0*Q(1,j,i-1)                 - Q(1,j,i-2) ) &
										  + coeff_yy*( -Q(1,j-2,i)			        + 16.0*Q(1,j-1,i)                  - 30.0*Q(1,j,i) + 16.0*GhostDown(1,1,i,thread_nr) - GhostDown(1,2,i,thread_nr) )
							
					CASE DEFAULT
						WRITE(*,*) 'Encounted unexpected value for order_advection! For Burgers equation, only order two and four are implemented. Now exiting.'
						STOP
						
				END SELECT 
				
				END SUBROUTINE AddDiffusion
				
	END SUBROUTINE GetRHS
	
	
	SUBROUTINE InitializeFiniteVolumes(Nx_max, Ny_max, Nthreads, mpi_init_thread, echo_on)
	
		INTEGER, INTENT(IN) :: Nthreads, Ny_max, Nx_max, mpi_init_thread
		LOGICAL, INTENT(IN) :: echo_on
		
		INTEGER :: i, thread_nr
		
		param%echo_on         = echo_on
		param%nthreads        = Nthreads
		param%mpi_init_thread_flag = mpi_init_thread
				
		! The index in the "thread-dimension" starts with zero, so that
		! it directly coincides with thread numbers
		ALLOCATE(GhostLeft(     1:nr_fields, Ny_max, 3, 0:Nthreads-1))
		ALLOCATE(GhostRight(    1:nr_fields, Ny_max, 3, 0:Nthreads-1))
		ALLOCATE(GhostFluxLeft( 1:nr_fields, Ny_max, 3, 0:Nthreads-1))
		ALLOCATE(GhostFluxRight(1:nr_fields, Ny_max, 3, 0:Nthreads-1))
		ALLOCATE(GhostCorners(1:nr_fields, 4, 0:Nthreads-1))

		ALLOCATE(GhostUp(       1:nr_fields, 3, Nx_max, 0:Nthreads-1))
		ALLOCATE(GhostDown(     1:nr_fields, 3, Nx_max, 0:Nthreads-1))
		ALLOCATE(GhostFluxUp(   1:nr_fields, 3, Nx_max, 0:Nthreads-1))
		ALLOCATE(GhostFluxDown( 1:nr_fields, 3, Nx_max, 0:Nthreads-1))
		
		! If there are Nx cells, there are Nx+1 interfaces
		ALLOCATE(FluxInt_hor( 1:nr_fields, Ny_max,   Nx_max+1, 0:Nthreads-1))
		ALLOCATE(FluxCell_hor(1:nr_fields, Ny_max,   Nx_max,   0:Nthreads-1))
		ALLOCATE(FluxInt_ver( 1:nr_fields, Ny_max+1, Nx_max,   0:Nthreads-1))
		ALLOCATE(FluxCell_ver(1:nr_fields, Ny_max,   Nx_max,   0:Nthreads-1))
		
		! Now perform first-touch initialization, i.e. every thread initializes its
		! part of the buffers
		
		!$OMP PARALLEL private(thread_nr)
		!$OMP DO schedule(static)
		DO i=1,Nthreads
			thread_nr = omp_get_thread_num()
			GhostLeft(     :,:,:,thread_nr) = 0.0
			GhostRight(    :,:,:,thread_nr) = 0.0
			GhostFluxLeft( :,:,:,thread_nr) = 0.0
			GhostFluxRight(:,:,:,thread_nr) = 0.0
			GhostUp(       :,:,:,thread_nr) = 0.0
			GhostDown(     :,:,:,thread_nr) = 0.0
			GhostFluxUp(   :,:,:,thread_nr) = 0.0
			GhostFluxDown( :,:,:,thread_nr) = 0.0
			FluxInt_hor(   :,:,:,thread_nr) = 0.0
			FluxCell_hor(  :,:,:,thread_nr) = 0.0
			FluxInt_ver(   :,:,:,thread_nr) = 0.0
			FluxCell_ver(  :,:,:,thread_nr) = 0.0
		END DO
		!$OMP END DO
		!$OMP END PARALLEL
		
	END SUBROUTINE InitializeFiniteVolumes
	
	SUBROUTINE UpdateHorizontalFlux(Qcell, max_vel)
		
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: Qcell
		DOUBLE PRECISION,                   INTENT(IN) :: max_vel
		
		DOUBLE PRECISION, DIMENSION(nr_fields, 6) :: Qcell_local, Fcell_local
		INTEGER :: Ny, Nx, i, j, thread_nr
		
		thread_nr = omp_get_thread_num()
		
		! Out of the global fields Qcell and FluxQcell, updated interface
		! values of the flux are computed
		Ny = SIZE(Qcell,2)
		Nx = SIZE(Qcell,3)
		
		DO j=1,Ny
			
			! Left boundary
			i=0
			
			! Initialize "moving window" arrays
			Qcell_local(:,1) = GhostLeft(:,j,3,thread_nr) 
			Qcell_local(:,2) = GhostLeft(:,j,2,thread_nr) 
			Qcell_local(:,3) = GhostLeft(:,j,1,thread_nr)
			Qcell_local(:,4) = Qcell(:,j,i+1)
			Qcell_local(:,5) = Qcell(:,j,i+2)
			Qcell_local(:,6) = Qcell(:,j,i+3)
			
			Fcell_local(:,1) = GhostFluxLeft(:,j,3,thread_nr)
			Fcell_local(:,2) = GhostFluxLeft(:,j,2,thread_nr)
			Fcell_local(:,3) = GhostFluxLeft(:,j,1,thread_nr)
			Fcell_local(:,4) = FluxCell_hor( :,j,i+1,thread_nr)
			Fcell_local(:,5) = FluxCell_hor( :,j,i+2,thread_nr)
			Fcell_local(:,6) = FluxCell_hor( :,j,i+3,thread_nr)
			
			! Reconstruct interface value from local stencil
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_hor(:,j,i+1,thread_nr))

			! Inner cells
			DO i=1,Nx-3
					! Update local "moving window" arrays
					CALL ShiftLocalValues(Qcell_local, Fcell_local, Qcell(:,j,i+3), FluxCell_hor(:,j,i+3,thread_nr))
					CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_hor(:,j,i+1,thread_nr))
			END DO
			
			! Right boundary
			i=Nx-2
			CALL ShiftLocalValues(Qcell_local, Fcell_local, GhostRight(:,j,1,thread_nr), GhostFluxRight(:,j,1,thread_nr))
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_hor(:,j,i+1,thread_nr))
					
			i=Nx-1
			CALL ShiftLocalValues(Qcell_local, Fcell_local, GhostRight(:,j,2,thread_nr), GhostFluxRight(:,j,2,thread_nr))
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_hor(:,j,i+1,thread_nr))
					
			i=Nx
			CALL ShiftLocalValues(Qcell_local, Fcell_local, GhostRight(:,j,3,thread_nr), GhostFluxRight(:,j,3,thread_nr))
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_hor(:,j,i+1,thread_nr))
		
		END DO
		
		CONTAINS 
		
			PURE SUBROUTINE ShiftLocalValues(Qcell_local, Fcell_local, Qcell_new, Fcell_new)
			
				DOUBLE PRECISION, DIMENSION(nr_fields, 6), INTENT(INOUT) :: Qcell_local, Fcell_local
				DOUBLE PRECISION, DIMENSION(nr_fields),    INTENT(IN)    :: Qcell_new, Fcell_new
				
				INTEGER :: i
								
				! Shift all values one to the left				
				DO i=1,5
					Qcell_local(:,i) = Qcell_local(:,i+1)
					Fcell_local(:,i) = Fcell_local(:,i+1)
				END DO
				
				! ..and fill in the new values to the right
				Qcell_local(:,6) = Qcell_new
				Fcell_local(:,6) = Fcell_new				
							
			END SUBROUTINE ShiftLocalValues
				
	END SUBROUTINE UpdateHorizontalFlux
	
	
	SUBROUTINE UpdateVerticalFlux(Qcell, max_vel)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: Qcell
		DOUBLE PRECISION,                   INTENT(IN) :: max_vel
		
		DOUBLE PRECISION, DIMENSION(nr_fields, 6) :: Qcell_local, Fcell_local
		INTEGER :: Ny, Nx, i, j, thread_nr	
		
		thread_nr = omp_get_thread_num()
		
		! Out of the global fields Qcell and FluxQcell, updated interface
		! values of the flux are computed
		Ny = SIZE(Qcell,2)
		Nx = SIZE(Qcell,3)
		
		DO i=1,Nx
		
			j=0

			Qcell_local(:,6) = GhostUp(:,3,i,thread_nr)
			Qcell_local(:,5) = GhostUp(:,2,i,thread_nr)
			Qcell_local(:,4) = GhostUp(:,1,i,thread_nr)
			Qcell_local(:,3) = Qcell(:,1,i)
			Qcell_local(:,2) = Qcell(:,2,i)
			Qcell_local(:,1) = Qcell(:,3,i)
			
			Fcell_local(:,6) = GhostFluxUp(:,3,i,thread_nr)
			Fcell_local(:,5) = GhostFluxUp(:,2,i,thread_nr)
			Fcell_local(:,4) = GhostFluxUp(:,1,i,thread_nr)
			Fcell_local(:,3) = FluxCell_ver(:,1,i,thread_nr)
			Fcell_local(:,2) = FluxCell_ver(:,2,i,thread_nr)
			Fcell_local(:,1) = FluxCell_ver(:,3,i,thread_nr)
			
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_ver(:,j+1,i,thread_nr))

			DO j=1,Ny-3
				CALL ShiftLocalValues(Qcell_local, Fcell_local, Qcell(:,j+3,i), FluxCell_ver(:,j+3,i,thread_nr))
				CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_ver(:,j+1,i,thread_nr))
			END DO
			
			j=Ny-2
			CALL ShiftLocalValues(Qcell_local, Fcell_local, GhostDown(:,1,i,thread_nr), GhostFluxDown(:,1,i,thread_nr))
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_ver(:,j+1,i,thread_nr))
							
			j=Ny-1
			CALL ShiftLocalValues(Qcell_local, Fcell_local, GhostDown(:,2,i,thread_nr), GhostFluxDown(:,2,i,thread_nr))
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_ver(:,j+1,i,thread_nr))
				
			j=Ny
			CALL ShiftLocalValues(Qcell_local, Fcell_local, GhostDown(:,3,i,thread_nr), GhostFluxDown(:,3,i,thread_nr))
			CALL ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel, FluxInt_ver(:,j+1,i,thread_nr))
				
		END DO

		CONTAINS 
		
			PURE SUBROUTINE ShiftLocalValues(Qcell_local, Fcell_local, Qcell_new, Fcell_new)
			
				DOUBLE PRECISION, DIMENSION(nr_fields, 6), INTENT(INOUT) :: Qcell_local, Fcell_local
				DOUBLE PRECISION, DIMENSION(nr_fields),    INTENT(IN)    :: Qcell_new, Fcell_new
				
				INTEGER :: j
								
				! Shift all values one down				
				DO j=6, 2, -1 ! count backwards
					Qcell_local(:,j) = Qcell_local(:,j-1)
					Fcell_local(:,j) = Fcell_local(:,j-1)
				END DO
				
				! ..and fill in the new values at the bottom
				Qcell_local(:,1) = Qcell_new
				Fcell_local(:,1) = Fcell_new
							
			END SUBROUTINE ShiftLocalValues	
	
	END SUBROUTINE UpdateVerticalFlux
		
	PURE SUBROUTINE ReconstructInterfaceValue(Qcell_local, Fcell_local, local_vel, Fint)
	
		DOUBLE PRECISION, DIMENSION(nr_fields, 6), INTENT(IN)  :: Qcell_local, Fcell_local
		DOUBLE PRECISION,                          INTENT(IN)  :: local_vel
		DOUBLE PRECISION, DIMENSION(nr_fields),    INTENT(OUT) :: Fint
		
		DOUBLE PRECISION, DIMENSION(nr_fields, 5) :: Fcell_plus, Fcell_minus
		DOUBLE PRECISION, DIMENSION(nr_fields)    :: Fint_plus,  Fint_minus
		
		! From the local values of Q and flux of Q, retrieve the positive and negative
		! components of the flux function Fcell_plus and Fcell_minus
		CALL GetLocalFluxSplit(Qcell_local, Fcell_local, local_vel, Fcell_plus, Fcell_minus)
	
		! Perform a WENO reconstruction of the interface value for the positive and negative
		! component
		CALL GetWenoInterfaceValue(Fcell_plus, Fcell_minus, Fint_plus, Fint_minus)
		
		! The final interface value is the sum of the values reconstructed from the positive
		! and negative flux component
		Fint = Fint_plus + Fint_minus
		
		CONTAINS
		
			PURE SUBROUTINE GetLocalFluxSplit(Qcell_local, Fcell_local, local_vel, Fcell_plus, Fcell_minus)
	
				DOUBLE PRECISION, DIMENSION(nr_fields, 6), INTENT(IN)  :: Qcell_local, Fcell_local
				DOUBLE PRECISION,                          INTENT(IN)  :: local_vel
				DOUBLE PRECISION, DIMENSION(nr_fields, 5), INTENT(OUT) :: Fcell_plus, Fcell_minus
			
				! Lax-Friedrichs flux-splitting
				Fcell_plus( :,1:5) = 0.5*( Fcell_local(:,1:5) + local_vel*Qcell_local(:,1:5) )
				Fcell_minus(:,1:5) = 0.5*( Fcell_local(:,2:6) - local_vel*Qcell_local(:,2:6) )
										
			END SUBROUTINE GetLocalFluxSplit
	
			PURE SUBROUTINE GetWenoInterfaceValue(Fcell_plus, Fcell_minus, Fint_plus, Fint_minus)
		
				DOUBLE PRECISION, DIMENSION(nr_fields,5), INTENT(IN)  :: Fcell_plus, Fcell_minus
				DOUBLE PRECISION, DIMENSION(nr_fields),   INTENT(OUT) :: Fint_plus, Fint_minus
				
				DOUBLE PRECISION, DIMENSION(nr_fields,3) :: pol_values, beta, alpha
				DOUBLE PRECISION                         :: alpha_sum_inv
				INTEGER i,j
				
				pol_values(:,1) = stencil_weights(1)*Fcell_plus(:,3) + stencil_weights(2)*Fcell_plus(:,4) + stencil_weights(3)*Fcell_plus(:,5)
				pol_values(:,2) = stencil_weights(3)*Fcell_plus(:,2) + stencil_weights(2)*Fcell_plus(:,3) + stencil_weights(1)*Fcell_plus(:,4)
				pol_values(:,3) = stencil_weights(1)*Fcell_plus(:,1) + stencil_weights(4)*Fcell_plus(:,2) + stencil_weights(5)*Fcell_plus(:,3) 
				
				! Second, compute smoothness measures
				beta(:,1) = coeff_1*(     Fcell_plus(:,3) - 2.0*Fcell_plus(:,4) +     Fcell_plus(:,5) )**2 &
						  + coeff_2*( 3.0*Fcell_plus(:,3) - 4.0*Fcell_plus(:,4) +     Fcell_plus(:,5) )**2
				beta(:,2) = coeff_1*(     Fcell_plus(:,2) - 2.0*Fcell_plus(:,3) +     Fcell_plus(:,4) )**2 &
						  + coeff_2*(     Fcell_plus(:,2)                       -     Fcell_plus(:,4) )**2
				beta(:,3) = coeff_1*(     Fcell_plus(:,1) - 2.0*Fcell_plus(:,2) +     Fcell_plus(:,3) )**2 &
						  + coeff_2*(     Fcell_plus(:,1) - 4.0*Fcell_plus(:,2) + 3.0*Fcell_plus(:,3) )**2
								
				! Third, compute weights out of the smoothness measures
				alpha(:,1) = weights_plus(1)/( beta(:,1) + weno_tol )**weno_n	
				alpha(:,2) = weights_plus(2)/( beta(:,2) + weno_tol )**weno_n
				alpha(:,3) = weights_plus(3)/( beta(:,3) + weno_tol )**weno_n
									
				! Fourth, normalize the weights.					
				DO j=1,nr_fields
					alpha_sum_inv = 1.0/SUM(alpha(j,:))
					alpha(j,:) = alpha_sum_inv*alpha(j,:)
				END DO

				! Finally, compute the superposition of the three candidate-stencil values using the computed weights
				Fint_plus(:) = alpha(:,1)*pol_values(:,1) + alpha(:,2)*pol_values(:,2) + alpha(:,3)*pol_values(:,3)
				
				! *** Now perform corresponding computation for Fcell_minus ****
				pol_values(:,1) = stencil_weights(5)*Fcell_minus(:,3) + stencil_weights(4)*Fcell_minus(:,4) &
					+ stencil_weights(1)*Fcell_minus(:,5) 
				pol_values(:,2) = stencil_weights(1)*Fcell_minus(:,2) + stencil_weights(2)*Fcell_minus(:,3) &
					+ stencil_weights(3)*Fcell_minus(:,4)
				pol_values(:,3) = stencil_weights(3)*Fcell_minus(:,1) + stencil_weights(2)*Fcell_minus(:,2) &
					+ stencil_weights(1)*Fcell_minus(:,3)
				
				! Smoothness measures
				beta(:,1) = coeff_1*(       Fcell_minus(:,3) - 2.0*Fcell_minus(:,4) +     Fcell_minus(:,5) )**2 &
						  + coeff_2*(   3.0*Fcell_minus(:,3) - 4.0*Fcell_minus(:,4) +     Fcell_minus(:,5) )**2
				beta(:,2) = coeff_1*(       Fcell_minus(:,2) - 2.0*Fcell_minus(:,3) +     Fcell_minus(:,4) )**2 &
						  + coeff_2*(       Fcell_minus(:,2)                        -     Fcell_minus(:,4) )**2
				beta(:,3) = coeff_1*(       Fcell_minus(:,1) - 2.0*Fcell_minus(:,2) +     Fcell_minus(:,3) )**2 &
						  + coeff_2*(       Fcell_minus(:,1) - 4.0*Fcell_minus(:,2) + 3.0*Fcell_minus(:,3) )**2
					
				! Compute weights	
				alpha(:,1) = weights_plus(3)/( beta(:,1) + weno_tol )**weno_n
				alpha(:,2) = weights_plus(2)/( beta(:,2) + weno_tol )**weno_n
				alpha(:,3) = weights_plus(1)/( beta(:,3) + weno_tol )**weno_n

				! Normalize weights
				DO j=1,nr_fields							
					alpha_sum_inv = 1.0/SUM(alpha(j,:))
					alpha(j,:) = alpha_sum_inv*alpha(j,:)									
				END DO
				
				Fint_minus(:) = alpha(:,1)*pol_values(:,1) + alpha(:,2)*pol_values(:,2) + alpha(:,3)*pol_values(:,3)
			
			END SUBROUTINE GetWenoInterfaceValue

	END SUBROUTINE ReconstructInterfaceValue

	SUBROUTINE GetFluxDivergence(RQ, dy, dx)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: dy, dx
		
		DOUBLE PRECISION :: coeff_x, coeff_y
		INTEGER          :: i,j, Ny, Nx, thread_nr
		
		Ny        = SIZE(RQ,2)
		Nx        = SIZE(RQ,3)
		coeff_x   = 1.0/dx
		coeff_y   = 1.0/dy
		thread_nr = omp_get_thread_num()
		
		DO i=1,Nx
			DO j=1,Ny
				RQ(:,j,i) = -coeff_x*( FluxInt_hor(:,j,i+1,thread_nr) - FluxInt_hor(:,j,i,thread_nr) ) &
							-coeff_y*( FluxInt_ver(:,j,i,thread_nr) - FluxInt_ver(:,j+1,i,thread_nr) )
			END DO
		END DO
		
	END SUBROUTINE GetFluxDivergence
	
	SUBROUTINE FillinGhostcells(Qleft, Qright, Qup, Qdown, Qupleft, Qupright, Qdownleft, Qdownright, Nx, Ny, Nghost)
	
		DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Qleft, Qright, Qup, Qdown, Qupleft, Qupright, Qdownleft, Qdownright
		INTEGER, INTENT(IN) :: Nx, Ny, Nghost

		INTEGER :: thread_nr, counter, i, j, k
		
		thread_nr = omp_get_thread_num()
		
		IF (param%mpi_init_thread_flag==1) THEN
		
			GhostCorners(:,1,0:param%nthreads-1) = RESHAPE( Qupleft(   1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			GhostCorners(:,2,0:param%nthreads-1) = RESHAPE( Qupright(  1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			GhostCorners(:,3,0:param%nthreads-1) = RESHAPE( Qdownright(1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			GhostCorners(:,4,0:param%nthreads-1) = RESHAPE( Qdownleft( 1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			
			! In MPI_THREAD_FUNNELED mode, this routine is only called by the master thread
			counter=1
			DO thread_nr=0,param%nthreads-1
				DO i=1,Nghost
					DO j=1,Ny
						DO k=1,nr_fields
							GhostLeft(k,j,Nghost-i+1,thread_nr)  = Qleft( counter)
							GhostRight(k,j,i,thread_nr)          = Qright(counter)
							counter=counter+1
						END DO
					END DO
				END DO
			END DO
			
			counter=1
			DO thread_nr=0,param%nthreads-1
				DO i=1,Nx
					DO j=1,Nghost
						DO k=1,nr_fields
							GhostUp(k,Nghost-j+1,i,thread_nr)   = Qup(  counter)
							GhostDown(k,j,i,thread_nr)          = Qdown(counter)
							counter=counter+1
						END DO
					END DO
				END DO	
			END DO	

		ELSE IF ( (param%mpi_init_thread_flag==2) .or. (param%mpi_init_thread_flag==3) ) THEN

			counter=1
			DO i=1,Nghost
				DO j=1,Ny
					DO k=1,nr_fields
						GhostLeft(k,j,Nghost-i+1,thread_nr)  = Qleft( counter)
						GhostRight(k,j,i,thread_nr)          = Qright(counter)
						counter=counter+1
					END DO
				END DO
			END DO	
			
			counter=1
			DO i=1,Nx
				DO j=1,Nghost
					DO k=1,nr_fields
						GhostUp(k,Nghost-j+1,i,thread_nr)   = Qup(  counter)
						GhostDown(k,j,i,thread_nr)          = Qdown(counter)
						counter=counter+1
					END DO
				END DO
			END DO								
		END IF
				
		IF (param%echo_on) WRITE(*,'(A)') ' Communication completed, ghost-cell buffers updated... '
				
	END SUBROUTINE FillinGhostcells	
	
	SUBROUTINE GetMpiDatatypePar(Nx, Ny, Nghost, blocklengths, indices, length, length_singleThread)
		
		INTEGER,                 INTENT(IN)  :: Nx, Ny, Nghost
		INTEGER, DIMENSION(:,:), INTENT(OUT) :: blocklengths, indices
		INTEGER, DIMENSION(:),   INTENT(OUT) :: length, length_singleThread
		
		INTEGER :: i
			
		! --- Upper and lower ghost-cells ---
		length(2)              = Nx*param%Nthreads
		length_singleThread(2) = Nx
		length(7)              = Nx*param%Nthreads
		length_singleThread(7) = Nx

		! Upper row datatype: For every column  one block of Nghost*nr_fields values has to be defined								
		blocklengths(1:length(2),2) = Nghost*nr_fields
		DO i=1,length(2)
			indices(i,2) = (i-1)*Ny*nr_fields
		END DO		

		! Lower row datatype
		blocklengths(1:length(7),7) = Nghost*nr_fields	
		DO i=1,length(7)
			indices(i,7) = i*Ny*nr_fields-Nghost*nr_fields
		END DO
		
		! --- Left and right ghost-cells ---
		length(4)              = param%Nthreads
		length_singleThread(4) = 1
		length(5)              = param%Nthreads
		length_singleThread(5) = 1
		
		! Left column datatype
		blocklengths(1:length(4),4) = Nghost*Ny*nr_fields
		DO i=1,length(4)
			indices(i,4) = (i-1)*Nx*Ny*nr_fields
		END DO

		! Right column datatype
		blocklengths(1:length(5),5) = Nghost*Ny*nr_fields
		DO i=1,length(5)
			indices(i,5) = (Nx-Nghost)*Ny*nr_fields + (i-1)*Nx*Ny*nr_fields
		END DO
			
		! --- Corner ghost-cells ---
		length(1)              = nr_fields*param%Nthreads
		length_singleThread(1) = nr_fields
		length(3)              = nr_fields*param%Nthreads
		length_singleThread(3) = nr_fields
		length(6)              = nr_fields*param%Nthreads
		length_singleThread(6) = nr_fields
		length(8)              = nr_fields*param%Nthreads
		length_singleThread(8) = nr_fields
																		
		! Upper left
		blocklengths(1:length(1),1) = 1
		DO i=1,length(1)
			indices(i,1) = (i-1)*Nx*Ny
		END DO				

		! Upper right
		blocklengths(1:length(3),3) = 1
		DO i=1,length(3)
			indices(i,3) = (Nx-1)*Ny + (i-1)*Nx*Ny
		END DO
				
		! Lower left
		blocklengths(1:length(6),6) = 1
		DO i=1,length(6)
			indices(i,6) = (Ny-1) + (i-1)*Nx*Ny
		END DO
		
		! Lower right
		blocklengths(1:length(8),8) = 1
		DO i=1,length(8)
			indices(i,8) = i*Nx*Ny-1
		END DO

	END SUBROUTINE GetMpiDatatypePar
	
	PURE SUBROUTINE UnpackSolution(Q, Y, nr_fields, Ny, Nx)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Q
		DOUBLE PRECISION, DIMENSION(:),     INTENT(IN)  :: Y
		INTEGER,                            INTENT(IN)  :: nr_fields, Ny, Nx
		
		INTEGER :: i,j, k, counter
		
		counter=1
		DO i=1,Nx
			DO j=1,Ny
				DO k=1,nr_fields
					Q(k,j,i) = Y(counter)
					counter = counter+1
				END DO
			END DO
		END DO
				
	END SUBROUTINE UnpackSolution


	PURE SUBROUTINE PackSolution(Y, Q, nr_fields, Ny, Nx)
		
		DOUBLE PRECISION, DIMENSION(:),     INTENT(OUT) :: Y
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		INTEGER,                            INTENT(IN)  :: nr_fields, Ny, Nx

		INTEGER :: i, j, k, counter
		
		counter=1
		
		DO i=1,Nx
			DO j=1,Ny
				DO k=1,nr_fields
					Y(counter) = Q(k,j,i)
					counter=counter+1
				END DO
			END DO
		END DO

		
	END SUBROUTINE PackSolution	
	
	!-------------------------------------------------------------------------------------------------|
	! This subroutine fills the arrays storing the ghost cell values depending on the 
	! prescribed boundary condition. Values are stored into the module arrays GhostLeft, GhostRight,
	! GhostUp, GhostDown for every indicated field in the solution array Q.
	!
	! That is, GhostUp(:,:,k,thread_nr) contains the upper ghost cell values for the solution component
	! Q(:,:,k,thread_nr) etc.
	!
	! In a future hybrid MPI / Open MP version, this routine will handle communcation with neighboring
	! subdomains.
	SUBROUTINE FillGhostCells(Q, Nghost, BC, thread_nr)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		INTEGER,						    INTENT(IN)  :: Nghost, thread_nr, BC
		
		INTEGER :: i,j,k,Nx, Ny
		
		!thread_nr = omp_get_thread_num()
		
		Ny = SIZE(Q, 2)
		Nx = SIZE(Q, 3)
		
		SELECT CASE (BC)
		
			! Periodic in all four directions
			CASE (1)
			
				! Fill horizontal ghost cells
				DO i=1,Nghost
					DO j=1,Ny
						GhostLeft( :,j,i,thread_nr) = Q(:,j, Nx-i+1)
						GhostRight(:,j,i,thread_nr) = Q(:,j,i)
					END DO
				END DO
				
				! Fill vertical ghost cells
				DO i=1,Nx
					DO j=1,Nghost
						GhostUp(  :,j,i,thread_nr) = Q(:,Ny-j+1,i)
						GhostDown(:,j,i,thread_nr) = Q(:,j,i)
					END DO
				END DO
				
				! Upper left corner
				GhostCorners(:,1,thread_nr) = Q(:,Ny,Nx)
				
				! Upper right corner
				GhostCorners(:,2,thread_nr) = Q(:,Ny,1)
				
				! Lower right corner
				GhostCorners(:,3,thread_nr) = Q(:,1,1)
				
				! Lower left corner
				GhostCorners(:,4,thread_nr) = Q(:,1,Nx)
				
			! Outflow boundary conditions in all four direction
			CASE (2)
			
				! Fill horizontal ghost cells
				DO i=1,Nghost
					DO j=1,Ny
						GhostLeft( :,j,i,thread_nr)  = Q(:,j,1)
						GhostRight(:,j,i,thread_nr)  = Q(:,j, Nx)
					END DO
				END DO
			
				! Fill vertical ghost cells	
				DO i=1,Nx
					DO j=1,Nghost
						GhostUp(  :,j,i,thread_nr) = Q(:,1,i)
						GhostDOwn(:,j,i,thread_nr) = Q(:,Ny,i)
					END DO
				END DO
						
			! Periodic in the horizontal, solid walls in the vertical
			CASE (3)
			
				! Iterate over the different fields of the solution
				DO k=1,nr_fields
				
					! Fill horizontal ghost cells
					DO i=1,Nghost
						DO j=1,Ny
							GhostLeft( j,i,k,thread_nr) = Q(j, Nx-i+1, k)
							GhostRight(j,i,k,thread_nr) = Q(j, i, k)
						END DO
					END DO
					
					! Fill vertical ghost cells: To realize solid wall BC, invert the sign for the velocities
					! but keep the sign for pressure and buoyancy; cf. eq. (7.17), p. 137 in LeVeque
					!
					! Fields one and two are the velocities.
					IF ( (k==1) .OR. (k==2) ) THEN
						
						DO i=1,Nx
							DO j=1,Nghost
								GhostUp(  j, i, k, thread_nr) = -Q(j, i, k)
								GhostDown(j, i, k, thread_nr) = -Q(Ny-j+1, i, k)
							END DO
						END DO
				
					ELSE
					
						DO i=1,Nx
							DO j=1,Nghost
								GhostUp(  j, i, k, thread_nr) = Q(j, i, k)
								GhostDown(j, i, k, thread_nr) = Q(Ny-j+1,i,k)
							END DO
						END DO
						
					END IF
					
				END DO
			
			CASE DEFAULT
			
				WRITE(*, *) ' No implementation available for selected boundary condition. Now exiting.'
				STOP
			
		END SELECT
		
	END SUBROUTINE FillGhostCells			
	
END MODULE FiniteVolumes