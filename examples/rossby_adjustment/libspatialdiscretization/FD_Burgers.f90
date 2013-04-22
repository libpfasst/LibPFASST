MODULE FiniteVolumes
!> This module provides a set of subroutines that can be used to evaluate the right hand side of a PDE in
!> flux form, that is
!>
!>					q_t + F(q)_x + G(q)_y = 0
!>
!> by a finite difference discretization in conservation form, leading to a semi-discrete equation
!>
!>					q_t = - (1/dx)*[ F_(i+1/2, j) - F_(i-1/2,j) ] - (1/dy)*[ G_(i, j+1/2) - G(i, j-1/2) ]		(1)
!> 
!> See LeVeque, "Finite Volume Methods for Hyperbolic Problems", Section 19.4, p.443ff for details. The weights for the stencils for
!> the interface fluxes are from Durran, "Numerical Methods for Fluid Dynamics", Table 5.2, p. 248.
!>
!> The geometry has to be a simple square, discretized by a set of rectangles with constant height dy and width dx. 
!>
!> Each solution component of q has to contain Ny times Nx cell-centered values. To define a problem, a routine returning the right hand
!> side of (1) for a given Q has to be provided. This can then be used in a generic ODE solver for problems of the form
!>
!>					q_t = L(q)
!>
!> Note that the spatial discretization will result in restrictions on the useable time steps. 
!>
!> Also note that the functions related to the spatial discretization work on arrays of dimension Ny times Ny, where indices (1,1) correspond 
!> to a value in the upper left corner, (1,Nx) to a value in the upper right corner, (Ny,1) to a value in the lower left and (Ny, Nx) to a value in 
!> the lower right corner.
!>
!> NOTE: The wrapper module 'RHSFunctions' reshapes a column vector q with a length equal to the total DoF into a Ny x Nx x nr_fields buffer that can
!>       be used in the right hand side functions provided below.
!>
!> Daniel Ruprecht
!> Institute of Computational Science, Lugano
!> February, 17th, 2011.

USE omp_lib,              only : omp_get_thread_num

IMPLICIT NONE

! Only the problem specific routines implementing the right hand side functions and the initialization routines are accessible, all
! other routines are private and can not be accessed from outside the module.
!
! Note that the right hand side functions should usually only be called by the wrapper module 'RHSFunctions' and not directly.
!

INTEGER, PARAMETER :: buffer_layout = 0
INTEGER, PARAMETER :: nr_fields     = 1

TYPE fdm_parameter
	INTEGER :: Nthreads, mpi_init_thread_flag
	LOGICAL :: echo_on
END TYPE

TYPE(fdm_parameter) :: param

PRIVATE
PUBLIC GetRHS, InitializeFiniteVolumes, FillinGhostcells, GhostLeft, GhostRight, GhostUp, GhostDown, GhostCorners, &
	GetMpiDatatypePar, PackSolution, UnpackSolution, buffer_layout, nr_fields, FillGhostCells

! Specifiy temporary buffers used inside this module. The last index of each buffer indicates the thread using this part of it.
! Every thread initializes its own part, so that hopefully the buffer is allocated in memory in a suitable way, allowing for
! quick access by the different threads (first touch policy).
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: GhostLeft, GhostRight, GhostUp, GhostDown, FluxHor, FluxVer
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: Divergence, Uadv, Vadv, GhostCorners

CONTAINS

	! ------ Pure advection -----|
	SUBROUTINE GetRHS(Q, order_advection, order_diffusion, RQ, dx, dy, nu)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		INTEGER,                            INTENT(IN)  :: order_advection, order_diffusion
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: dx, dy, nu
		
		DOUBLE PRECISION :: coeff_xx, coeff_yy
		INTEGER          :: thread_nr, u_field_nr, v_field_nr, Nx, Ny, i, j
					

		CALL UpdateAdvectionVelocityNonlinear( Q(:,:,1), Q(:,:,1), u_field_nr = 1, v_field_nr = 1)
	
		thread_nr = omp_get_thread_num()
		Uadv(:,:,thread_nr) = 0.0
		Vadv(:,:,thread_nr) = 0.0
		
		Ny = SIZE(Q,1)
		Nx = SIZE(Q,2)

		RQ(:,:,1)                = DBLE(0.0)
		FluxHor(:,:,1,thread_nr) = DBLE(0.0)
		FluxVer(:,:,1,thread_nr) = DBLE(0.0)
			
		! Add flux from horizontal advection
		CALL UpdateHorizontalAdvectionFlux( 1, Q(:,:,1), 1, order_advection)
				
		! Add flux from vertical advection
		CALL UpdateVerticalAdvectionFlux(   1, Q(:,:,1), 1, order_advection)
				
		CALL UpdateFluxDivergence(1, 1, RQ, dx, dy)
				
		! REMARK:: For Burger's equation, second or fourth order centered finite difference stencils are hardcoded here. The stencils are not
		! in conservation form and there is only one solution field, so this involves simply runnin two nested loops over the two spatial
		! directions.
		SELECT CASE (order_diffusion)
		
			CASE (2)
	
			coeff_xx = nu/(dx*dx)
			coeff_yy = nu/(dy*dy)

		
			! -- Stencils using left ghost cells --
			i=1
			j=1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( Q(j,i+1,1) - 2.0*Q(j,i,1) + GhostLeft(j,1,1,thread_nr) ) + coeff_yy*( GhostUp(1,i,1,thread_nr) - 2.0*Q(j,i,1) + Q(j+1,i,1) )
			
			DO j=2,Ny-1
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( Q(j,i+1,1) - 2.0*Q(j,i,1) + GhostLeft(j,1,1,thread_nr) ) + coeff_yy*( Q(j-1,i,1) - 2.0*Q(j,i,1) + Q(j+1,i,1) )
			END DO
			
			j=Ny
			RQ(j,i,1) = RQ(j,i,1) +coeff_xx*( Q(j,i+1,1) - 2.0*Q(j,i,1) + GhostLeft(j,1,1,thread_nr) ) + coeff_yy*( Q(j-1,i,1) - 2.0*Q(j,i,1) + GhostDown(1,i,1,thread_nr) )
			
			! -- Inner stencils --
			DO i=2,Nx-1
				
				j=1
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( Q(j,i+1,1) - 2.0*Q(j,i,1) + Q(j,i-1,1) ) + coeff_yy*( GhostUp(1,i,1,thread_nr) - 2.0*Q(j,i,1) + Q(j+1,i,1) )
				
				DO j=2,Ny-1
					RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( Q(j,i+1,1) - 2.0*Q(j,i,1) + Q(j,i-1,1) ) + coeff_yy*( Q(j-1,i,1) - 2.0*Q(j,i,1) + Q(j+1,i,1) )
				END DO
				
				j=Ny
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( Q(j,i+1,1) - 2.0*Q(j,i,1) + Q(j,i-1,1) ) + coeff_yy*( Q(j-1,i,1) - 2.0*Q(j,i,1) + GhostDown(1,i,1,thread_nr) )
				
			END DO
			
			! -- Stencils using right ghost cell --
			
			i=Nx
			j=1
			RQ(j,i,1) = RQ(j,i,1) +coeff_xx*( GhostRight(j,1,1,thread_nr) - 2.0*Q(j,i,1) + Q(j,i-1,1) ) +coeff_yy*( GhostUp(1,i,1,thread_nr) - 2.0*Q(j,i,1) + Q(j+1,i,1) )
			
			DO j=2,Ny-1
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( GhostRight(j,1,1,thread_nr) - 2.0*Q(j,i,1) + Q(j,i-1,1) ) + coeff_yy*( Q(j-1,i,1) - 2.0*Q(j,i,1) + Q(j+1,i,1) )
			END DO
			
			j=Ny
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( GhostRight(j,1,1,thread_nr) - 2.0*Q(j,i,1) + Q(j,i-1,1) ) + coeff_yy*( Q(j-1,i,1) - 2.0*Q(j,i,1) + GhostDown(1,i,1,thread_nr) )
			
					
			CASE (4)
			
			coeff_xx = nu/(12.0*dx*dx)
			coeff_yy = nu/(12.0*dy*dy) 
			
			! -- Stencils using left ghost cell --
			
			i=1 ! ----- !
			j=1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -Q(j,i+2,1)               + 16.0*Q(j,i+1,1)               - 30.0*Q(j,i,1) + 16.0*GhostLeft(j,1,1,thread_nr) - GhostLeft(j,2,1,thread_nr) ) &
							      + coeff_yy*( -GhostUp(2,i,1,thread_nr) + 16.0*GhostUp(1,i,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1)                 - Q(j+2,i,1) )
		
			j=2
			RQ(j,i,1) = RQ(j,i,1)+ coeff_xx*( -Q(j,i+2,1)               + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Ghostleft(j,1,1,thread_nr) - GhostLeft(j,2,1,thread_nr) ) &
   					              +coeff_yy*( -GhostUp(1,i,1,thread_nr) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1)                 - Q(j+2,i,1) )
			
			DO j=3,Ny-2
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*GhostLeft(j,1,1,thread_nr) - GhostLeft(j,2,1,thread_nr) ) &
							          + coeff_yy*( -Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1)                 - Q(j+2,i,1) )
			END DO
			
			j=Ny-1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*GhostLeft(j,1,1,thread_nr) - GhostLeft(j,2,1,thread_nr) ) &
						          + coeff_yy*(-Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1)                 - GhostDown(1,i,1,thread_nr) ) 
			
			j=Ny
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*GhostLeft(j,1,1,thread_nr) - GhostLeft(j,2,1,thread_nr) ) &
						          + coeff_yy*(-Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*GhostDown(1,i,1,thread_nr) - GhostDown(2,i,1,thread_nr) )
			 
			i=2 ! ----- !
			j=1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1)               + 16.0*Q(j,i+1,1)               - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - GhostLeft(j,1,1,thread_nr) )&
								  + coeff_yy*(-GhostUp(2,i,1,thread_nr) + 16.0*GhostUp(1,i,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) )
			
			j=2
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1)               + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - GhostLeft(j,1,1,thread_nr) ) &
						          + coeff_yy*(-GhostUp(1,i,1,thread_nr) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) )
			
			DO j=3,Ny-2
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - GhostLeft(j,1,1,thread_nr) ) &
							          + coeff_yy*(-Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) )
			END DO
			
			j=Ny-1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - GhostLeft(j,1,1,thread_nr) ) &
						          + coeff_yy*(-Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - GhostDown(1,i,1,thread_nr) )
			
			j=Ny
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1)                 - GhostLeft(j,1,1,thread_nr) ) &
						          + coeff_yy*(-Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*GhostDown(1,i,1,thread_nr) - GhostDown(2,i,1,thread_nr) )
			
			! -- Inner stencils --
			DO i=3,Nx-2
			
				j=1
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1)               + 16.0*Q(j,i+1,1)               - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
									  + coeff_yy*(-GhostUp(2,i,1,thread_nr) + 16.0*GhostUp(1,i,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) ) 
				
				j=2
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1)               + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
							          + coeff_yy*(-GhostUp(1,i,1,thread_nr) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) )	
				
				DO j=3,Ny-2
					RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
								          + coeff_yy*( -Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) )
				END DO
				
				j=Ny-1
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
									  + coeff_yy*(-Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - GhostDown(1,i,1,thread_nr) )
				
				j=Ny
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -Q(j,i+2,1) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1)                 - Q(j,i-2,1) ) &
							          + coeff_yy*( -Q(j-2,i,1) + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*GhostDown(1,i,1,thread_nr) - GhostDown(2,i,1,thread_nr) )
							
			END DO
			
			! -- Stencils using right ghost cell --
			
			i=Nx-1 ! ---- !
			j=1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,1,1,thread_nr) + 16.0*Q(j,i+1,1)               - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
								  + coeff_yy*( -GhostUp(2,i,1,thread_nr)    + 16.0*GhostUp(1,i,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) )
			
			j=2
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,1,1,thread_nr) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
						          + coeff_yy*( -GhostUp(1,i,1,thread_nr)    + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) ) 
			
			DO j=3,Ny-2
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,1,1,thread_nr) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
									  + coeff_yy*( -Q(j-2,i,1)                  + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) ) 
			END DO
			
			j=Ny-1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,1,1,thread_nr) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
								  + coeff_yy*( -Q(j-2,i,1)                  + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - GhostDown(1,i,1,thread_nr) )
			 
			j=Ny
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,1,1,thread_nr) + 16.0*Q(j,i+1,1) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1)                 - Q(j,i-2,1) ) &
								  + coeff_yy*( -Q(j-2,i,1)                  + 16.0*Q(j-1,i,1) - 30.0*Q(j,i,1) + 16.0*GhostDown(1,i,1,thread_nr) - GhostDown(2,i,1,thread_nr) )
			
			i=Nx ! ---- !
			j=1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,2,1,thread_nr) + 16.0*GhostRight(j,1,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
						          + coeff_yy*( -GhostUp(2,i,1,thread_nr)    + 16.0*GhostUp(1,i,1,thread_nr)    - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) ) 
						
			j=2
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*(-GhostRight(j,2,1,thread_nr) + 16.0*GhostRight(j,1,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
								  + coeff_yy*(-GhostUp(1,i,1,thread_nr)    + 16.0*Q(j-1,i,1)                  - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) )
			
			DO j=3,Ny-2
				RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,2,1,thread_nr) + 16.0*GhostRight(j,1,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
									  + coeff_yy*( -Q(j-2,i,1)                  + 16.0*Q(j-1,i,1)                  - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - Q(j+2,i,1) ) 
			END DO
			
			j=Ny-1
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,2,1,thread_nr) + 16.0*GhostRight(j,1,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1) - Q(j,i-2,1) ) &
								  + coeff_yy*( -Q(j-2,i,1)                  + 16.0*Q(j-1,i,1)                  - 30.0*Q(j,i,1) + 16.0*Q(j+1,i,1) - GhostDown(1,i,1,thread_nr) )
			
			j=Ny
			RQ(j,i,1) = RQ(j,i,1) + coeff_xx*( -GhostRight(j,2,1,thread_nr) + 16.0*GhostRight(j,1,1,thread_nr) - 30.0*Q(j,i,1) + 16.0*Q(j,i-1,1)                 - Q(j,i-2,1) ) &
								  + coeff_yy*(-Q(j-2,i,1)				   + 16.0*Q(j-1,i,1)                  - 30.0*Q(j,i,1) + 16.0*GhostDown(1,i,1,thread_nr) - GhostDown(2,i,1,thread_nr) )
					
			CASE DEFAULT
				!WRITE(*,*) 'Encounted unexpected value for order_advection! For Burgers equation, only order two and four are implemented. Now exiting.'
				!STOP
				
		END SELECT 
		
		IF (param%echo_on) WRITE(*,*) ' Finished evaluation GetBurgersRHS ... '
		
	END SUBROUTINE GetRHS
		
	! ===== Initialization Routines ===== |
	
	SUBROUTINE InitializeFiniteVolumes(Nx_max, Ny_max, Nthreads, mpi_init_thread, echo_on)
		
		INTEGER, INTENT(IN) :: Nx_max,     &!
							   Ny_max,	   &!
							   Nthreads,   &!
							   mpi_init_thread
							   
		LOGICAL, INTENT(IN) :: echo_on
							   	
		INTEGER :: i, thread_nr
		
		param%echo_on         = echo_on
		param%nthreads        = Nthreads
		param%mpi_init_thread_flag = mpi_init_thread
		
		! Ghost cell values of divergence are stored in Ghost(:,:,nr_fields+1,:)
		ALLOCATE(GhostLeft( Ny_max, 3,  1:nr_fields+1, 0:Nthreads-1))
		ALLOCATE(GhostRight(Ny_max, 3,  1:nr_fields+1, 0:Nthreads-1))
		ALLOCATE(GhostUp(    3, Nx_max, 1:nr_fields+1, 0:Nthreads-1))
		ALLOCATE(GhostDown(  3, Nx_max, 1:nr_fields+1, 0:Nthreads-1))
		ALLOCATE(GhostCorners(4, 1:nr_fields, 0:Nthreads-1))
		
		ALLOCATE(Divergence(   Ny_max, Nx_max, 0:Nthreads-1))
		
		ALLOCATE(Uadv(Ny_max, Nx_max+1, 0:Nthreads-1))
		ALLOCATE(Vadv(Ny_max+1, Nx_max, 0:Nthreads-1))
		
		ALLOCATE(FluxHor(Ny_max, Nx_max+1, 1:nr_fields, 0:Nthreads-1))
		ALLOCATE(FluxVer(Ny_max+1, Nx_max, 1:nr_fields, 0:Nthreads-1))
						
		! Now let every thread initialize its part of each buffer so that 
		! by first touch policy the allocated memory is "close" the the CPU
		! handling the thread.
		
		!$OMP PARALLEL private(thread_nr)
		!$OMP DO
		DO i=0,Nthreads-1
		
			thread_nr = omp_get_thread_num()
			
			GhostLeft( :,:,:,thread_nr)  = 0.0
			GhostRight(:,:,:,thread_nr)  = 0.0
			GhostUp(   :,:,:,thread_nr)  = 0.0
			GhostDown( :,:,:,thread_nr)  = 0.0
			GhostCorners(:,:,thread_nr)  = 0.0
			
			Divergence(:,:, thread_nr) = 0.0
			
			Uadv(:,:,thread_nr) = 0.0
			Vadv(:,:,thread_nr) = 0.0
			
			FluxHor(:,:,:,thread_nr) = 0.0
			FluxVer(:,:,:,thread_nr) = 0.0
						
		END DO
		!$OMP END DO
		!$OMP END PARALLEL
	
		IF (param%echo_on) WRITE(*, '(A, I2)') ' Module FiniteVolumes successfully initialized. Solution components : ', nr_fields
	
	END SUBROUTINE InitializeFiniteVolumes
	
	!------------------------------------------------------------------|
	! This function evaluates the flux divergence
	!
	! - [ (1/dx)*( F_(i+1/2, j) - F_(i-1/2, j) ) + (1/dy)*( G_(i, j+1/2) - G_(i, j-1/2) ) ]
	!
	! for components with indices field_nr_start : field_nr_end. This allows to update, e.g.
	! only the velocity components in a partially split scheme.
	!
	SUBROUTINE UpdateFluxDivergence(field_nr_start, field_nr_end, RQ, dx, dy)
	
		INTEGER, INTENT(IN)                               :: field_nr_start, &!
							                                 field_nr_end     !
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: RQ               ! Buffer for RHS
		DOUBLE PRECISION,                   INTENT(IN)    :: dx, dy
	
		DOUBLE PRECISION :: coeff_x, coeff_y
		INTEGER :: i, j, field_nr, thread_nr, Nx, Ny
		
		Ny = SIZE(RQ, 1)
		Nx = SIZE(RQ, 2)
		
		thread_nr = omp_get_thread_num()
		
		coeff_x = DBLE(1.0)/dx
		coeff_y = DBLE(1.0)/dy
		
		DO field_nr=field_nr_start,field_nr_end
			DO i=1,Nx
				DO j=1,Ny
					RQ(j,i,field_nr) = RQ(j,i,field_nr) - coeff_x*( FluxHor(j,i+1,field_nr,thread_nr) - FluxHor(j,i,field_nr,thread_nr)) &
						- coeff_y*( FluxVer(j,i,field_nr,thread_nr) - FluxVer(j+1,i,field_nr,thread_nr) )					
				END DO
			END DO
		END DO
		
	END SUBROUTINE UpdateFluxDivergence
	
	! ---------------------------------------------------------------- |
	! For nonlinear advection, the interface values of the advective
	! velocity have to be recomputed for every evaluation of the
	! advective fluxes.
	!
	! Interface values are generated by this routine using linear
	! interpolating the provided cell-centered values. The results are 
	! stored in the module-wide buffers Uadv, Vadv.
	!
	! Note that as dx and dz, that is the height and width of the cells,
	! are assumed to be constant over the mesh, the linear interpolation
	! simplifies to just averaging.
	SUBROUTINE UpdateAdvectionVelocityNonlinear(U_cell, V_cell, u_field_nr, v_field_nr)
	
		DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: U_cell, V_cell
		INTEGER,					      INTENT(IN) :: u_field_nr, &!
		                                                v_field_nr   !		
		INTEGER i, j, thread_nr, Nx, Ny
		
		DOUBLE PRECISION :: coeff_x, coeff_y
		
		Ny = SIZE(U_cell, 1)
		Nx = SIZE(U_cell, 2)
		thread_nr = omp_get_thread_num()
		
		! First, compute horizontal advective velocities
		
		coeff_x = 0.5
		coeff_y = 0.5
		
		! Values using left ghost cells
		i=1
		DO j=1,Ny
			!coeff_x = 0.25*( SIGN( dble(1.0), GhostLeft(j,1,u_field_nr, thread_nr)) + SIGN( dble(1.0),U_cell(j,i)) + 2.0)
			Uadv(j,i,thread_nr) = coeff_x*( GhostLeft(j,1,u_field_nr, thread_nr) + U_cell(j,i) )
		END DO
		
		! Values using no ghost cells
		DO i=2,Nx
			DO j=1,Ny
				!coeff_x = 0.25*( SIGN( dble(1.0), U_cell(j,i-1)) + SIGN( dble(1.0), U_cell(j,i)) + 2.0 )
				Uadv(j,i,thread_nr) = coeff_x*( U_cell(j,i-1) + U_cell(j,i) )
			END DO
		END DO
		
		! Values using right ghost cells
		i=Nx+1
		DO j=1,Ny
			!coeff_x = 0.25*( SIGN( dble(1.0), U_cell(j,i-1)) + SIGN( dble(1.0), GhostRight(j,1,u_field_nr, thread_nr)) + 2.0 )
			Uadv(j,i,thread_nr) = coeff_x*( U_cell(j,i-1) + GhostRight(j,1,u_field_nr,thread_nr) )
		END DO
		
		! Second, compute vertical advective velocities
		DO i=1,Nx
			
			j=1
			!coeff_y = 0.25*( SIGN( dble(1.0), GhostUp(1,i,v_field_nr,thread_nr)) + SIGN( dble(1.0), V_cell(j,i)) + 2.0 )
			Vadv(j,i,thread_nr) = coeff_y*( GhostUp(1,i,v_field_nr,thread_nr) + V_cell(j,i) )
		
			DO j=2,Ny
				!coeff_y = 0.25*( SIGN( dble(1.0), V_cell(j-1,i)) + SIGN( dble(1.0),V_cell(j,i)) + 2.0 )
				Vadv(j,i,thread_nr) = coeff_y*( V_cell(j-1,i) + V_cell(j,i) )
			END DO
			
			j=Ny+1
			!coeff_y = 0.25*( SIGN( dble(1.0), V_cell(j-1,i)) + SIGN(dble(1), GhostDown(1,i,v_field_nr, thread_nr)) + 2.0 )
			Vadv(j,i,thread_nr) = coeff_y*( V_cell(j-1,i) + GhostDown(1,i,v_field_nr,thread_nr) )
			
		END DO
			
	END SUBROUTINE UpdateAdvectionVelocityNonlinear
		
	! -------------------------------------------------------------------|
	! Computes Ny x Nx+1 interface values for the horizontal fluxes
	! F_(i+1/2, j), i=0...Nx given Ny x Nx cell values of the unknown Q.
	!
	! The order of the approximation can be set to 1,2,3,4,5,6
	!
	! The implemented formulats can be found for example in Wicker, Skamarock MWR(2002)
	! or in "Numerical Methods for Fluid Dynamics", Durran, Springer(2010)	
	SUBROUTINE UpdateHorizontalAdvectionFlux(F_field_nr, Q, field_nr, order)
	
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Q           ! Array of Ny x Nx cell values of advected quantity
			INTEGER,						  INTENT(IN) :: F_field_nr,&! Index of FluxHor(:,:,i,:) to write into, that is, solution component to which this flux contributes
															field_nr,  &! The number of the field that Q represents. Required to access the correct ghost cell values.
															order       ! Order of accuracy of flux computation.
			DOUBLE PRECISION, DIMENSION(3) :: weights
			DOUBLE PRECISION :: coeff
			INTEGER :: i, j, thread_nr, Nx, Ny
			
			Ny = SIZE(Q, 1)
			Nx = SIZE(Q, 2)
			thread_nr = omp_get_thread_num()
			
			SELECT CASE (order)
			
				CASE (1)
				
					weights(1) = 1.0/2.0
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(1)*Q(j,i) ) - 0.5*ABS(Uadv(j,i,thread_nr))*( Q(j,i) - GhostLeft(j,1,field_nr,thread_nr) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=2,Nx
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-1) + weights(1)*Q(j,i) ) &
								- 0.5*ABS(Uadv(j,i,thread_nr))*( Q(j,i) - Q(j,i-1) )
						END DO
					END DO
					
					! Fluxes accessing a ghost cell to the right
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-1) &
							+ weights(1)*GhostRight(j,1,field_nr,thread_nr) ) - 0.5*ABS(Uadv(j,i,thread_nr))*( GhostRight(j,1,field_nr,thread_nr) &
							- Q(j,i-1) )
					END DO
					
				CASE (2)
				
					weights(1) = 1.0/2.0
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(1)*Q(j,i) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=2,Nx
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-1) + weights(1)*Q(j,i) )
						END DO
					END DO
					
					! Fluxes accessing a ghost cell to the right
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-1) &
							+ weights(1)*GhostRight(j,1,field_nr,thread_nr) ) 
					END DO
									
				CASE (3)
				
					weights(1) = -1.0
					weights(2) =  7.0
					coeff = DBLE(1.0)/DBLE(12.0)
					
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) ) &
							- coeff*ABS(Uadv(j,i,thread_nr))*( -Q(j,i+1) + DBLE(3.0)*Q(j,i) & 
								- DBLE(3.0)*GhostLeft(j,1,field_nr,thread_nr) + GhostLeft(j,2,field_nr,thread_nr) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(2)*Q(j,i-1) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) ) &
							- coeff*ABS(Uadv(j,i,thread_nr))*( -Q(j,i+1) + DBLE(3.0)*Q(j,i) - DBLE(3.0)*Q(j,i-1) + GhostLeft(j,1,field_nr,thread_nr) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=3,Nx-1
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) &
								+ weights(2)*Q(j,i) + weights(1)*Q(j,i+1) ) &
								- coeff*ABS(Uadv(j,i,thread_nr))*( -Q(j,i+1) + DBLE(3.0)*Q(j,i) - DBLE(3.0)*Q(j,i-1) + Q(j,i-2) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) &
							+ weights(2)*Q(j,i) + weights(1)*GhostRight(j,1,field_nr,thread_nr)) &
							- coeff*ABS(Uadv(j,i,thread_nr))*( -GhostRight(j,1,field_nr,thread_nr) + DBLE(3.0)*Q(j,i) & 
							- DBLE(3.0)*Q(j,i-1) + Q(j,i-2) )
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) + &
							weights(2)*GhostRight(j,1,field_nr,thread_nr) + weights(1)*GhostRight(j,2,field_nr,thread_nr) ) &
							- coeff*ABS(Uadv(j,i,thread_nr))*( -GhostRight(j,2,field_nr,thread_nr) &
								+ DBLE(3.0)*GhostRight(j,1,field_nr,thread_nr) - DBLE(3.0)*Q(j,i-1) + Q(j,i-2) )
					END DO
					
				CASE (4)

					weights(1) = -1.0
					weights(2) =  7.0
					coeff = DBLE(1.0)/DBLE(12.0)
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(2)*Q(j,i-1) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=3,Nx-1
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-2) &
								+ weights(2)*Q(j,i-1) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) &
							+ weights(2)*Q(j,i) + weights(1)*GhostRight(j,1,field_nr,thread_nr))
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) + &
							weights(2)*GhostRight(j,1,field_nr,thread_nr) + weights(1)*GhostRight(j,2,field_nr,thread_nr) )
					END DO
					
				CASE (5)
				
					weights(1) = 1.0
					weights(2) = -8.0
					weights(3) = 37.0
					coeff = DBLE(1.0)/DBLE(60.0)
					
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,3,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,2,field_nr,thread_nr) + weights(3)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) ) &
							-coeff*ABS(Uadv(j,i,thread_nr))*( -GhostLeft(j,3,field_nr,thread_nr) + DBLE(5.0)*GhostLeft(j,2,field_nr,thread_nr) &
								- DBLE(10.0)*GhostLeft(j,1,field_nr,thread_nr) + DBLE(10.0)*Q(j,i) - DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) + weights(3)*Q(j,i-1) + weights(3)*Q(j,i) &
							+ weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) ) &
							-coeff*ABS(Uadv(j,i,thread_nr))*( -GhostLeft(j,2,field_nr,thread_nr) + DBLE(5.0)*GhostLeft(j,1,field_nr,thread_nr) &
								-DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) - DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
					END DO
					i=3
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(2)*Q(j,i-2) + weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) ) &
							-coeff*ABS(Uadv(j,i,thread_nr))*( -GhostLeft(j,1,field_nr,thread_nr) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) &
								+ DBLE(10.0)*Q(j,i) - DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=4,Nx-2
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )    &
								-coeff*ABS(Uadv(j,i,thread_nr))*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) &
								- DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx-1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*GhostRight(j,1,field_nr,thread_nr) ) &
								-coeff*ABS(Uadv(j,i,thread_nr))*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) &
								 - DBLE(5.0)*Q(j,i+1) + GhostRight(j,1,field_nr,thread_nr) )
					END DO
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
							+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*GhostRight(j,1,field_nr,thread_nr) &
							+ weights(1)*GhostRight(j,2,field_nr,thread_nr) ) - coeff*ABS(Uadv(j,i,thread_nr))*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) &
							- DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) - DBLE(5.0)*GhostRight(j,1,field_nr,thread_nr) &
							+ GhostRight(j,2,field_nr,thread_nr) )
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
							+ weights(3)*Q(j,i-1) + weights(3)*GhostRight(j,1,field_nr,thread_nr) + weights(2)*GhostRight(j,2,field_nr,thread_nr) &
							+ weights(1)*GhostRight(j,3,field_nr,thread_nr) ) &
							-coeff*ABS(Uadv(j,i,thread_nr))*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) &
								+ DBLE(10.0)*GhostRight(j,1,field_nr,thread_nr) - DBLE(5.0)*GhostRight(j,2,field_nr,thread_nr) &
								+ GhostRight(j,3,field_nr,thread_nr) )
					END DO					
												
				CASE (6)
				
					weights(1) = 1.0
					weights(2) = -8.0
					weights(3) = 37.0
					coeff = DBLE(1.0)/DBLE(60.0)
					
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,3,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(3)*GhostLeft(j,1,field_nr,thread_nr) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
					END DO
					i=3
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(2)*Q(j,i-2) + weights(3)*Q(j,i-1) &
							+ weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=4,Nx-2
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx-1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) &
								+ weights(2)*Q(j,i+1) + weights(1)*GhostRight(j,1,field_nr,thread_nr) )
					END DO
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) &
								+ weights(2)*GhostRight(j,1,field_nr,thread_nr) + weights(1)*GhostRight(j,2,field_nr,thread_nr) )
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*Uadv(j,i,thread_nr)*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*GhostRight(j,1,field_nr,thread_nr) &
								+ weights(2)*GhostRight(j,2,field_nr,thread_nr) + weights(1)*GhostRight(j,3,field_nr,thread_nr) )
					END DO
			
				CASE DEFAULT
					WRITE(*,*) 'No implementation available for flux of selected order. Now exiting.'
					STOP
			END SELECT
			
	END SUBROUTINE UpdateHorizontalAdvectionFlux
	
	! -------------------------------------------------------------------|
	! Computes Ny+1 x Nx interface values for the vertical fluxes
	! G_(i, j+1/2), j=0...Ny given Ny x Nx cell values of the unknown Q.
	!
	! IMPORTANT: Note that the index j=1 corresponds to the upper row of values
	!		     while j=Ny is the lowest row. That is, Q(1,1) is the value in the
	!            upper, left cell for example.
	!
	! The order of the approximation can be set to 1,2,3,4,6
	!
	! The implemented formulats can be found for example in Wicker, Skamarock MWR(2002)
	! or in "Numerical Methods for Fluid Dynamics", Durran, Springer(2010)
	SUBROUTINE UpdateVerticalAdvectionFlux(F_field_nr, Q, field_nr, order)
	
		DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Q			 ! Ny x Nx cell values of one solution component
		INTEGER,						  INTENT(IN) :: F_field_nr,& !
														field_nr, &  ! Index of solution component, required to access correct 
																	 ! ghost cells.
														order        ! Order of approximation, can be 1, 2,3,4 or 6
														
		DOUBLE PRECISION, DIMENSION(3) :: weights
		DOUBLE PRECISION :: coeff
		INTEGER :: i, j, thread_nr, Nx, Ny
		
		Ny = SIZE(Q, 1)
		Nx = SIZE(Q, 2)
		thread_nr = omp_get_thread_num()
		
		SELECT CASE (order)
		
			CASE (1)
			
				weights(1) = 1.0/2.0
				DO i=1,Nx
					! Fluxes accessing an upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + Vadv(j,i,thread_nr)*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(1)*Q(j,i) ) &
						-0.5*ABS(Vadv(j,i,thread_nr))*( GhostUp(1,i,field_nr,thread_nr) - Q(j,i) )
					
					! Fluxes not accessing any ghost cell
					DO j=2,Ny
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + Vadv(j,i,thread_nr)*( weights(1)*Q(j-1,i) + weights(1)*Q(j,i) ) &
							-0.5*ABS(Vadv(j,i,thread_nr))*( Q(j-1,i) - Q(j,i) )
					END DO
					
					! Fluxes accessing a lower ghost cell
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + Vadv(j,i,thread_nr)*( weights(1)*Q(j-1,i) &
						+ weights(1)*GhostDown(1,i,field_nr,thread_nr) ) &
						-0.5*ABS(Vadv(j,i,thread_nr))*( Q(j-1,i) - GhostDown(1,i,field_nr,thread_nr) )
				END DO
			
			CASE (2)
			
				weights(1) = 1.0/2.0
				DO i=1,Nx
					! Fluxes accessing an upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + Vadv(j,i,thread_nr)*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(1)*Q(j,i) ) 					
					! Fluxes not accessing any ghost cell
					DO j=2,Ny
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + Vadv(j,i,thread_nr)*( weights(1)*Q(j-1,i) + weights(1)*Q(j,i) ) 
					END DO
					
					! Fluxes accessing a lower ghost cell
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + Vadv(j,i,thread_nr)*( weights(1)*Q(j-1,i) &
						+ weights(1)*GhostDown(1,i,field_nr,thread_nr) ) 					
				END DO				
			
			CASE (3)
			
				weights(1) = -1.0
				weights(2) =  7.0
				coeff = DBLE(1.0)/DBLE(12.0)
				
				DO i=1,Nx
					! Fluxes accessing at least one upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(2,i,field_nr, thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr, thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) &
						- coeff*ABS(Vadv(j,i,thread_nr))*( -GhostUp(2,i,field_nr, thread_nr) + DBLE(3.0)*GhostUp(1,i,field_nr,thread_nr) &
							- DBLE(3.0)*Q(j,i) + Q(j+1,i) )
					j=2  
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-1,i) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i)) &
						- coeff*ABS(Vadv(j,i,thread_nr))*( -GhostUp(1,i,field_nr,thread_nr) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*Q(j,i) + Q(j+1,i) )
					
					! Fluxes accessing no ghost cell
					DO j=3,Ny-1
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
							+ weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) &
							- coeff*ABS(Vadv(j,i,thread_nr))*( -Q(j-2,i) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*Q(j,i) + Q(j+1,i) )
					END DO
					
					! Fluxes accessing at least one lower ghost cell
					j=Ny
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*Q(j,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr)) &
						- coeff*ABS(Vadv(j,i,thread_nr))*( -Q(j-2,i) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*Q(j,i) + GhostDown(1,i,field_nr,thread_nr))
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*GhostDown(1,i,field_nr,thread_nr) + weights(1)*GhostDown(2,i,field_nr,thread_nr) ) &
						- coeff*ABS(Vadv(j,i,thread_nr))*( -Q(j-2,i) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*GhostDown(1,i,field_nr,thread_nr) &
						+ GhostDown(2,i,field_nr,thread_nr))
					
				END DO
			
			CASE (4)
			
				weights(1) = -1.0
				weights(2) =  7.0
				coeff = DBLE(1.0)/DBLE(12.0)
				
				DO i=1,Nx
					! Fluxes accessing at least one upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(2,i,field_nr, thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr, thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) 
					j=2  
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-1,i) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i)) 					
					! Fluxes accessing no ghost cell
					DO j=3,Ny-1
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
							+ weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) 
					END DO
					
					! Fluxes accessing at least one lower ghost cell
					j=Ny
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*Q(j,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr)) 
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*GhostDown(1,i,field_nr,thread_nr) + weights(1)*GhostDown(2,i,field_nr,thread_nr) )					
				END DO			
				
			CASE (5)
				
				weights(1) =  1.0
				weights(2) = -8.0
				weights(3) = 37.0
				coeff = DBLE(1.0)/DBLE(60.0)
								
				DO i=1,Nx
					! Fluxes accessing upper ghost cells
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(3,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(2,i,field_nr,thread_nr) + weights(3)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) ) &
						-coeff*ABS(Vadv(j,i,thread_nr))*( GhostUp(3,i,field_nr,thread_nr) - DBLE(5.0)*GhostUp(2,i,field_nr,thread_nr) &
							+ DBLE(10.0)*GhostUp(1,i,field_nr,thread_nr) - DBLE(10.0)*Q(j,i) + DBLE(5.0)*Q(j+1,i) - Q(j+2,i)  )
						
					j=2
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(2,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr,thread_nr) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) &
						+ weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) ) &
						-coeff*ABS(Vadv(j,i,thread_nr))*( GhostUp(2,i,field_nr,thread_nr) - DBLE(5.0)*GhostUp(1,i,field_nr,thread_nr) &
								+ DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) + DBLE(5.0)*Q(j+1,i) - Q(j+2,i) )
								
					j=3
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-2,i) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) ) &
						-coeff*ABS(Vadv(j,i,thread_nr))*( GhostUp(1,i,field_nr,thread_nr) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) &
							- DBLE(10.0)*Q(j,i) + DBLE(5.0)*Q(j+1,i) - Q(j+2,i) )
					
					! Fluxes accessing no ghost cell
					DO j=4,Ny-2
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )    &
							-coeff*ABS(Vadv(j,i,thread_nr))*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) &
								+ DBLE(5.0)*Q(j+1,i) - Q(j+2,i) )
					END DO
					
					! Fluxes accessing lower ghost cells
					j=Ny-1
					FluxVer(j,i,F_field_nr, thread_nr) = FluxVer(j,i,F_field_nr, thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr) ) &
							-coeff*ABS(Vadv(j,i,thread_nr))*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) &
								+ DBLE(5.0)*Q(j+1,i) - GhostDown(1,i,field_nr,thread_nr) )
					j=Ny
					FluxVer(j,i,F_field_nr, thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) & 
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*GhostDown(1,i,field_nr,thread_nr) &
						    + weights(1)*GhostDown(2,i,field_nr,thread_nr) ) &
							-coeff*ABS(Vadv(j,i,thread_nr))*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) &
								+ DBLE(5.0)*GhostDown(1,i,field_nr,thread_nr) - GhostDown(2,i,field_nr,thread_nr) )
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*GhostDown(1,i,field_nr,thread_nr) &
						    + weights(2)*GhostDown(2,i,field_nr,thread_nr) + weights(1)*GhostDown(3,i,field_nr,thread_nr) ) &
							-coeff*ABS(Vadv(j,i,thread_nr))*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) &
								- DBLE(10.0)*GhostDown(1,i,field_nr,thread_nr) + DBLE(5.0)*GhostDown(2,i,field_nr,thread_nr) &
								- GhostDown(3,i,field_nr,thread_nr) )
				
				END DO				
			
			CASE (6)
				
				weights(1) =  1.0
				weights(2) = -8.0
				weights(3) = 37.0
				coeff = DBLE(1.0)/DBLE(60.0)
								
				DO i=1,Nx
					! Fluxes accessing upper ghost cells
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(3,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(2,i,field_nr,thread_nr) + weights(3)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
						
					j=2
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(2,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr,thread_nr) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) &
							+ weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
					j=3
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-2,i) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
					
					! Fluxes accessing no ghost cell
					DO j=4,Ny-2
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
					END DO
					
					! Fluxes accessing lower ghost cells
					j=Ny-1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr) )
					j=Ny
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) & 
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*GhostDown(1,i,field_nr,thread_nr) &
						    + weights(1)*GhostDown(2,i,field_nr,thread_nr) )				
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*Vadv(j,i,thread_nr)*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*GhostDown(1,i,field_nr,thread_nr) &
						    + weights(2)*GhostDown(2,i,field_nr,thread_nr) + weights(1)*GhostDown(3,i,field_nr,thread_nr) )
				
				END DO
				
			CASE DEFAULT
				WRITE(*,*) 'No implementation available for flux of selected order. Now exiting.'
				STOP			
		END SELECT
	END SUBROUTINE UpdateVerticalAdvectionFlux
	
	SUBROUTINE FillinGhostcells(Qleft, Qright, Qup, Qdown, Qupleft, Qupright, Qdownleft, Qdownright, Nx, Ny, Nghost)
		DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Qleft, Qright, Qup, Qdown, Qupleft, Qupright, Qdownleft, Qdownright
		INTEGER, INTENT(IN) :: Nx, Ny, Nghost

		INTEGER :: thread_nr, counter, i, j, k
		
		thread_nr = omp_get_thread_num()

		IF (param%mpi_init_thread_flag==1) THEN
		
			GhostCorners(1,:,0:param%nthreads-1) = RESHAPE( Qupleft(   1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			GhostCorners(2,:,0:param%nthreads-1) = RESHAPE( Qupright(  1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			GhostCorners(3,:,0:param%nthreads-1) = RESHAPE( Qdownright(1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			GhostCorners(4,:,0:param%nthreads-1) = RESHAPE( Qdownleft( 1:nr_fields*param%nthreads), (/ nr_fields, param%nthreads /) )
			
			! In MPI_THREAD_FUNNELED mode, this routine is only called by the master thread
			counter=1
			DO thread_nr=0,param%nthreads-1
				DO k=1,nr_fields
					DO i=1,Nghost
						DO j=1,Ny
							GhostLeft(j,Nghost-i+1,k,thread_nr)  = Qleft( counter)
							GhostRight(j,i,k,thread_nr)          = Qright(counter)
							counter=counter+1
						END DO
					END DO
				END DO
			END DO
			
			counter=1
			DO thread_nr=0,param%nthreads-1
				DO k=1,nr_fields
					DO i=1,Nx
						DO j=1,Nghost
							GhostUp(Nghost-j+1,i,k,thread_nr)   = Qup(  counter)
							GhostDown(j,i,k,thread_nr)          = Qdown(counter)
							counter=counter+1
						END DO
					END DO
				END DO	
			END DO	

		ELSE IF ( (param%mpi_init_thread_flag==2) .or. (param%mpi_init_thread_flag==3) ) THEN
			! NEEDS TO BE IMPLEMENTED
		END IF
				
		IF (param%echo_on) WRITE(*,'(A)') ' Communication completed, ghost-cell buffers updated... '
				
	END SUBROUTINE FillinGhostcells
	
	SUBROUTINE GetMpiDatatypePar(Nx, Ny, Nghost, blocklengths, indices, length, length_singleThread)
		
		INTEGER,                 INTENT(IN)  :: Nx, Ny, Nghost
		INTEGER, DIMENSION(:,:), INTENT(OUT) :: blocklengths, indices
		INTEGER, DIMENSION(:),   INTENT(OUT) :: length, length_singleThread
		
		INTEGER :: i
			
		! --- Upper and lower ghost-cells ---
		length(2)              = Nx*nr_fields*param%Nthreads
		length_singleThread(2) = Nx*nr_fields
		length(7)              = Nx*nr_fields*param%Nthreads
		length_singleThread(7) = Nx*nr_fields
		
		! Upper row datatype: For every column and every field, one block of Nghost values has to be defined								
		blocklengths(1:length(2),2) = Nghost
		DO i=1,length(2)
			indices(i,2) = (i-1)*Ny
		END DO		

		! Lower row datatype
		blocklengths(1:length(7),7) = Nghost	
		DO i=1,length(7)
			indices(i,7) = i*Ny-Nghost
		END DO
		
		! --- Left and right ghost-cells ---
		length(4)              = nr_fields*param%Nthreads
		length_singleThread(4) = nr_fields
		length(5)              = nr_fields*param%Nthreads
		length_singleThread(5) = nr_fields
		
		! Left column datatype
		blocklengths(1:length(4),4) = Nghost*Ny
		DO i=1,length(4)
			indices(i,4) = (i-1)*Nx*Ny
		END DO

		! Right column datatype
		blocklengths(1:length(5),5) = Nghost*Ny
		DO i=1,length(5)
			indices(i,5) = (Nx-Nghost)*Ny + (i-1)*Nx*Ny
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
		DO k=1,nr_fields
			DO i=1,Nx
				DO j=1,Ny
					Q(j,i,k) = Y(counter)
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
		
		DO k=1,nr_fields
			DO i=1,Nx
				DO j=1,Ny
					Y(counter) = Q(j,i,k)
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
		
		Ny = SIZE(Q, 1)
		Nx = SIZE(Q, 2)
		
		SELECT CASE (BC)
		
			! Periodic in all four directions
			CASE (1)
			
				! Iterate over the different fields of the solution.
				DO k=1,nr_fields
			
					! Fill horizontal ghost cells
					DO i=1,Nghost
						DO j=1,Ny
							GhostLeft( j,i,k,thread_nr) = Q(j, Nx-i+1, k)
							GhostRight(j,i,k,thread_nr) = Q(j,i, k)
						END DO
					END DO
					
					! Fill vertical ghost cells
					DO i=1,Nx
						DO j=1,Nghost
							GhostUp(  j,i,k,thread_nr) = Q(Ny-j+1,i, k)
							GhostDown(j,i,k,thread_nr) = Q(j,i, k)
						END DO
					END DO
				
				END DO
				
				! Upper left corner
				GhostCorners(1,:,thread_nr) = Q(Ny,Nx,:)
				
				! Upper right corner
				GhostCorners(2,:,thread_nr) = Q(Ny,1,:)
				
				! Lower right corner
				GhostCorners(3,:,thread_nr) = Q(1,1,:)
				
				! Lower left corner
				GhostCorners(4,:,thread_nr) = Q(1,Nx,:)
				
			! Outflow boundary conditions in all four direction
			CASE (2)
			
				! Iterate over the different fields of the solution
				DO k=1,nr_fields
			
					! Fill horizontal ghost cells
					DO i=1,Nghost
						DO j=1,Ny
							GhostLeft( j,i,k,thread_nr)  = Q(j,1,k)
							GhostRight(j,i,k,thread_nr)  = Q(j, Nx, k)
						END DO
					END DO
				
					! Fill vertical ghost cells	
					DO i=1,Nx
						DO j=1,Nghost
							GhostUp(  j,i,k,thread_nr) = Q(1,i,k)
							GhostDOwn(j,i,k,thread_nr) = Q(Ny, i,k)
						END DO
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