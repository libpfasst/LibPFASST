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
INTEGER, PARAMETER :: nr_fields = 1, buffer_layout = 0

TYPE fdm_parameter
	INTEGER :: Nthreads, mpi_init_thread_flag
	DOUBLE PRECISION :: c_s, stabFreq
	LOGICAL :: echo_on
END TYPE

TYPE(fdm_parameter) :: param

PRIVATE
PUBLIC GetRHS, GetAcousticAdvectionRHS, GetBoussinesqRHS, InitializeFiniteVolumes, nr_fields, GetMpiDatatypePar, FillinGhostcells, &
	PackSolution, UnpackSolution, buffer_layout, FillGhostCells

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
		INTEGER :: thread_nr, u_field_nr, v_field_nr, Nx, Ny, i, j
					
			! If Uadv and Vadv contain only a single number, these numbers are the indices of the solution's velocity components (usually 1 and 2)
		u_field_nr = 1
		v_field_nr = 1
		CALL UpdateAdvectionVelocityNonlinear( Q(:,:,u_field_nr), Q(:,:,v_field_nr), u_field_nr, v_field_nr)
		!ELSE
		!	CALL UpdateLinearAdvectionVelocity(Uadv, Vadv)
		!END IF

		thread_nr = omp_get_thread_num()
		
		Ny = SIZE(Q,1)
		Nx = SIZE(Q,2)

		RQ(:,:,1) = DBLE(0.0)
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
	
	! ---- Acoustic-advection ----------------|
	! The flags "fluxes" and "fields" allow to determine which fields are updated with
	! which physical effects.
	!
	! F, F : Full RHS; tendencies of velocity and pressure are evaluated from advection, acoustics and divergence damping
	! F, V : returns tendencies of velocity only, tendencies of pressure are set to zero (useful for forward-backward integration, e.g.)
	! F, P : returns tendencies of pressure only, tendenccies of velocity are set to zero
	! A, V : tendencies of all fields due to advection (useful for operator splitting, e.g.)
	! A, S : tenencies for all fields due to acoustics and divergence damping
	SUBROUTINE GetAcousticAdvectionRHS(Uadv, Vadv, Q, order_advection, order_sound, RQ, dx, dy, nu, fluxes, fields)
	
		DOUBLE PRECISION, DIMENSION(:,:),   INTENT(IN)  :: Uadv, Vadv ! see comment in GetAdvectionRHS
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		INTEGER,                            INTENT(IN)  :: order_advection, &!
		                                                   order_sound       !
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: nu, dx, dy
		CHARACTER(len=1),                   INTENT(IN)  :: fluxes, &! F = all fluxes, A = advection fluxes only, S = sound fluxes only
														   fields   ! F = all fields, V = velocity only, P = pressure only
																	
																	! New possibility: fluxes=C, provides acoustic fluxes hardcoded 2nd order stencils. Fast, as no
																	! flux divergence is computed.
		INTEGER :: thread_nr, i, i1, i2, u_field_nr, v_field_nr
		
		! Components of Q and their physical interpretation
		! Q(:,:,1) = horizontal velocity u
		! Q(:,:,2) = horizontal velocity v
		! Q(:,:,3) = pressure p
		
		thread_nr = -1 ! Dummy line to avoid problems with TAU
		
		thread_nr = omp_get_thread_num()
	
		! Initialize RHS
		RQ(:,:,:) = DBLE(0.0)
						
		IF (fields=='F') THEN
			i1 = 1
			i2 = 3
		ELSEIF (fields=='V') THEN
			i1 = 1
			i2 = 2
		ELSEIF (fields=='P') THEN
			i1 = 3
			i2 = 3
		ELSE
			! This should not happen
			i1 = -1
			i2 = -1
			PRINT *,'Unknown entry for flag <fields> !'
		END IF
		
		! Initialize buffers storing contributions to fluxes
		IF (fluxes=='C') THEN
			! Do nothing
		ELSE
			FluxHor(:,:,i1:i2,thread_nr) = DBLE(0.0)
			FluxVer(:,:,i1:i2,thread_nr) = DBLE(0.0)
		END IF	
			
		IF (fluxes == 'F' .OR. fluxes == 'A') THEN
			
			IF ( SIZE(Uadv, 1) == 1 ) THEN
				u_field_nr = INT(Uadv(1,1))
				v_field_nr = INT(Vadv(1,1))
				CALL UpdateAdvectionVelocityNonlinear( Q(:,:, u_field_nr), Q(:,:,v_field_nr), u_field_nr, v_field_nr)
			ELSE
				CALL UpdateLinearAdvectionVelocity(Uadv, Vadv)
			END IF
					
			!# Compute Advective fluxes
			DO i=i1,i2
				! Add flux from horizontal advection
				CALL UpdateHorizontalAdvectionFlux( i, Q(:,:,i), i, order_advection)
				
				! Add flux from vertical advection
				CALL UpdateVerticalAdvectionFlux(   i, Q(:,:,i), i, order_advection)
												
			END DO	
		
		ELSEIF (fluxes=='S') THEN
			! If only update of sound fluxes is required, no advection update is necessary
			
		ELSEIF (fluxes=='C') THEN	
		
			CALL FillCellDivergence(Q(:,:,1), Q(:,:,2), 1, 2, dx, dy)

			! Provides coarse but fast evaluation of acoustic fluxes
			CALL FastModesAcousticAdvection( Q(:,:,:), RQ(:,:,:), fields, dx, dy, nu)
			
		ELSE
			WRITE(*,*) 'Unknown entry for flag <fluxes> !'									
		END IF
				
		! If full fluxes are requested, compute all fluxes with set order in conservation form		
		IF (fluxes=='F' .OR. fluxes=='S') THEN
		
			CALL FillCellDivergence(Q(:,:,1), Q(:,:,2), 1, 2, dx, dy)
		
			! # Second : Acoustic fluxes and divergence damping for the velocity
			IF (fields=='F' .OR. fields=='V') THEN
			
				! Fill velocity tendencies if either "full" fields (F) or "velocity" (V) is requested
							
				! Acoustic fluxes
				CALL UpdateHorizontalFlux(param%c_s, 1, Q(:,:,3), 3, order_sound)
				! Fluxes from divergence damping	
				CALL UpdateHorizontalFlux(-nu, 1, Divergence(:,:,thread_nr), nr_fields+1, order_sound)
				
				! Add contribution from acoustic flux
				CALL UpdateVerticalFlux(param%c_s, 2, Q(:,:,3), 3, order_sound)	
				! Add contribution from divergence damping				
				CALL UpdateVerticalFlux( -nu, 2, Divergence(:,:,thread_nr), nr_fields+1, order_sound)
				! Compute flux divergence and add to tendencies

			END IF
			
			! # Third : Acoustic flux for the pressure
			IF (fields=='F' .OR. fields=='P') THEN
			
				! Fill pressure tendencies if either "full" fields (F) or "pressure" (P) is requested
				
				CALL UpdateHorizontalFlux(param%c_s, 3, Q(:,:,1), 1, order_sound)
				CALL UpdateVerticalFlux(  param%c_s, 3, Q(:,:,2), 2, order_sound)
				
			END IF
					
		END IF
		
		! Finally, update all RHS with flux divergences
		IF (fluxes=='C') THEN
			! For coarse acoustic fluxes, no flux divergence has to be computed			
		ELSE
			CALL UpdateFluxDivergence(i1, i2, RQ, dx, dy)
		END IF
		
	END SUBROUTINE GetAcousticAdvectionRHS
	
	! ----- Linear, compressible Boussinesq equations -----------------------------------------|
	SUBROUTINE GetBoussinesqRHS(Uadv, Vadv, Q, order_advection, order_sound, RQ, dx, dy, nu, fluxes, fields)
	
		DOUBLE PRECISION, DIMENSION(:,:),   INTENT(IN)  :: Uadv, Vadv ! see comments in GetAdvectionRHS
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		INTEGER,                            INTENT(IN)  :: order_advection, &!
		                                                   order_sound       !
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: nu, dx, dy
		CHARACTER(len=1),                   INTENT(IN)  :: fluxes, &! F = all fluxes, A = advection fluxes only, S = sound and buoyancy fluxes only
														   fields   ! F = all fields, V = velocity only, P = pressure and buoyancy only	
	
		INTEGER :: thread_nr, i, i1, i2, u_field_nr, v_field_nr
		
		!> Components of Q and their physical interpretation:
		!> Q(:,:,1) = horizontal velocity u
		!> Q(:,:,2) = vertical velocity w
		!> Q(:,:,3) = buoyancy b
		!> Q(:,:,4) = pressure p
			
		thread_nr = -1 ! Dummy line to avoid problems with TAU		
			
		thread_nr = omp_get_thread_num()

		IF (fields=='F') THEN
			i1 = 1
			i2 = 4
		ELSEIF (fields=='V') THEN
			i1 = 1
			i2 = 2
		ELSEIF (fields=='P') THEN
			i1 = 3
			i2 = 4
		ELSE
			i1 = -1
			i2 = -1
			PRINT *,'Unknown entry for flag <fields> !'
		END IF
				
		! Initialize RHS
		RQ(:,:,:)                    = DBLE(0.0)
		
		IF (fluxes=='C') THEN
			! No initialization of flux buffers required
		ELSE
			FluxHor(:,:,i1:i2,thread_nr) = DBLE(0.0)
			FluxVer(:,:,i1:i2,thread_nr) = DBLE(0.0)						
		END IF
		
		IF (fluxes == 'F' .OR. fluxes == 'A') THEN
		
			IF ( SIZE(Uadv, 1) == 1 ) THEN
				u_field_nr = INT(Uadv(1,1))
				v_field_nr = INT(Vadv(1,1))
				CALL UpdateAdvectionVelocityNonlinear( Q(:,:,u_field_nr), Q(:,:,v_field_nr), u_field_nr, v_field_nr)
			ELSE
				CALL UpdateLinearAdvectionVelocity(Uadv, Vadv)
			END IF		
						
			!# Compute Advective fluxes
			DO i=i1,i2

				CALL UpdateHorizontalAdvectionFlux( i, Q(:,:,i), i, order_advection)
				CALL UpdateVerticalAdvectionFlux(   i, Q(:,:,i), i, order_advection)

			END DO	
		
		ELSEIF (fluxes=='S') THEN
			! If only update of sound fluxes is required, no advection update is necessary
			
		ELSEIF (fluxes=='C') THEN
		
			CALL FillCellDivergence(Q(:,:,1), Q(:,:,2), 1, 2, dx, dy)

			! Provides coarse but fast evaluation of acoustic fluxes
			CALL FastModesAcousticAdvection( Q(:,:,:), RQ(:,:,:), fields, dx, dy, nu)
						
		ELSE
			WRITE(*,*) 'Unknown entry for flag <fluxes> !'									
		END IF
				
		IF (fluxes=='F' .OR. fluxes=='S') THEN
		
			! Compute divergence of velocity field, required for applying divergence damping 
			CALL FillCellDivergence(Q(:,:,1), Q(:,:,2), 1, 2, dx, dy)
		
			! # Second : Acoustic and buoyancy fluxes and divergence damping			
			IF (fields=='F' .OR. fields=='V') THEN
			
				! Fill velocity tendencies if either "full" fields (F) or "velocity" (V) is requested
				
				! Acoustic fluxes
				CALL UpdateHorizontalFlux(param%c_s, 1, Q(:,:,4), 4, order_sound)
				! Add contribution from acoustic flux
				CALL UpdateVerticalFlux(  param%c_s, 2, Q(:,:,4), 4, order_sound)

				! Fluxes from divergence damping	
				CALL UpdateHorizontalFlux(-nu, 1, Divergence(:,:,thread_nr), nr_fields+1, 2)
				
				! Add contribution from divergence damping				
				CALL UpdateVerticalFlux(  -nu, 2, Divergence(:,:,thread_nr), nr_fields+1, 2)
				
				! Bouyancy flux
				RQ(:,:,2) = Q(:,:,3)
			
			END IF
			
			IF (fields=='F' .OR. fields=='P') THEN
			
				! Fill pressure and buoyancy tendencies if either "full" fields (F) or "pressure" (P) is requested
				RQ(:,:,3) = -param%stabFreq*param%stabFreq*Q(:,:,2)
				
				CALL UpdateHorizontalFlux(param%c_s, 4, Q(:,:,1), 1, order_sound)
				CALL UpdateVerticalFlux(  param%c_s, 4, Q(:,:,2), 2, order_sound)
				
			END IF
			
		END IF
		
		! Finally, update all RHS with flux divergences
		IF (fluxes=='C') THEN
			! Do nothing
		ELSE
			CALL UpdateFluxDivergence(i1, i2, RQ, dx, dy)
		END IF
		
	END SUBROUTINE GetBoussinesqRHS
	
	! ===== Initialization Routines ===== |
	
	SUBROUTINE InitializeFiniteVolumes(Nx_max, Ny_max, Nthreads, mpi_init_thread, echo_on, c_s, stabFreq)
		
		INTEGER, INTENT(IN) :: Nx_max,     &!
							   Ny_max,	   &!
							   Nthreads, mpi_init_thread
							   
		LOGICAL, INTENT(IN) :: echo_on
		DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: c_s, stabFreq
							   	
		INTEGER :: i, thread_nr
		
		param%echo_on = echo_on
		param%Nthreads = Nthreads
		
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
		
		IF (present(c_s)) param%c_s = c_s
		IF (present(stabFreq)) param%stabFreq = stabFreq
	
		IF (param%echo_on) WRITE(*, '(A, I2)') ' Module FiniteVolumes successfully initialized. Solution components : ', nr_fields
	
	END SUBROUTINE InitializeFiniteVolumes
	
	! ----------------------------------------------------------------
	! For linear advection, Ny x Nx+1 interface values for the horizontal
	! advection speed and Ny+1 x Nx interface values for the vertical
	! advection speed have to be prescribed.
	!
	! This subroutine is called only once, during initialization.
	! Nonlinear advection needs an update of the interface values by
	! linear interpolation for every evaluation of the flux functions.
	SUBROUTINE UpdateLinearAdvectionVelocity(U_interfaces, V_interfaces)
	
		DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: U_interfaces
		DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: V_interfaces
		
		INTEGER thread_nr, Nx, Ny
		
		Ny = SIZE(U_interfaces, 1)
		Nx = SIZE(V_interfaces, 2)
				
		thread_nr = omp_get_thread_num()		

		Uadv(1:Ny,1:Nx+1, thread_nr) = U_interfaces
		Vadv(1:Ny+1,1:Nx, thread_nr) = V_interfaces

				
	END SUBROUTINE UpdateLinearAdvectionVelocity
	
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
	
	!-----------------------------------------------------|
	! Update non-advective horizontal flux; no upwinding. |
	SUBROUTINE UpdateHorizontalFlux(alpha, F_field_nr, Q, field_nr, order)
	
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Q           ! Array of Ny x Nx cell values of advected quantity
			INTEGER,						  INTENT(IN) :: F_field_nr,&! The number of the field for which the current flux is to be added to the overall tendencies
															field_nr,  &! The number of the field that Q represents. Required to access the correct ghost cell values.
															order		! Order of accuracy of flux computation.
			DOUBLE PRECISION, INTENT(IN) :: alpha
			
			DOUBLE PRECISION, DIMENSION(3) :: weights
			DOUBLE PRECISION :: coeff
			INTEGER :: i, j, thread_nr, Nx, Ny
					
			Ny = SIZE(Q, 1)
			Nx = SIZE(Q, 2)
			thread_nr = omp_get_thread_num()
			
			SELECT CASE (order)
			
				CASE (1)
				
					weights(1) = 1.0
					coeff = alpha*DBLE(1.0)/DBLE(2.0)
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr, thread_nr) = FluxHor(j,i,F_field_nr, thread_nr) + coeff*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) + weights(1)*Q(j,i) ) &
							- ABS(coeff)*( Q(j,i) - GhostLeft(j,1,field_nr,thread_nr) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=2,Nx
						DO j=1,Ny
							FluxHor(j,i,F_field_nr, thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-1) + weights(1)*Q(j,i) ) &
								- ABS(coeff)*( Q(j,i) - Q(j,i-1) )
						END DO
					END DO
					
					! Fluxes accessing a ghost cell to the right
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr, thread_nr) = FluxHor(j,i,F_field_nr, thread_nr) + coeff*( weights(1)*Q(j,i-1) &
							+ weights(1)*GhostRight(j,1,field_nr,thread_nr) ) - ABS(coeff)*( GhostRight(j,1,field_nr,thread_nr) - Q(j,i-1) )
					END DO
					
				CASE (2)
				
					weights(1) = 1.0
					coeff = alpha*DBLE(1.0)/DBLE(2.0)
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr, thread_nr) = FluxHor(j,i,F_field_nr, thread_nr) + coeff*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) + weights(1)*Q(j,i) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=2,Nx
						DO j=1,Ny
							FluxHor(j,i,F_field_nr, thread_nr) = FluxHor(j,i,F_field_nr, thread_nr) + coeff*( weights(1)*Q(j,i-1) + weights(1)*Q(j,i) )
						END DO
					END DO
					
					! Fluxes accessing a ghost cell to the right
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-1) &
							+ weights(1)*GhostRight(j,1,field_nr,thread_nr) )
					END DO
									
				CASE (3)
				
					weights(1) = -1.0
					weights(2) =  7.0
					coeff = alpha*DBLE(1.0)/DBLE(12.0)
					
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) ) &
							- coeff*( -Q(j,i+1) + DBLE(3.0)*Q(j,i) - DBLE(3.0)*GhostLeft(j,1,field_nr,thread_nr) + GhostLeft(j,2,field_nr,thread_nr) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) + weights(2)*Q(j,i-1) &
							+ weights(2)*Q(j,i) + weights(1)*Q(j,i+1) ) &
							- coeff*( -Q(j,i+1) + DBLE(3.0)*Q(j,i) - DBLE(3.0)*Q(j,i-1) + GhostLeft(j,1,field_nr,thread_nr) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=3,Nx-1
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) &
								+ weights(2)*Q(j,i) + weights(1)*Q(j,i+1) ) - coeff*( -Q(j,i+1) + DBLE(3.0)*Q(j,i) - DBLE(3.0)*Q(j,i-1) + Q(j,i-2) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) + weights(2)*Q(j,i) &
							+ weights(1)*GhostRight(j,1,field_nr,thread_nr)) &
							- coeff*( -GhostRight(j,1,field_nr,thread_nr) + DBLE(3.0)*Q(j,i) - DBLE(3.0)*Q(j,i-1) + Q(j,i-2) )
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr, thread_nr) + coeff*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) + &
							weights(2)*GhostRight(j,1,field_nr,thread_nr) + weights(1)*GhostRight(j,2,field_nr,thread_nr) ) &
							- coeff*( -GhostRight(j,2,field_nr,thread_nr) + DBLE(3.0)*GhostRight(j,1,field_nr,thread_nr) - DBLE(3.0)*Q(j,i-1) &
								+ Q(j,i-2) )
					END DO
					
				CASE (4)

					weights(1) = -1.0
					weights(2) =  7.0
					coeff = alpha*DBLE(1.0)/DBLE(12.0)
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) + weights(2)*Q(j,i-1) &
							+ weights(2)*Q(j,i) + weights(1)*Q(j,i+1) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=3,Nx-1
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-2) &
								+ weights(2)*Q(j,i-1) + weights(2)*Q(j,i) + weights(1)*Q(j,i+1) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) + weights(2)*Q(j,i) &
							+ weights(1)*GhostRight(j,1,field_nr,thread_nr))
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-2) + weights(2)*Q(j,i-1) + &
							weights(2)*GhostRight(j,1,field_nr,thread_nr) + weights(1)*GhostRight(j,2,field_nr,thread_nr) )
					END DO
					
				CASE (5)
				
					weights(1) = 1.0
					weights(2) = -8.0
					weights(3) = 37.0
					coeff = alpha*DBLE(1.0)/DBLE(60.0)
					
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,3,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(3)*GhostLeft(j,1,field_nr,thread_nr) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) ) &
							-coeff*( -GhostLeft(j,3,field_nr,thread_nr) + DBLE(5.0)*GhostLeft(j,2,field_nr,thread_nr) &
								- DBLE(10.0)*GhostLeft(j,1,field_nr,thread_nr) + DBLE(10.0)*Q(j,i) - DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) ) &
							-coeff*( -GhostLeft(j,2,field_nr,thread_nr) + DBLE(5.0)*GhostLeft(j,1,field_nr,thread_nr) - DBLE(10.0)*Q(j,i-1) &
								+ DBLE(10.0)*Q(j,i) - DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
					END DO
					i=3
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(2)*Q(j,i-2) + weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) ) &
							-coeff*( -GhostLeft(j,1,field_nr,thread_nr) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) &
								-DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=4,Nx-2
						DO j=1,Ny
							FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) ) &
								-coeff*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) - DBLE(5.0)*Q(j,i+1) + Q(j,i+2) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx-1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) + weights(3)*Q(j,i-1) &
								+ weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*GhostRight(j,1,field_nr,thread_nr) ) &
								-coeff*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) - DBLE(5.0)*Q(j,i+1) &
									+ GhostRight(j,1,field_nr, thread_nr) )
					END DO
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) + weights(3)*Q(j,i-1) &
								+ weights(3)*Q(j,i) + weights(2)*GhostRight(j,1,field_nr,thread_nr) + weights(1)*GhostRight(j,2,field_nr,thread_nr) )  &
								-coeff*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*Q(j,i) &
									- DBLE(5.0)*GhostRight(j,1,field_nr, thread_nr) + GhostRight(j,2,field_nr, thread_nr) )
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) + weights(3)*Q(j,i-1) &
							+ weights(3)*GhostRight(j,1,field_nr,thread_nr) + weights(2)*GhostRight(j,2,field_nr,thread_nr) &
							+ weights(1)*GhostRight(j,3,field_nr,thread_nr) ) &
							-coeff*( -Q(j,i-3) + DBLE(5.0)*Q(j,i-2) - DBLE(10.0)*Q(j,i-1) + DBLE(10.0)*GhostRight(j,1,field_nr, thread_nr) &
								- DBLE(5.0)*GhostRight(j,2,field_nr,thread_nr) + GhostRight(j,3,field_nr, thread_nr) )
					END DO					
												
				CASE (6)
				
					weights(1) = 1.0
					weights(2) = -8.0
					weights(3) = 37.0
					coeff = alpha*DBLE(1.0)/DBLE(60.0)
					
					! Fluxes accessing a ghost cell to the left
					i=1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,3,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(3)*GhostLeft(j,1,field_nr,thread_nr) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
					END DO
					i=2
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,2,field_nr,thread_nr) &
							+ weights(2)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
					END DO
					i=3
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostLeft(j,1,field_nr,thread_nr) &
							+ weights(2)*Q(j,i-2) + weights(3)*Q(j,i-1) &
							+ weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
					END DO
					
					! Fluxes not accessing any ghost cell
					DO i=4,Nx-2
						DO j=1,Ny
							FluxHor(j,i,F_field_nr, thread_nr) = FluxHor(j,i,F_field_nr, thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) + weights(2)*Q(j,i+1) + weights(1)*Q(j,i+2) )
						END DO
					END DO
					
					! Fluxes accessing ghost cells to the right
					i=Nx-1
					DO j=1,Ny
						FluxHor(j,i, F_field_nr, thread_nr) = FluxHor(j,i, F_field_nr, thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) &
								+ weights(2)*Q(j,i+1) + weights(1)*GhostRight(j,1,field_nr,thread_nr) )
					END DO
					i=Nx
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*Q(j,i) &
								+ weights(2)*GhostRight(j,1,field_nr,thread_nr) + weights(1)*GhostRight(j,2,field_nr,thread_nr) )
					END DO
					i=Nx+1
					DO j=1,Ny
						FluxHor(j,i,F_field_nr,thread_nr) = FluxHor(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j,i-3) + weights(2)*Q(j,i-2) &
								+ weights(3)*Q(j,i-1) + weights(3)*GhostRight(j,1,field_nr,thread_nr) &
								+ weights(2)*GhostRight(j,2,field_nr,thread_nr) + weights(1)*GhostRight(j,3,field_nr,thread_nr) )
					END DO
			
				CASE DEFAULT
					WRITE(*,*) 'No implementation available for flux of selected order. Now exiting.'
					STOP
			END SELECT
									
	END SUBROUTINE UpdateHorizontalFlux

	! --------------------------------------------------
	! Update non-advective vertical flux; no upwinding |
	! Here, alpha is a constant, scalar, positive(??)
	! coefficient. 
	SUBROUTINE UpdateVerticalFlux(alpha, F_field_nr, Q, field_nr, order)
	
		DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Q			 ! Cell values of the solution component for which the flux is to be computed.
		INTEGER,						  INTENT(IN) :: F_field_nr, &! Index of solution component to which tendency the computed flux contributes
														field_nr,   &! Index of solution component, required to access correct ghost cells.
														order        ! Order of approximation, can be 1, 2,3,4 or 6
		DOUBLE PRECISION,                 INTENT(IN) :: alpha
		
		DOUBLE PRECISION, DIMENSION(3) :: weights  ! Weights of the stencil, filled depending on the chosen order
		DOUBLE PRECISION               :: coeff    ! auxiliary buffer 
		INTEGER                        :: i, j, thread_nr, Nx, Ny
		
		Ny = SIZE(Q, 1)
		Nx = SIZE(Q, 2)
		thread_nr = omp_get_thread_num()
		
		SELECT CASE (order)
		
			CASE (1)
			
				weights(1) = 1.0
				coeff      = alpha*(DBLE(1.0)/DBLE(2.0))
				
				DO i=1,Nx
					! Fluxes accessing an upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(1,i,field_nr,thread_nr) + weights(1)*Q(j,i) )&
						-ABS(coeff)*( GhostUp(1,i,field_nr,thread_nr) - Q(j,i) )
					
					! Fluxes not accessing any ghost cell
					DO j=2,Ny
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-1,i) + weights(1)*Q(j,i) )&
							-ABS(coeff)*( Q(j-1,i) - Q(j,i) )
					END DO
					
					! Fluxes accessing a lower ghost cell
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-1,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr) ) &
						-ABS(coeff)*( Q(j-1,i) - GhostDown(1,i,field_nr,thread_nr) )
				END DO
			
			CASE (2)
			
				weights(1) = 1.0
				coeff      = alpha*(DBLE(1.0)/DBLE(2.0))
				DO i=1,Nx
					! Fluxes accessing an upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(1,i,field_nr,thread_nr) + weights(1)*Q(j,i) )				
					! Fluxes not accessing any ghost cell
					DO j=2,Ny
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-1,i) + weights(1)*Q(j,i) )
					END DO
					
					! Fluxes accessing a lower ghost cell
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-1,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr) )					
				END DO				
			
			CASE (3)
			
				weights(1) = -1.0
				weights(2) =  7.0
				coeff = alpha*DBLE(1.0)/DBLE(12.0)
				
				DO i=1,Nx
					! Fluxes accessing at least one upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(2,i,field_nr, thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr, thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) &
						- coeff*( -GhostUp(2,i,field_nr, thread_nr) + DBLE(3.0)*GhostUp(1,i,field_nr,thread_nr) &
							- DBLE(3.0)*Q(j,i) + Q(j+1,i) )
					j=2  
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-1,i) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i)) &
						- coeff*( -GhostUp(1,i,field_nr,thread_nr) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*Q(j,i) + Q(j+1,i) )
					
					! Fluxes accessing no ghost cell
					DO j=3,Ny-1
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
							+ weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) &
							- coeff*( -Q(j-2,i) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*Q(j,i) + Q(j+1,i) )
					END DO
					
					! Fluxes accessing at least one lower ghost cell
					j=Ny
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*Q(j,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr)) &
						- coeff*( -Q(j-2,i) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*Q(j,i) + GhostDown(1,i,field_nr,thread_nr))
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*GhostDown(1,i,field_nr,thread_nr) + weights(1)*GhostDown(2,i,field_nr,thread_nr) ) &
						- coeff*( -Q(j-2,i) + DBLE(3.0)*Q(j-1,i) - DBLE(3.0)*GhostDown(1,i,field_nr,thread_nr) &
						+ GhostDown(2,i,field_nr,thread_nr))
					
				END DO
			
			CASE (4)
			
				weights(1) = -1.0
				weights(2) =  7.0
				coeff = alpha*DBLE(1.0)/DBLE(12.0)
				
				DO i=1,Nx
					! Fluxes accessing at least one upper ghost cell
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(2,i,field_nr, thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr, thread_nr) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) 
					j=2  
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-1,i) + weights(2)*Q(j,i) + weights(1)*Q(j+1,i)) 					
					! Fluxes accessing no ghost cell
					DO j=3,Ny-1
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
							+ weights(2)*Q(j,i) + weights(1)*Q(j+1,i) ) 
					END DO
					
					! Fluxes accessing at least one lower ghost cell
					j=Ny
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*Q(j,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr)) 
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr, thread_nr) + coeff*( weights(1)*Q(j-2,i) + weights(2)*Q(j-1,i) &
						+ weights(2)*GhostDown(1,i,field_nr,thread_nr) + weights(1)*GhostDown(2,i,field_nr,thread_nr) )					
				END DO			
				
			CASE (5)
				
				weights(1) =  1.0
				weights(2) = -8.0
				weights(3) = 37.0
				coeff = alpha*DBLE(1.0)/DBLE(60.0)
								
				DO i=1,Nx
					! Fluxes accessing upper ghost cells
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(3,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(2,i,field_nr,thread_nr) + weights(3)*GhostUp(1,i,field_nr,thread_nr) + weights(3)*Q(j,i) &
						+ weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) ) - coeff*( GhostUp(3,i,field_nr,thread_nr) &
						- DBLE(5.0)*GhostUp(2,i,field_nr,thread_nr) + DBLE(10.0)*GhostUp(1,i,field_nr,thread_nr) &
						- DBLE(10.0)*Q(j,i) + DBLE(5.0)*Q(j+1,i) - Q(j+2,i) )
					j=2
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(2,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr,thread_nr) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) &
						+ weights(1)*Q(j+2,i) ) - coeff*( GhostUp(2,i,field_nr,thread_nr) - DBLE(5.0)*GhostUp(1,i,field_nr,thread_nr) &
						+ DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) + DBLE(5.0)*Q(j+1,i) - Q(j+2,i) )
					j=3
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-2,i) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) ) &
						-coeff*( GhostUp(1,i,field_nr,thread_nr) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) &
						+ DBLE(5.0)*Q(j+1,i) - Q(j+2,i) )
					
					! Fluxes accessing no ghost cell
					DO j=4,Ny-2
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) ) &
							-coeff*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) + DBLE(5.0)*Q(j+1,i) - Q(j+2,i) )
					END DO
					
					! Fluxes accessing lower ghost cells
					j=Ny-1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) + weights(3)*Q(j-1,i) &
						+ weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr) ) &
						-coeff*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) + DBLE(5.0)*Q(j+1,i) &
						- GhostDown(1,i,field_nr,thread_nr) )
						
					j=Ny
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) + weights(3)*Q(j-1,i) &
						+ weights(3)*Q(j,i) + weights(2)*GhostDown(1,i,field_nr,thread_nr) + weights(1)*GhostDown(2,i,field_nr,thread_nr) )	&
						-coeff*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*Q(j,i) &
						+ DBLE(5.0)*GhostDown(1,i,field_nr,thread_nr) - GhostDown(2,i,field_nr,thread_nr) )
						
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) + weights(3)*Q(j-1,i)  &
						+ weights(3)*GhostDown(1,i,field_nr,thread_nr) + weights(2)*GhostDown(2,i,field_nr,thread_nr) &
						+ weights(1)*GhostDown(3,i,field_nr,thread_nr) ) &
						-coeff*( Q(j-3,i) - DBLE(5.0)*Q(j-2,i) + DBLE(10.0)*Q(j-1,i) - DBLE(10.0)*GhostDown(1,i,field_nr,thread_nr) &
							+ DBLE(5.0)*GhostDown(2,i,field_nr,thread_nr) - GhostDown(3,i,field_nr,thread_nr) )
				
				END DO				
			
			CASE (6)
				
				weights(1) =  1.0
				weights(2) = -8.0
				weights(3) = 37.0
				coeff = alpha*DBLE(1.0)/DBLE(60.0)
								
				DO i=1,Nx
					! Fluxes accessing upper ghost cells
					j=1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(3,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(2,i,field_nr,thread_nr) + weights(3)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
						
					j=2
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(2,i,field_nr,thread_nr) &
						+ weights(2)*GhostUp(1,i,field_nr,thread_nr) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) &
							+ weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
					j=3
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*GhostUp(1,i,field_nr,thread_nr) &
						+ weights(2)*Q(j-2,i) + weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
					
					! Fluxes accessing no ghost cell
					DO j=4,Ny-2
						FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*Q(j+2,i) )
					END DO
					
					! Fluxes accessing lower ghost cells
					j=Ny-1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*Q(j+1,i) + weights(1)*GhostDown(1,i,field_nr,thread_nr) )
					j=Ny
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) & 
							+ weights(3)*Q(j-1,i) + weights(3)*Q(j,i) + weights(2)*GhostDown(1,i,field_nr,thread_nr) &
						    + weights(1)*GhostDown(2,i,field_nr,thread_nr) )				
					j=Ny+1
					FluxVer(j,i,F_field_nr,thread_nr) = FluxVer(j,i,F_field_nr,thread_nr) + coeff*( weights(1)*Q(j-3,i) + weights(2)*Q(j-2,i) &
							+ weights(3)*Q(j-1,i) + weights(3)*GhostDown(1,i,field_nr,thread_nr) &
						    + weights(2)*GhostDown(2,i,field_nr,thread_nr) + weights(1)*GhostDown(3,i,field_nr,thread_nr) )
				
				END DO
				
			CASE DEFAULT
				WRITE(*,*) 'No implementation available for flux of selected order. Now exiting.'
				STOP
		END SELECT
				
	END SUBROUTINE UpdateVerticalFlux
	
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
	
	! -----
	SUBROUTINE FastModesAcousticAdvection(Q, RQ, fields, dx, dy, nu)
		
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)    :: Q
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: RQ
		CHARACTER(len=1),                   INTENT(IN)    :: fields
		DOUBLE PRECISION,					INTENT(IN)    :: dx, dy, nu
		
		DOUBLE PRECISION :: coeff_x, coeff_y
		INTEGER          :: Nx, Ny, i, j, thread_nr, p_index
		
		Ny     = SIZE(Q,1)
		Nx     = SIZE(Q,2)
		coeff_x = (0.5/dx)
		coeff_y = (0.5/dy)
		
		p_index = -1
		IF (nr_fields==3) THEN
			p_index = 3
		ELSEIF (nr_fields==4) THEN
			p_index = 4
		END IF
		
		thread_nr = omp_get_thread_num()
	
		IF (fields=='V') THEN

			! ### First field ###	
			

			! Left cells
			i=1
			j=1
			RQ(j,i,1) = -param%c_s*coeff_x*( Q(j,i+1,p_index) - GhostLeft(j,1,p_index,thread_nr) ) + nu*coeff_x*( Divergence(j,i+1,thread_nr) - GhostLeft(j,1,nr_fields+1,thread_nr) )
			
			DO j=2,Ny-1
				RQ(j,i,1) = -param%c_s*coeff_x*( Q(j,i+1,p_index) - GhostLeft(j,1,p_index,thread_nr) ) + nu*coeff_x*( Divergence(j,i+1,thread_nr) - GhostLeft(j,1,nr_fields+1,thread_nr) )
			END DO
			
			i=1
			j=Ny
			RQ(j,i,1) = -param%c_s*coeff_x*( Q(j,i+1,p_index) - GhostLeft(j,1,p_index,thread_nr) ) + nu*coeff_x*( Divergence(j,i+1,thread_nr) - GhostLeft(j,1,nr_fields+1,thread_nr) )
			
			DO i=2,Nx-1
				
				j=1
				RQ(j,i,1) = -param%c_s*coeff_x*( Q(j,i+1,p_index) - Q(j,i-1,p_index) )               + nu*coeff_x*( Divergence(j,i+1,thread_nr) - Divergence(j,i-1,thread_nr) )
				
				DO j=2,Ny-1
					RQ(j,i,1) = -param%c_s*coeff_x*( Q(j,i+1,p_index) - Q(j,i-1,p_index) ) + nu*coeff_x*( Divergence(j,i+1,thread_nr) - Divergence(j,i-1,thread_nr) )
				END DO
				
				j=Ny
				RQ(j,i,1) = -param%c_s*coeff_x*( Q(j,i+1,p_index) - Q(j,i-1,p_index) )                 + nu*coeff_x*( Divergence(j,i+1,thread_nr) - Divergence(j,i-1,thread_nr) )
				
			END DO
			
			! Right cells
			i=Nx
			j=1
			RQ(j,i,1) = -param%c_s*coeff_x*( GhostRight(j,1,p_index,thread_nr) - Q(j,i-1,p_index) ) + nu*coeff_x*( GhostRight(j,1,nr_fields+1,thread_nr) - Divergence(j,i-1,thread_nr) )
			
			DO j=2,Ny-1
				RQ(j,i,1) = -param%c_s*coeff_x*( GhostRight(j,1,p_index,thread_nr) - Q(j,i-1,p_index) ) + nu*coeff_x*( GhostRight(j,1,nr_fields+1,thread_nr) - Divergence(j,i-1,thread_nr) )
			END DO
			
			j=Ny
			RQ(j,i,1) = -param%c_s*coeff_x*( GhostRight(j,1,p_index,thread_nr) - Q(j,i-1,p_index) ) + nu*coeff_x*( GhostRight(j,1,nr_fields+1,thread_nr) - Divergence(j,i-1,thread_nr) )
			
			! ### Second field ###
			
			i=1
			j=1
			RQ(j,i,2) = -param%c_s*coeff_y*( GhostUp(1,i,p_index,thread_nr) - Q(j+1,i,p_index) )   + nu*coeff_y*( GhostUp(j,1,nr_fields+1,thread_nr) - Divergence(j+1,i,thread_nr) )

			DO j=2,Ny-1
				RQ(j,i,2) = -param%c_s*coeff_y*( Q(j-1,i,p_index) - Q(j+1,i,p_index) ) + nu*coeff_y*( Divergence(j-1,i,thread_nr) - Divergence(j+1,i,thread_nr) )
			END DO
			
			i=1
			j=Ny
			RQ(j,i,2) = -param%c_s*coeff_y*( Q(j-1,i,p_index) - GhostDown(1,i,p_index,thread_nr) )      + nu*coeff_y*( Divergence(j-1,i,thread_nr) - GhostDown(1,i,nr_fields+1,thread_nr) )


			DO i=2,Nx-1
				j=1
				RQ(j,i,2) = -param%c_s*coeff_y*( GhostUp(1,i,p_index,thread_nr) - Q(j+1,i,p_index) ) + nu*coeff_y*( GhostUp(1,i,nr_fields+1,thread_nr) - Divergence(j+1,i,thread_nr) )
				
				DO j=2,Ny-1
					RQ(j,i,2) = -param%c_s*coeff_y*( Q(j-1,i,p_index) - Q(j+1,i,p_index) ) + nu*coeff_y*( Divergence(j-1,i,thread_nr) - Divergence(j+1,i,thread_nr) )
				END DO
			
				j=Ny
				RQ(j,i,2) = -param%c_s*coeff_y*( Q(j-1,i,p_index) - GhostDown(1,i,p_index,thread_nr) ) + nu*coeff_y*( Divergence(j-1,i,thread_nr) - GhostDown(1,i,nr_fields+1,thread_nr) )				
			END DO
			
			i=Nx
			j=1
			RQ(j,i,2) = -param%c_s*coeff_y*( GhostUp(1,i,p_index,thread_nr) - Q(j+1,i,p_index) )    + nu*coeff_y*( GhostUp(1,i,nr_fields+1,thread_nr) - Divergence(j+1,i,thread_nr) )

			DO j=2,Ny-1
				RQ(j,i,2) = -param%c_s*coeff_y*( Q(j-1,i,p_index) - Q(j+1,i,p_index) ) + nu*coeff_y*( Divergence(j-1,i,thread_nr) - Divergence(j+1,i,thread_nr) )			
			END DO
			
			j=Ny
			RQ(j,i,2) = -param%c_s*coeff_y*( Q(j-1,i,p_index) - GhostDown(1,i,p_index,thread_nr) )  + nu*coeff_y*( Divergence(j-1,i,thread_nr) - GhostDown(1,i,nr_fields+1,thread_nr) )

			IF (nr_fields==4) THEN
				! Add buoyancy induced tendency to vertical velocity tendency
				RQ(:,:,2) = RQ(:,:,2) + Q(:,:,3)
			END IF			
			
			
			
		ELSEIF (fields=='P') THEN
				
			! Left cells
		
			i=1
			j=1
			RQ(j,i,p_index) = -param%c_s*coeff_x*( Q(j,i+1,1) - GhostLeft(j,1,1,thread_nr) ) - param%c_s*coeff_y*( GhostUp(1,i,2,thread_nr) - Q(j+1,i, 2) )
			
			DO j=2,Ny-1
				RQ(j,i,p_index) = -param%c_s*coeff_x*( Q(j,i+1,1) - GhostLeft(j,1,1,thread_nr) ) - param%c_s*coeff_y*( Q(j-1,i,2) - Q(j+1,i,2) )
			END DO
			
			i=1
			j=Ny
			RQ(j,i,p_index) = -param%c_s*coeff_x*( Q(j,i+1,1) - GhostLeft(j,1,1,thread_nr) ) - param%c_s*coeff_y*(  Q(j-1,i,2) - GhostDown(1,i,2,thread_nr) )
			
			! Inner cells
		
			DO i=2,Nx-1
				
				j=1
				RQ(j,i,p_index) = -param%c_s*coeff_x*( Q(j,i+1,1) - Q(j,i-1,1) ) - param%c_s*coeff_y*( GhostUp(1,i,2,thread_nr) - Q(j+1,i,2) )
				
				DO j=2,Ny-1
					RQ(j,i,p_index) = -param%c_s*coeff_x*( Q(j,i+1,1) - Q(j,i-1,1) ) - param%c_s*coeff_y*( Q(j-1,i,2) - Q(j+1,i,2) )
				END DO
				
				j=Ny
				RQ(j,i,p_index) = -param%c_s*coeff_x*( Q(j,i+1,1) - Q(j,i-1,1) ) - param%c_s*coeff_y*( Q(j-1,i,2) - GhostDown(1,i,2,thread_nr) )
				
			END DO
			
			! Right cells
			i=Nx
			j=1
			RQ(j,i,p_index) = -param%c_s*coeff_x*( GhostRight(j,1,1,thread_nr) - Q(j,i-1,1) ) - param%c_s*coeff_y*( GhostUp(1,i,2,thread_nr) - Q(j+1,i,2) )
			
			DO j=2,Ny-1
				RQ(j,i,p_index) = -param%c_s*coeff_x*( GhostRight(j,1,1,thread_nr) - Q(j,i-1,1) ) - param%c_s*coeff_y*( Q(j-1,i,2) - Q(j+1,i,2) )
			END DO
			
			j=Ny
			RQ(j,i,p_index) = -param%c_s*coeff_x*( GhostRight(j,1,1,thread_nr) - Q(j,i-1,1) ) - param%c_s*coeff_y*( Q(j-1,i,2) - GhostDown(1,i,2,thread_nr) )
			
			IF (nr_fields==4) THEN
				! Stability
				RQ(:,:,3) = - param%stabFreq*param%stabFreq*Q(:,:,2)
			END IF					
										
		END IF
		
	END SUBROUTINE FastModesAcousticAdvection
	
	! -----------------------------------------------------------|
	! Compute cell centered values of the divergence u_x + v_z
	! using a second order, centered stencil. Values are written
	! into Divergence(:,:,thread_nr)
	!
	! Note that for the acoustic-advection equation, this limits
	! the accuracy of the acoustic terms to 2.
	SUBROUTINE FillCellDivergence(U, V, u_field_nr, v_field_nr, dx, dy)
		
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: U, &          ! Horizontal velocity field , Ny x Nx cell values
													        V             ! Vertical velocity field,    Ny x Nx cell values
			INTEGER,						  INTENT(IN) :: u_field_nr, & ! Index of U in the global array Q, that is U = Q(:,:,u_field_nr,thread_nr)
																		  ! Required to access correct ghost cells.
															v_field_nr    ! Index of V in the global array Q, that is V = Q(:,:,v_field_nr,thread_nr)
			DOUBLE PRECISION,                 INTENT(IN) :: dx, dy
			
			DOUBLE PRECISION :: coeff_x, coeff_y
			INTEGER :: i, j, thread_nr, Nx, Ny
			
			Ny = SIZE(V, 1)
			Nx = SIZE(U, 2)
			thread_nr = omp_get_thread_num()

			coeff_x	= DBLE((2*dx)**(-1))
			coeff_y = DBLE((2*dy)**(-1))
			
			! ### First, fill cells in leftmost column (hence i=1)
			i=1 
			j=1 ! Left upper cell
			Divergence(j,i,thread_nr) = coeff_x*( U(j,i+1) - GhostLeft(j,1,u_field_nr, thread_nr) ) &
				+ coeff_y*( GhostUp(1,i,v_field_nr, thread_nr) - V(j+1,i) )
	
			! other cells in leftmost column
			DO j=2,Ny-1
				Divergence(j,i,thread_nr) = coeff_x*( U(j,i+1) - GhostLeft(j,1,u_field_nr, thread_nr) ) + coeff_y*( V(j-1,i) - V(j+1,i) )
			END DO
			
			j=Ny ! left lower cell
			Divergence(j,i,thread_nr) = coeff_x*( U(j,i+1) - GhostLeft(j,1,u_field_nr, thread_nr) ) &
				+ coeff_y*( V(j-1,i) - GhostDown(1, 1, v_field_nr, thread_nr) )
			
			! ### Second, fill cells in columns between leftmost and rightmost column (hence i=2,..,Nx-1)
			DO i=2,Nx-1
				j=1 ! Uppermost cell in current column
				Divergence(j,i,thread_nr) = coeff_x*( U(j,i+1) - U(j,i-1) ) + coeff_y*( GhostUp(1,i,v_field_nr,thread_nr) - V(j+1,i) )
				
				! Inner cells in current column
				DO j=2,Ny-1
					Divergence(j,i,thread_nr) = coeff_x*( U(j,i+1) - U(j,i-1) ) + coeff_y*( V(j-1,i) - V(j+1,i) )
				END DO
				
				j=Ny ! Lowermost cell in current column
				Divergence(j,i,thread_nr) = coeff_x*( U(j,i+1) - U(j,i-1) ) + coeff_y*( V(j-1,i) - GhostDown(1,i,v_field_nr, thread_nr) )
				
			END DO
			
			! ### Third, fill cells in rightmost column (hence i=Nx)
			i=Nx
			
			j=1 ! Upper right cell
			Divergence(j,i,thread_nr) = coeff_x*( GhostRight(j,1,u_field_nr,thread_nr) - U(j,i-1) ) &
				+ coeff_y*( GhostUp(1,i,v_field_nr,thread_nr) - V(j+1,i) )
			
			! other cells in rightmost column
			DO j=2,Ny-1
				Divergence(j,i,thread_nr) = coeff_x*( GhostRight(j,1,u_field_nr, thread_nr) - U(j,i-1) ) + coeff_y*( V(j-1,i) - V(j+1,i) )
			END DO
			
			j=Ny ! Lower right cell
			Divergence(j,i,thread_nr) = coeff_x*( GhostRight(j,1,u_field_nr,thread_nr) - U(j,i-1) ) &
				+ coeff_y*( V(j-1,i) - GhostDown(1,i,v_field_nr, thread_nr) )
			
			! ### Last step: Fill one ghost cell in every direction for use in divergence damping term flux computation
			!	  Ghost cell values of divergence are stored in the ghost cell arrays at position field_nr + 1
			
			! As the divergence is only required to apply to anyhow artificially added divergence damping, the ghost cells
			! for the divergence are simply filled with zeros and no effort is made to obtain an accurate representation of the
			! divergence at the boundaries.
			
			! FOR NOW, SET CORNER GHOST CELLS TO ZERO
			!GhostCorners(:,:,thread_nr) = DBLE(0.0)
			
			j=1 ! Upper left and right ghost cell
			
			GhostLeft(1, j, nr_fields+1, thread_nr) =  coeff_x*( U(j,1) - GhostLeft(j,2,u_field_nr, thread_nr) ) &
				+ coeff_y*( GhostCorners(1,v_field_nr,thread_nr) - GhostLeft(j+1,1,v_field_nr,thread_nr) )
			
			GhostRight(1, j, nr_fields+1, thread_nr) = coeff_x*( GhostRight(j,2,u_field_nr,thread_nr) - U(j,Nx) ) &
				+ coeff_y*( GhostCorners(2,v_field_nr,thread_nr) - GhostRight(j+1,1,v_field_nr,thread_nr) )
				
			DO j=2,Ny-1 ! Middle left and right ghost cells
				
				GhostLeft(j, 1, nr_fields+1, thread_nr) = coeff_x*( U(j,1) - GhostLeft(j,2,u_field_nr,thread_nr) ) &
					+ coeff_y*( GhostLeft(j-1,1,v_field_nr,thread_nr) - GhostLeft(j+1,1,v_field_nr,thread_nr) )
				
				GhostRight(j, 1, nr_fields+1,thread_nr) = coeff_x*( GhostRight(j,2,u_field_nr,thread_nr) - U(j,Nx) ) &
					+ coeff_y*( GhostRight(j-1,1,v_field_nr,thread_nr) - GhostRight(j+1,1,v_field_nr,thread_nr) )
			END DO
			
			j=Ny ! Lower left and right ghost cell
			
			GhostLeft(j,1,nr_fields+1,thread_nr) = coeff_x*( U(j,1) - GhostLeft(j,2,u_field_nr,thread_nr) ) &
					+ coeff_y*( GhostLeft(j-1,1,v_field_nr,thread_nr) - GhostCorners(4,v_field_nr,thread_nr) )
			
			GhostRight(j,1,nr_fields+1,thread_nr) = coeff_x*( GhostRight(j,2,u_field_nr,thread_nr) - U(j,Nx) ) &
					+ coeff_y*( GhostRight(j-1,1,v_field_nr,thread_nr) - GhostCorners(3,v_field_nr,thread_nr) )
			
			! Divergence in upper and lower ghost cells.
			
			i=1
			
			GhostUp(1,i,nr_fields+1,thread_nr) = coeff_x*( GhostUp(1,i+1,u_field_nr,thread_nr) - GhostCorners(1,u_field_nr,thread_nr)) &
					+ coeff_y*( GhostUp(2,i,v_field_nr,thread_nr) - V(1,i) )
			
			GhostDown(1,i,nr_fields+1,thread_nr) = coeff_x*( GhostDown(1,i+1,u_field_nr,thread_nr) - GhostCorners(4,u_field_nr,thread_nr) ) &
					+ coeff_y*( V(Ny,i) - GhostDown(2,i,v_field_nr,thread_nr) )
			DO i=2,Nx-1
			
				GhostUp(1,i,nr_fields+1,thread_nr) = coeff_x*( GhostUp(1,i+1,u_field_nr,thread_nr) - GhostUp(1,i-1,u_field_nr,thread_nr) ) &
					+ coeff_y*( GhostUp(2,i,v_field_nr,thread_nr) - V(1,i) )
				
				GhostDown(1,i,nr_fields+1,thread_nr) = coeff_x*( GhostDown(1,i+1,u_field_nr,thread_nr) &
					- GhostDown(1,i-1,u_field_nr,thread_nr) ) + coeff_y*( V(Ny,i) - GhostDown(2,i,v_field_nr,thread_nr) )
			END DO
			
			i=Nx
			GhostUp(1,i,nr_fields+1,thread_nr) = coeff_x*( GhostCorners(2, u_field_nr,	thread_nr) - GhostUp(1,i-1,u_field_nr,thread_nr) ) &
					+ coeff_y*( GhostUp(2,i,v_field_nr,thread_nr) - V(1,i) )
					
			GhostDown(1,i,nr_fields+1,thread_nr) = coeff_x*( GhostCorners(3, u_field_nr, thread_nr) &
					- GhostDown(1,i-1,u_field_nr,thread_nr) ) + coeff_y*( V(Ny,i) - GhostDown(2,i,v_field_nr,thread_nr) )
			
			! TESTING: OVERWRITE ALL VALUES WITH ZEROS
			!GhostLeft( :, :, nr_fields+1, thread_nr) = DBLE(0.0)
			!GhostRight(:, :, nr_fields+1, thread_nr) = DBLE(0.0)
			!GhostUp(   :, :, nr_fields+1, thread_nr) = DBLE(0.0)
			!GhostDown( :, :, nr_fields+1, thread_nr) = DBLE(0.0)
			
	END SUBROUTINE FillCellDivergence
	
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