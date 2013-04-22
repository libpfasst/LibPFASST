MODULE WENO5_1D
! This module provides the implementation of a one-dimensional, finite difference WENO-5 scheme.
! Note that in order to ensure compatibility, all data structures do posses the same layout as
! the two-dimesional version, but the y-coordinate is always assumed to have only length one.
!
! Daniel Ruprecht, 19.1.2012
! ICS Lugano

!USE FVMParameters, only : nr_fields, BC, c_s, stabFreq, coriolisPar, grav
USE omp_lib,       only : omp_get_thread_num

IMPLICIT NONE

INTEGER, PARAMETER :: nr_fields = 1, buffer_layout = 1, stabFreq = 0.01, coriolisPar = 0.0, grav = 1.0

! Define buffers storing ghost-cell values. NOTE: In the 1-D module, GhostUp and GhostDown are only listed
! to ensure compatibility, they are neither used nor allocated.
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: GhostLeft, GhostRight, GhostUp, GhostDown, GhostFluxLeft, GhostFluxRight

! Define buffers storing the horizontal cell and interface flux values
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FluxInt_hor, FluxCell_hor

! Define fixed parameters used by the WENO-5 method (see Shu .... e.g.)
DOUBLE PRECISION, PARAMETER, DIMENSION(3) :: weights_plus = (/ 0.3, 0.6, 0.1 /)
DOUBLE PRECISION, PARAMETER, DIMENSION(5) :: stencil_weights = (/ 2.0, 5.0, -1.0, -7.0, 11.0 /)*(1.0/6.0)
DOUBLE PRECISION, PARAMETER               :: coeff_1 = 13.0/12.0, coeff_2 = 1.0/4.0
DOUBLE PRECISION, PARAMETER               :: weno_tol = 1.0e-6
INTEGER,          PARAMETER               :: weno_n   = 2

CONTAINS
	
	!
	! Nonlinear, rotating 1-D shallow water equations in conservation form with
	!
	! h = height
	! q = horizontal momentum
	! r = vertical momentum
	!
	! read
	!
	! h_t + f1(h,q,r) = 0
	! q_t + f2(h,q,r) = f_cor*r
	! r_t + f3(h,q,r) = -f_cor*q
	! 
	! with flux functions
	!
	! f1(h,q,r) = q
	! f2(h,q,r) = q^2/h + 0.5*g*h^2
	! f3(h,q,r) = q*r/h
	!
	SUBROUTINE GetRHS(Q, RQ, dt, dx, order)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: dt, dx
		INTEGER,                            INTENT(IN)  :: order
		
		INTEGER :: thread_nr
		
		thread_nr = omp_get_thread_num()
		
		IF (order==5) THEN
				
			! Flux function f1
			FluxCell_hor(  1,:,:,thread_nr) = Q(2,:,:)
			GhostFluxLeft( 1,:,:,thread_nr) = GhostLeft( 2,:,:,thread_nr)
			GhostFluxRight(1,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)
			
			! Flux function f2
			FluxCell_hor( 2,:,:,thread_nr)  = ( Q(2,:,:)*Q(2,:,:)/Q(1,:,:) ) + 0.5*grav*Q(1,:,:)*Q(1,:,:)
			
			GhostFluxLeft(2,:,:,thread_nr)  = ( GhostLeft(2,:,:,thread_nr)*GhostLeft(2,:,:,thread_nr)/GhostLeft(1,:,:,thread_nr) ) &
				+ 0.5*grav*GhostLeft(1,:,:,thread_nr)*GhostLeft(1,:,:,thread_nr)
				
			GhostFluxRight(2,:,:,thread_nr) = ( GhostRight(2,:,:,thread_nr)*GhostRight(2,:,:,thread_nr)/GhostRight(1,:,:,thread_nr) ) &
				+ 0.5*grav*GhostRight(1,:,:,thread_nr)*GhostRight(1,:,:,thread_nr)
				
			! FLux function f3	
			FluxCell_hor(3,:,:,thread_nr)   = Q(2,:,:)*Q(3,:,:)/Q(1,:,:)
			
			GhostFluxLeft(3,:,:,thread_nr)  = GhostLeft(2,:,:,thread_nr)*GhostLeft(3,:,:,thread_nr)/GhostLeft(1,:,:,thread_nr)
			
			GhostFluxRight(3,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)*GhostRight(3,:,:,thread_nr)/GhostRight(1,:,:,thread_nr)
		
			CALL UpdateHorizontalFlux(Q, dble(1.0))
		
		ELSE IF (order==1) THEN
					
			! Flux function f1
			FluxCell_hor(  1,:,:,thread_nr) = Q(2,:,:)
			GhostFluxLeft( 1,:,:,thread_nr) = GhostLeft( 2,:,:,thread_nr)
			GhostFluxRight(1,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)
			
			! Flux function f2
			FluxCell_hor( 2,:,:,thread_nr)  = ( Q(2,:,:)*Q(2,:,:)/Q(1,:,:) ) + 0.5*grav*Q(1,:,:)*Q(1,:,:)
			
			GhostFluxLeft(2,:,:,thread_nr)  = ( GhostLeft(2,:,:,thread_nr)*GhostLeft(2,:,:,thread_nr)/GhostLeft(1,:,:,thread_nr) ) &
				+ 0.5*grav*GhostLeft(1,:,:,thread_nr)*GhostLeft(1,:,:,thread_nr)
				
			GhostFluxRight(2,:,:,thread_nr) = ( GhostRight(2,:,:,thread_nr)*GhostRight(2,:,:,thread_nr)/GhostRight(1,:,:,thread_nr) ) &
				+ 0.5*grav*GhostRight(1,:,:,thread_nr)*GhostRight(1,:,:,thread_nr)
				
			! FLux function f3	
			FluxCell_hor(3,:,:,thread_nr)   = Q(2,:,:)*Q(3,:,:)/Q(1,:,:)
			
			GhostFluxLeft(3,:,:,thread_nr)  = GhostLeft(2,:,:,thread_nr)*GhostLeft(3,:,:,thread_nr)/GhostLeft(1,:,:,thread_nr)
			
			GhostFluxRight(3,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)*GhostRight(3,:,:,thread_nr)/GhostRight(1,:,:,thread_nr)
		
			CALL UpdateLaxFriedrichFlux(Q, dt, dx)
					
		END IF
		
		CALL GetFluxDivergence(RQ, dx)
		
		RQ(2,:,:) = RQ(2,:,:) + coriolisPar*Q(3,:,:)
		RQ(3,:,:) = RQ(3,:,:) - coriolisPar*Q(2,:,:)
		
	END SUBROUTINE GetRHS
	
	SUBROUTINE InitializeWENO5(nr_threads, Nx_max)
	
		INTEGER, INTENT(IN) :: nr_threads, Nx_max
		
		INTEGER :: i, thread_nr
		
		! In the 1-D module, Ny is always assumed to be one.
		! The index in the "thread-dimension" starts with zero, so that
		! it directly coincides with thread numbers
		ALLOCATE(GhostLeft(     1:nr_fields, 1, 3, 0:nr_threads-1))
		ALLOCATE(GhostRight(    1:nr_fields, 1, 3, 0:nr_threads-1))
		ALLOCATE(GhostFluxLeft( 1:nr_fields, 1, 3, 0:nr_threads-1))
		ALLOCATE(GhostFluxRight(1:nr_fields, 1, 3, 0:nr_threads-1))
		! NOTE: In the 1-D module, GhostUp and GhostDown are NOT ALLOCATED
		
		! If there are Nx cells, there are Nx+1 interfaces
		ALLOCATE(FluxInt_hor( 1:nr_fields, 1, 1:Nx_max+1, 0:nr_threads-1))
		ALLOCATE(FluxCell_hor(1:nr_fields, 1, 1:Nx_max,   0:nr_threads-1))
		
		! Now perform first-touch initialization, i.e. every thread initializes its
		! part of the buffers
		
		!$OMP PARALLEL DO SCHEDULE(static) private(thread_nr)
		DO i=0,nr_threads-1
			thread_nr = omp_get_thread_num()
			GhostLeft(     :,:,:,thread_nr) = 0.0
			GhostRight(    :,:,:,thread_nr) = 0.0
			GhostFluxLeft( :,:,:,thread_nr) = 0.0
			GhostFluxRight(:,:,:,thread_nr) = 0.0
			FluxInt_hor(   :,:,:,thread_nr) = 0.0
			FluxCell_hor(  :,:,:,thread_nr) = 0.0
		END DO
		!$OMP END PARALLEL DO
		
	END SUBROUTINE InitializeWENO5
	
	! -------------------------------------------
	
	SUBROUTINE UpdateLaxFriedrichFlux(Qcell, dt, dx)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Qcell
		DOUBLE PRECISION,                   INTENT(IN)  :: dt, dx
		
		INTEGER :: Ny, Nx, i, j, thread_nr
		DOUBLE PRECISION :: coeff
		
		thread_nr = omp_get_thread_num()
		
		Ny = SIZE(Qcell, 2)
		Nx = SIZE(Qcell, 3)
		
		coeff = dx/(2.0*dt)
		
		j=1 ! Remains fixed in 1-D case
		
		! Left boundary
		i=0
		FluxInt_hor(:,j,i+1,thread_nr) = DBLE(0.5)*( GhostFluxLeft(:,j,1,thread_nr) + FluxCell_hor(:,j,i+1,thread_nr) ) - coeff*( GhostLeft(:,j,1,thread_nr) + Qcell(:,j,i+1) )
		
		DO i=1,Nx-1
			FluxInt_hor(:,j,i+1,thread_nr) = DBLE(0.5)*( FluxCell_hor(:,j,i,thread_nr) + FluxCell_hor(:,j,i+1,thread_nr) ) - coeff*( Qcell(:,j,i) + Qcell(:,j,i+1) )
		END DO
		
		! Right boundary
		i=Nx
		FluxInt_hor(:,j,i+1,thread_nr) = DBLE(0.5)*( FluxCell_hor(:,j,i,thread_nr) + GhostFluxRight(:,j,1,thread_nr) ) - coeff*( Qcell(:,j,i) + GhostRight(:,j,1,thread_nr) )
		
	END SUBROUTINE UpdateLaxFriedrichFlux
	
	! ---------------------------------------------
	
	SUBROUTINE UpdateHorizontalFlux(Qcell, max_vel)
		
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: Qcell
		DOUBLE PRECISION,                   INTENT(IN) :: max_vel
		
		DOUBLE PRECISION, DIMENSION(nr_fields, 6) :: Qcell_local, Fcell_local
		INTEGER :: Ny, Nx, i, j, thread_nr
		
		thread_nr = omp_get_thread_num()
		
		! Out of the global fields Qcell and FluxQcell, updated interface
		! values of the flux are computed
		Ny = SIZE(Qcell,2) ! Should be =1 in this 1-D module
		Nx = SIZE(Qcell,3)
		
		j=1 ! Remains fixed here
		
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
				! And reconstruct interface value from local stencil
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
		
		! Put the two following subroutines into a CONTAINS block to make it easier for the compiler to inline them.
		! Also, declare them as PURE.
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
				INTEGER j
				
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

	SUBROUTINE GetFluxDivergence(RQ, dx)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: dx
		
		DOUBLE PRECISION :: coeff
		INTEGER          :: i,Nx, thread_nr
		
		Nx        = SIZE(RQ,3)
		coeff     = 1.0/dx
		thread_nr = omp_get_thread_num()
		
		DO i=1,Nx
			RQ(:,1,i) = -coeff*( FluxInt_hor(:,1,i+1,thread_nr) - FluxInt_hor(:,1,i,thread_nr) )
		END DO
		
	END SUBROUTINE GetFluxDivergence
	
	SUBROUTINE FillGhostCells(Q, Nghost, BC, thread_nr)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		INTEGER,						    INTENT(IN)  :: Nghost, thread_nr, BC
		
		INTEGER :: i,Nx
		
		Nx = SIZE(Q, 3)
		
		SELECT CASE (BC)
		
			CASE (1)			
			
				! Fill horizontal ghost cells
				DO i=1,Nghost
					GhostLeft( :,1,i,thread_nr) = Q(:,1,Nx-i+1)
					GhostRight(:,1,i,thread_nr) = Q(:,1,i)
				END DO
									
			CASE (2)
			
				! Fill horizontal ghost cells
				DO i=1,Nghost
					GhostLeft( :,1,i,thread_nr)  = Q(:,1,i)
					GhostRight(:,1,i,thread_nr)  = Q(:,1,Nx-i+1)
				END DO
												
			CASE DEFAULT
			
				WRITE(*, *) ' No implementation available for selected boundary condition'
			
		END SELECT
		
	END SUBROUTINE FillGhostCells

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

        SUBROUTINE FillinGhostcells(Qleft, Qright, Qup, Qdown, Qupleft, Qupright, Qdownleft, Qdownright, Nx, Ny, Nghost)
	
		DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Qleft, Qright, Qup, Qdown, Qupleft, Qupright, Qdownleft, Qdownright
		INTEGER, INTENT(IN) :: Nx, Ny, Nghost

        END SUBROUTINE FillinGhostcells
        
END MODULE WENO5_1D
