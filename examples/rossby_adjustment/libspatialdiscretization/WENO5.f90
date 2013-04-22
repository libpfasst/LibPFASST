MODULE WENO5_2D
! This module provides the implementation of a two-dimensional, finite difference WENO-5 scheme.
!
! Daniel Ruprecht, 8.2.2012
! ICS Lugano

USE FVMParameters, only : nr_fields, BC, c_s, stabFreq, coriolisPar, grav
USE omp_lib,       only : omp_get_thread_num

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Pressure, PressureLeft, PressureRight, PressureUp, PressureDown

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

	SUBROUTINE GetAdvectionRHS(Uadv, Vadv, Q, RQ, dy, dx, order_advection)
	
		DOUBLE PRECISION,                   INTENT(IN)  :: Uadv, Vadv
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: dy, dx
		INTEGER,                            INTENT(IN)  :: order_advection

		INTEGER :: thread_nr

		thread_nr = omp_get_thread_num()
		
		! Evaluate flux function f(u) = U_adv*u on inner cells and ghost cells
		!FluxCell_hor(  :,:,:,thread_nr) = Uadv*Q		
		!GhostFluxLeft( :,:,:,thread_nr) = Uadv*GhostLeft( :,:,:,thread_nr)
		!GhostFluxRight(:,:,:,thread_nr) = Uadv*GhostRight(:,:,:,thread_nr)
		
		!FluxCell_ver(  :,:,:,thread_nr) = Vadv*Q
		!GhostFluxUp(   :,:,:,thread_nr) = Vadv*GhostUp(  :,:,:,thread_nr)
		!GhostFluxDown( :,:,:,thread_nr) = Vadv*GhostDown(:,:,:,thread_nr)
		
		! Nonlinear advection
		FluxCell_hor(  1,:,:,thread_nr) = 0.5*Q(1,:,:)*Q(1,:,:)
		GhostFluxLeft( 1,:,:,thread_nr) = 0.5*GhostLeft( 1,:,:,thread_nr)*GhostLeft( 1,:,:,thread_nr)
		GhostFluxRight(1,:,:,thread_nr) = 0.5*GhostRight(1,:,:,thread_nr)*GhostRight(1,:,:,thread_nr)
		
		FluxCell_ver( 1,:,:,thread_nr)  = 0.5*Q(1,:,:)*Q(1,:,:)
		GhostFluxUp(  1,:,:,thread_nr)  = 0.5*GhostUp(  1,:,:,thread_nr)*GhostUp(  1,:,:,thread_nr)
		GhostFluxDown(1,:,:,thread_nr)  = 0.5*GhostDown(1,:,:,thread_nr)*GhostDown(1,:,:,thread_nr)
		
		
		! Now update interface values of horizontal flux
		CALL UpdateHorizontalFlux(Q, maxval(abs(Q)))
		CALL UpdateVerticalFlux(Q, maxval(abs(Q)))
		
		!CALL UpdateHorizontalFlux(Q, abs(Uadv))
		!CALL UpdateVerticalFlux(  Q, abs(Vadv))
		
		! Compute flux divergence
		CALL GetFluxDivergence(RQ, dy, dx)
		
	END SUBROUTINE GetAdvectionRHS
	
	!
	! Nonlinear, two-dimensional Euler equations
	!
	!
	SUBROUTINE GetEulerRHS(Q, RQ, dy, dx, order_advection)
	
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
		DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: RQ
		DOUBLE PRECISION,                   INTENT(IN)  :: dy, dx
		INTEGER,                            INTENT(IN)  :: order_advection

		DOUBLE PRECISION :: sound_speed, u_max
		INTEGER :: thread_nr
			
		thread_nr = omp_get_thread_num()
		
		CALL FillGhostCells(Q, 3, BC)
		
		! Update the pressure in the domain as well as pressure ghost-cell values
		CALL EquationOfState(Q, sound_speed)
		
		! For Euler equations in conservative form see for example
		! Jebens et al in MWR 2009.

		IF (minval(abs(Q(1,:,:))) < 1e-10) THEN
			WRITE(*,*) ' Near zero density encountered. Now exiting. Dumping solution into SolutionAtExit.dat'
			OPEN(UNIT=1,FILE='SolutionAtExit.dat')
			WRITE(1,'(F25.20)') Q
			CLOSE(1)
			STOP
		ELSE IF ( ANY(ISNAN(Q)) ) THEN
		    WRITE(*,*) ' NaN encounterd, Simulation is most likely diverging. Now exiting'
			OPEN(UNIT=1,FILE='SolutionAtExit.dat')
			WRITE(1,'(F25.20)') Q
			CLOSE(1)
			STOP
		END IF

		! Fluxes for first compoment: density
		FluxCell_hor(  1,:,:,thread_nr) = Q(2,:,:)
		GhostFluxLeft( 1,:,:,thread_nr) = GhostLeft( 2,:,:,thread_nr)
		GhostFluxRight(1,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)
		
		FluxCell_ver(  1,:,:,thread_nr) = Q(3,:,:)
		GhostFluxUp(   1,:,:,thread_nr) = GhostUp(  3,:,:,thread_nr)
		GhostFluxDown( 1,:,:,thread_nr) = GhostDown(3,:,:,thread_nr)
		
		! Fluxes for second compoment: horizontal momentum
		FluxCell_hor(  2,:,:,thread_nr) = Q(2,:,:)*Q(2,:,:)/Q(1,:,:) + Pressure(:,:,thread_nr)
		GhostFluxLeft( 2,:,:,thread_nr) = GhostLeft(2,:,:,thread_nr)*GhostLeft(2,:,:,thread_nr)/GhostLeft(1,:,:,thread_nr) + PressureLeft(:,:,thread_nr) 
		GhostFluxRight(2,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)*GhostRight(2,:,:,thread_nr)/GhostRight(1,:,:,thread_nr) + PressureRight(:,:,thread_nr)

		FluxCell_ver(  2,:,:,thread_nr) = Q(2,:,:)*Q(3,:,:)/Q(1,:,:)
		GhostFluxUp(   2,:,:,thread_nr) = GhostUp(2,:,:,thread_nr)*GhostUp(3,:,:,thread_nr)/GhostUp(1,:,:,thread_nr)
		GhostFluxDown( 2,:,:,thread_nr) = GhostDown(2,:,:,thread_nr)*GhostDown(3,:,:,thread_nr)/GhostDown(1,:,:,thread_nr)
		
		! Fluxes for first compoment: vertical momentum
		FluxCell_hor(  3,:,:,thread_nr) = Q(2,:,:)*Q(3,:,:)/Q(1,:,:)
		GhostFluxLeft( 3,:,:,thread_nr) = GhostLeft(2,:,:,thread_nr)*GhostLeft(3,:,:,thread_nr)/GhostLeft(1,:,:,thread_nr)
		GhostFluxRight(3,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)*GhostRight(3,:,:,thread_nr)/GhostRight(1,:,:,thread_nr)

		FluxCell_ver(  3,:,:,thread_nr) = Q(3,:,:)*Q(3,:,:)/Q(1,:,:) + Pressure(:,:,thread_nr)
		GhostFluxUp(   3,:,:,thread_nr) = GhostUp(3,:,:,thread_nr)*GhostUp(3,:,:,thread_nr)/GhostUp(1,:,:,thread_nr) + PressureUp(:,:,thread_nr)
		GhostFluxDown( 3,:,:,thread_nr) = GhostDown(3,:,:,thread_nr)*GhostDown(3,:,:,thread_nr)/GhostDown(1,:,:,thread_nr) + PressureDown(:,:,thread_nr)
		
		! Fluxes for first compoment: potential temperature times density
		FluxCell_hor(  4,:,:,thread_nr) = Q(2,:,:)*Q(4,:,:)/Q(1,:,:)
		GhostFluxLeft( 4,:,:,thread_nr) = GhostLeft(2,:,:,thread_nr)*GhostLeft(4,:,:,thread_nr)/GhostLeft(1,:,:,thread_nr)
		GhostFluxRight(4,:,:,thread_nr) = GhostRight(2,:,:,thread_nr)*GhostRight(4,:,:,thread_nr)/GhostRight(1,:,:,thread_nr)

		FluxCell_ver(  4,:,:,thread_nr) = Q(3,:,:)*Q(4,:,:)/Q(1,:,:)
		GhostFluxUp(   4,:,:,thread_nr) = GhostUp(3,:,:,thread_nr)*GhostUp(4,:,:,thread_nr)/GhostUp(1,:,:,thread_nr)
		GhostFluxDown( 4,:,:,thread_nr)	= GhostDown(3,:,:,thread_nr)*GhostDown(4,:,:,thread_nr)/GhostDown(1,:,:,thread_nr)													 
		
		u_max = MAX( MAXVAL(ABS(Q(2,:,:))) , MAXVAL(ABS(Q(3,:,:))) )
		
		CALL UpdateHorizontalFlux(Q, u_max + sound_speed)
		CALL UpdateVerticalFlux(  Q, u_max + sound_speed)
				
		CALL GetFluxDivergence(RQ, dy, dx)

		RQ(3,:,:) = RQ(3,:,:) - grav*Q(1,:,:)
		
		CONTAINS 
			
			SUBROUTINE EquationOfState(Q, sound_speed)
			
				DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: Q
				DOUBLE PRECISION,                   INTENT(OUT) :: sound_speed
				
				DOUBLE PRECISION, PARAMETER :: spec_heat=1000.0, gas_const = 273.0
				DOUBLE PRECISION, PARAMETER :: kappa=gas_const/spec_heat, isentrop_exp = 1.0/(1.0-kappa), rho_ref = 1.0, theta_ref = 300.0, p_ground=1000.0,&
					coeff = (1.0/p_ground)*(gas_const*rho_ref*theta_ref/p_ground**kappa)**isentrop_exp, &
					sound_speed_coeff = sqrt(kappa*p_ground/rho_ref)
				
				INTEGER :: thread_nr

			
				thread_nr = omp_get_thread_num()
				
				! Evaluate pressure P from solution components rho and theta				
				Pressure(     :,:,thread_nr) = coeff*(Q(         4,:,:))**isentrop_exp
				PressureLeft( :,:,thread_nr) = coeff*(GhostLeft( 4,:,:,thread_nr))**isentrop_exp
				PressureRight(:,:,thread_nr) = coeff*(GhostRight(4,:,:,thread_nr))**isentrop_exp
				PressureUp(   :,:,thread_nr) = coeff*(GhostUp(   4,:,:,thread_nr))**isentrop_exp
				PressureDown( :,:,thread_nr) = coeff*(GhostDown( 4,:,:,thread_nr))**isentrop_exp
	
				sound_speed = SQRT(MAXVAL( Pressure(:,:,thread_nr)/Q(1,:,:) ))*sound_speed_coeff
				
			END SUBROUTINE EquationOfState
		
	END SUBROUTINE GetEulerRHS
	
	SUBROUTINE InitializeWENO5(nr_threads, Ny_max, Nx_max)
	
		INTEGER, INTENT(IN) :: nr_threads, Ny_max, Nx_max
		
		INTEGER :: i, thread_nr
		
		! In the 1-D module, Ny is always assumed to be one.
		! The index in the "thread-dimension" starts with zero, so that
		! it directly coincides with thread numbers
		ALLOCATE(GhostLeft(     1:nr_fields, Ny_max, 3, 0:nr_threads-1))
		ALLOCATE(GhostRight(    1:nr_fields, Ny_max, 3, 0:nr_threads-1))
		ALLOCATE(GhostFluxLeft( 1:nr_fields, Ny_max, 3, 0:nr_threads-1))
		ALLOCATE(GhostFluxRight(1:nr_fields, Ny_max, 3, 0:nr_threads-1))

		ALLOCATE(GhostUp(       1:nr_fields, 3, Nx_max, 0:nr_threads-1))
		ALLOCATE(GhostDown(     1:nr_fields, 3, Nx_max, 0:nr_threads-1))
		ALLOCATE(GhostFluxUp(   1:nr_fields, 3, Nx_max, 0:nr_threads-1))
		ALLOCATE(GhostFluxDown( 1:nr_fields, 3, Nx_max, 0:nr_threads-1))
		
		! If there are Nx cells, there are Nx+1 interfaces
		ALLOCATE(FluxInt_hor( 1:nr_fields, Ny_max,   Nx_max+1, 0:nr_threads-1))
		ALLOCATE(FluxCell_hor(1:nr_fields, Ny_max,   Nx_max,   0:nr_threads-1))
		ALLOCATE(FluxInt_ver( 1:nr_fields, Ny_max+1, Nx_max,   0:nr_threads-1))
		ALLOCATE(FluxCell_ver(1:nr_fields, Ny_max,   Nx_max,   0:nr_threads-1))
		
		ALLOCATE(Pressure(Ny_max, Nx_max, 0:nr_threads-1))
		ALLOCATE(PressureLeft(Ny_max, 3, 0:nr_threads-1))
		ALLOCATE(PressureRight(Ny_max, 3, 0:nr_threads-1))
		ALLOCATE(PressureUp(3, Nx_max, 0:nr_threads-1))
		ALLOCATE(PressureDown(3, Nx_max, 0:nr_threads-1))
		
		! Now perform first-touch initialization, i.e. every thread initializes its
		! part of the buffers
		
		!$OMP PARALLEL DO SCHEDULE(static) private(thread_nr)
		DO i=0,nr_threads-1
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
			Pressure(      :,:,thread_nr)   = 0.0
			PressureLeft(  :,:,thread_nr)   = 0.0
			PressureRight( :,:,thread_nr)   = 0.0
			PressureUp(    :,:,thread_nr)   = 0.0
			PressureDown(  :,:,thread_nr)   = 0.0
		END DO
		!$OMP END PARALLEL DO
		
	END SUBROUTINE InitializeWENO5
	
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
	
	END SUBROUTINE UpdateVerticalFlux
		
	SUBROUTINE ReconstructInterfaceValue(Qcell_local, Fcell_local, local_vel, Fint)
	
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
	
	
END MODULE WENO5_2D