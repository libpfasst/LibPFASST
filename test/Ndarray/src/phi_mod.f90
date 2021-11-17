! ============================ Module Description ===========================
! PHI_MOD: Initilizes Phi functions for scaler arguments. Contains:
!   1. phi - initilizes phi functions using contour integral & recursion
!            relation.
! ===========================================================================

module phi_mod

    ! Module Parameters
    use pf_mod_dtype, only:  pfdp
    implicit none
    real(pfdp), parameter :: PI = 3.1415926535897932_pfdp
    complex(pfdp), parameter :: II = cmplx(0.0_pfdp,1.0_pfdp)
     
    ! Cauchy Integral Settings
    real(pfdp),   parameter :: c_tol = 1.0_pfdp         ! Smallest Lambda for contour integral
    real(pfdp),   parameter :: c_R   = 2.0_pfdp*c_tol   ! Contour Radius
    integer,      parameter :: c_M   = 32               ! Number of Points

    ! Interfaces For phi and phiw
    interface phi
        module procedure phi_zscalar, phi_zvector, phi_zmatrix,phi_zmatrix3d, phi_rscalar, phi_rvector, phi_rmatrix
    end interface

    contains

    ! === FUNCTIONS FOR REAL-VALUED L ===========================
    
    subroutine phi_rvector(L, n, P)
        ! Arguments
        real(pfdp),  intent(in)  :: L(:)
        integer,     intent(in) :: n
        real(pfdp),  intent(out) :: P(n+1,size(L))
        ! Local Variables
        integer :: i, m
        complex(pfdp), allocatable :: LC(:)
        complex(pfdp), allocatable :: PC(:, :)
        m = size(L)
        allocate(LC(m))
        allocate(PC(n+1, m))
        do i = 1, m
            LC(i) = CMPLX(L(i), 0, pfdp)
        end do
        call phi_zvector(LC, n, PC)
        P = real(PC)
    end subroutine phi_rvector

    subroutine phi_rscalar(L, n, P)
        ! Arguments
        real(pfdp), intent(in)  :: L
        integer,    intent(in)  :: n
        real(pfdp), intent(out) :: P(n+1)
        ! Local Variables
        real(pfdp) :: PM(n+1, 1)
        real(pfdp) :: LV(1)
        ! Pack Arguments for call to phi_vector
        LV(1) = L
        call phi_rvector(LV, n, PM)
        ! Unpack Result
        P = PM(:,1)
    end subroutine phi_rscalar
    
    subroutine phi_rmatrix(L, n, P)
        ! Arguments
        real(pfdp), intent(in)  :: L(:,:)
        integer,       intent(in)  :: n
        real(pfdp), intent(out) :: P(n + 1, size(L, 1), size(L, 2))
        ! Local Variables
        integer       :: Ln, Lm, d
        real(pfdp), allocatable :: PM(:, :)
        real(pfdp), allocatable :: LV(:)
        ! Pack Arguments for call to phi_vector
        Lm = size(L, 1)
        Ln = size(L, 2)
        d  = Lm * Ln
        allocate(LV(d))
        allocate(PM(n+1, d))
        LV = reshape(L, (/ d /))
        call phi_rvector(LV, n, PM)
        ! Unpack Result
        P = reshape(PM, (/n+1, Ln, Lm/))
    end subroutine phi_rmatrix
    
    ! === FUNCTIONS FOR COMPLEX-VALUED L ===========================

    subroutine phi_zscalar(L, n, P)
        ! Arguments
        complex(pfdp), intent(in)  :: L
        integer,       intent(in)  :: n
        complex(pfdp), intent(out) :: P(n+1)
        ! Local Variables
        complex(pfdp) :: PM(n+1, 1)
        complex(pfdp) :: LV(1)
        ! Pack Arguments for call to phi_vector
        LV(1) = L
        call phi_zvector(LV, n, PM)
        ! Unpack Result
        P = PM(:,1)
    end subroutine phi_zscalar

    subroutine phi_zmatrix(L, n, P)
        ! Arguments
        complex(pfdp), intent(in)  :: L(:,:)
        integer,       intent(in)  :: n
        complex(pfdp), intent(out) :: P(n + 1, size(L, 1), size(L, 2))
        ! Local Variables
        integer       :: Ln, Lm, d
        complex(pfdp), allocatable :: PM(:, :)
        complex(pfdp), allocatable :: LV(:)
        ! Pack Arguments for call to phi_vector
        Lm = size(L, 1)
        Ln = size(L, 2)
        d  = Lm * Ln
        allocate(LV(d))
        allocate(PM(n+1, d))
        LV = reshape(L, (/ d /))
        call phi_zvector(LV, n, PM)
        ! Unpack Result
        P = reshape(PM, (/n+1, Ln, Lm/))
    end subroutine phi_zmatrix

    subroutine phi_zmatrix3d(L, n, P)
        ! Arguments
        complex(pfdp), intent(in)  :: L(:,:,:)
        integer,       intent(in)  :: n
        complex(pfdp), intent(out) :: P(n + 1, size(L, 1), size(L, 2),size(L,3))
        ! Local Variables
        integer       :: Ln, Lm,Lp, d
        complex(pfdp), allocatable :: PM(:, :)
        complex(pfdp), allocatable :: LV(:)
        ! Pack Arguments for call to phi_vector
        Lm = size(L, 1)
        Ln = size(L, 2)
        Lp = size(L, 3)
        d  = Lm * Ln *Lp
        allocate(LV(d))
        allocate(PM(n+1, d))
        LV = reshape(L, (/ d /))
        call phi_zvector(LV, n, PM)
        ! Unpack Result
        P = reshape(PM, (/n+1, Ln, Lm,Lp/))
    end subroutine phi_zmatrix3d


    ! =======================================================================
    ! PHI   Evaluates \phi_i(L) for i=0,...,n and scaler/vector L using
    !       Recursion relation and Cauchy Integral Formula.
    !
    ! Arguments
    !
    !   L   (input) COMPLEX*16 array, dimensions(n)
    !       array of Lambda values cooresponding to PDE linear component
    !
    !   n   (input) INTEGER
    !       highest phi function to initialize.
    !
    !   P   (output) COMPLEX*16, dimensions(n,size(L))
    !       array where on exit P(i,j) = \phi_{i-1}(L(j))
    ! =======================================================================

    subroutine phi_zvector(L, n, P)
        ! Arguments
        complex(pfdp), intent(in)  :: L(:)
        integer,     intent(in)  :: n
        complex(pfdp), intent(out) :: P(n+1,size(L))
        ! Local Variables
        integer     :: i,j,k,nL
        real(pfdp)    :: f
        complex(pfdp) :: z(c_M)
        complex(pfdp) :: Li,Lzi,pp

        nL = size(L)
        P = 0;
        ! Set contour points
        do i=1,c_M
            z(i) = c_R * exp(2.0_pfdp*PI*II*(i-1.0_pfdp)/real(c_M,pfdp))
        enddo
        ! Compute Phi
        do i=1,nL
            Li = L(i)
            if(abs(Li) >= c_tol) then
                ! Direct Formula
                P(1,i) = exp(Li)
                f = 1.0_pfdp
                do j=2,n+1
                    P(j,i) = (P(j-1,i) - 1.0_pfdp/f)/Li
                    f = f*(j-1)
                enddo
            else
                ! Cauchy Integral Formula
                do k=1,c_M
                    Lzi = Li + z(k)
                    pp = exp(Lzi)
                    P(1,i) = P(1,i) + pp/c_M
                    f = 1.0_pfdp;
                    do j=2,n+1
                        pp = (pp - 1.0_pfdp/f)/Lzi
                        P(j,i) = P(j,i) + pp/c_M;
                        f = f*(j-1)
                    enddo
                enddo
                ! remove imaginary roundoff if L(i) is real
                if(aimag(Li) == 0.0_pfdp) then
                  P(:,i) = REAL(P(:,i))
                endif
            end if
        end do
end subroutine phi_zvector
  ! =======================================================================
  ! WEIGHTS   Compute coefficients for finite difference approximation for
  !           the derivatives 1 to m at point z assuming data is known at
  !           points in array x. Based on the program "weights" in
  !           B. Fornberg, "Calculation of weights in finite difference
  !           formulas", SIAM Review 40 (1998), pp. 685-691.
  ! Arguments
  !
  !   z   (input) DOUBLE
  !       location where approximations are to be accurate
  !
  !   x   (input) DOUBLE Array
  !       array containing interpolation points
  !
  !   m   (input) INTEGER
  !       highest derivative for which weights are sought
  !
  !   W   (output) DOUBLE array, dimension(size(x),m+1)
  !       matrix that gives weights at grid locations x for
  !       derivative of order j<=m are found in c(:,j)
  ! =======================================================================

  subroutine weights(z, x, m, W)
    
    ! Arguments
    real(pfdp), intent(in)    :: z
    real(pfdp), intent(inout)    :: x(:)
    integer,    intent(in)    :: m
    real(pfdp), intent(out)   :: W(m+1,m+1)

    ! Variable Declarations
    real(pfdp) :: c1, c2, c3, c4, c5
    integer  :: ii,i,j,k,n,mn

    c1 = 1.0_pfdp
    c4 = x(1) - z
    W  = 0.0_pfdp
    W(1,1) = 1.0_pfdp

    n = SIZE(x)
    do i=2,n
      mn = min(i,m+1)
      c2 = 1.0_pfdp
      c5 = c4
      c4 = x(i) - z
      do j=1,i-1
        c3 = x(i) - x(j)
        c2 = c2*c3;
        if(j == i-1) then
          do k=mn,2,-1
            W(i,k) = c1*(REAL(k-1,pfdp)*W(i-1,k-1) - c5*W(i-1,k))/c2;
          enddo

          W(i,1) = -c1*c5*W(i-1,1)/c2;
        endif
        do k=mn,2,-1
          W(j,k) = (c4*W(j,k) - REAL(k-1,pfdp)*W(j,k-1))/c3;
        enddo
        W(j,1) = c4*W(j,1)/c3;
      enddo
      c1 = c2;
    enddo
    !        end do

  end subroutine weights


end module phi_mod
