module utils
  use pf_mod_dtype
  use probin
  use factory

  implicit none

  complex(pfdp), parameter :: &
       z0 = (0.0_pfdp, 0.0_pfdp), &
       z1 = (1.0_pfdp, 0.0_pfdp), &
       zm1 = (-1.0_pfdp, 0.0_pfdp), &
       zi = (0.0_pfdp, 1.0_pfdp), &
       zmi = (0.0_pfdp, -1.0_pfdp)

 contains
  !> Set initial condition.
  subroutine initial(L)

    class(pf_encap_t), intent(inout) :: L
    class(zndarray), pointer :: L_p

    L_p => cast_as_zndarray(L)

    select case(nprob)
    case (1)
       call init_toda(L_p)
    case (2)
       call init_Facke(L_p)
    case default
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',nprob)
    end select

    nullify(L_p)
  end subroutine initial
  
  subroutine init_toda(L)
    use probin, only:  toda_periodic

    type(zndarray), intent(inout) :: L

    complex(pfdp), allocatable :: u(:,:) !> definition of initial state
    complex(pfdp),      pointer :: L_array(:,:)
    real(pfdp), allocatable :: q(:), p(:)
    real(pfdp) :: tol=1d-15, alpha_toda
    integer :: Nmat, i

    Nmat = L%shape(1)  !  Assume square matrix
    allocate(u(Nmat, Nmat), q(Nmat), p(Nmat))

    u = z0
    q = 0.0_pfdp
    p = 0.0_pfdp
    L_array=>get_array2d(L)
    
    ! See Zanna thesis on "On the Numerical Solution of Isospectral Flows"
    if (toda_periodic .eqv. .true.) then
       q = 0
    else
       do i = 1, Nmat
           q(i) = i - Nmat/2 - 1 ! everyone assumes a unique position on the number line
       enddo
    endif

    do i = 1, Nmat
       if (i <= Nmat/2) then
          p(i) = 4
       elseif (i == Nmat/2+1) then
          p(i) = 0
       else
          p(i) = -4
       endif
       u(i,i) = p(i) ! everyone gets a different initial momentum
    enddo

    do i = 1, Nmat-1
       alpha_toda = exp(-0.5_pfdp * (q(i+1) - q(i)))
       if (i == 1 .and. toda_periodic .eqv. .true.) then
          u(i, Nmat) = alpha_toda
          u(Nmat, i) = alpha_toda
       endif
       u(i, i+1) = alpha_toda
       u(i+1, i) = alpha_toda
    enddo

    u(:,:) = 0.5_pfdp * u(:,:)
    L_array = u
    deallocate(u, q, p)
  end subroutine init_toda

  !>  Initial condition  for Facke
  subroutine init_Facke(L)
    implicit none

    type(zndarray), intent(inout) :: L

    complex(pfdp), allocatable :: u(:,:) !> definition of initial state
    real(pfdp), allocatable :: q(:), x(:)
    real(pfdp) :: dx,ip
    integer :: Nmat, i
    complex(pfdp),      pointer :: L_array(:,:)
    
    L_array=>get_array2d(L)

    Nmat = L%shape(1) !  Assumes a square matrix
    allocate(u(Nmat, Nmat), q(Nmat), x(Nmat))

    u = z0

    dx=2.0_pfdp*Xmax/real(Nmat-1,pfdp)
    do i = 1, Nmat
       x(i)=-Xmax+real(i-1,pfdp)*dx
    end do
    
    !  Give the q's strengths
    do i=1,Nmat
       ip = x(i)*3.0_pfdp
       q(i)=(0.5_pfdp+0.5_pfdp*cos(ip)*cos(ip))
    end do

    !  the matrix form of this is going to be the outer product of q with itself?
    do i=1,Nmat
       u(:,i) = q*q(i)
    end do
    
    L_array = u
    deallocate( u, q, x)

  end subroutine init_Facke
  
  
  subroutine compute_F_toda(L,B,t,level)
    use probin, only: toda_periodic
    ! RHS for Toda lattice problem
    type(zndarray), intent(in) :: L
    type(zndarray), intent(inout) :: B        
    real(pfdp), intent(in) :: t
    integer, intent(in) :: level
    complex(pfdp),      pointer :: L_array(:,:), B_array(:,:)
    
    integer :: i,Nmat
    L_array=>get_array2d(L)
    B_array=>get_array2d(B)    
    Nmat = L%shape(1)  !  Assume square matrix


    do i = 1, Nmat
       B_array(i,i) = 0.0_pfdp
    enddo
    
    do i = 1, Nmat-1
       B_array(i, i+1) = -1.0_pfdp * L_array(i, i+1)
       B_array(i+1, i) = L_array(i, i+1)
    enddo
    
    if (toda_periodic .eqv. .true.) then
       B_array(1, Nmat) = L_array(1, Nmat)
       B_array(Nmat, 1) = -1.0_pfdp * L_array(Nmat, 1)
    endif
    
  end subroutine compute_F_toda
  
  
  subroutine compute_Facke(L,B,t,level)
    use probin, only: Znuc,E0,Xmax

    !  RHS for fake Fock matrix example
    type(zndarray), intent(in) :: L
    type(zndarray), intent(inout) :: B        
    real(pfdp), intent(in) :: t
    integer, intent(in) :: level
    
    integer :: i,j,n,m,Nmat
    real(pfdp) :: xi,xj,xn,xm,cst,dx
    real(pfdp),allocatable :: x(:)
    complex(pfdp),      pointer :: L_array(:,:), B_array(:,:)
    L_array=>get_array2d(L)
    B_array=>get_array2d(B)    

    Nmat = L%shape(1)  !  Assume square matrix
    allocate(x(Nmat))

    dx=2.0_pfdp*Xmax/real(Nmat-1,pfdp)
    do i = 1, Nmat
       x(i)=-Xmax+real(i-1,pfdp)*dx
    end do
    
    !  Here we are going to loop over each particle and do a local
    !  and nonlocal interaction
    
    do i = 1, Nmat
       do j = 1, i
          cst = -Znuc*L_array(i,j)*conjg(L_array(j,i))
          B_array(i,j) = cst*one_electron(x(i),x(j))
          do m = 1, Nmat
             do n = 1, Nmat
                cst = E0*L_array(m,n)*conjg(L_array(n,m))*L_array(i,j)*conjg(L_array(j,i))
                B_array(i,j) = B_array(i,j) + cst*two_electron(x(i),x(j),x(m),x(n))
                cst = E0*conjg(L_array(m,n))*L_array(n,m)*L_array(i,j)*conjg(L_array(j,i))                   
                B_array(i,j) =  B_array(i,j) + cst*(-0.5_pfdp*two_electron(x(i),x(n),x(m),x(j))) 
             enddo
          enddo
          if (j < i) then
             B_array(j,i) = conjg(B_array(i,j))
          end if
       enddo
    end do

    B_array=cmplx(0.0,1.0)*B_array
    deallocate(x)

    
  end subroutine compute_Facke

  function one_electron(xi,xj) result(E)
   real(pfdp), intent(in) :: xi,xj
   real(pfdp) :: E

   real(pfdp) :: dij,dmn,aij,amn,r2

   E=0.0_pfdp
   !  Return some phony single integral
   dij = xi-xj
   aij = 0.5_pfdp*(xi+xj)
   r2 = aij*aij
   if (r2 .gt. 1.0e-13) then
      E= exp(-0.5_pfdp*dij*dij)/r2
   end if
   
   
 end function one_electron
 
  function two_electron(xi,xj,xm,xn) result(E)
    real(pfdp), intent(in) :: xi,xj,xm,xn
    real(pfdp) :: E
    
    real(pfdp) :: dij,dmn,aij,amn,r2
    
    E = 0.0_pfdp
    !  Return some phony double integral
    dij = xi-xj
    dmn = xm-xn
    aij = 0.5_pfdp*(xi+xj)
    amn = 0.5_pfdp*(xm+xn)
    r2 = (aij-amn)*(aij-amn)
    if (r2 .gt. 1.0d-13) then
      E=exp(-0.5_pfdp*dij*dij)*exp(-0.5_pfdp*dmn*dmn)/r2
   end if
   
 end function two_electron

  
end module utils
