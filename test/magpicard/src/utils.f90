module utils
  use pf_mod_dtype
  use pf_mod_zndarray
  use probin


  implicit none

  complex(pfdp), parameter :: &
       z1 = (1.0_pfdp, 0.0_pfdp), &
       zm1 = (-1.0_pfdp, 0.0_pfdp), &
       zmi = (0.0_pfdp, -1.0_pfdp)

 contains
  !> Set initial condition.
  subroutine initial(L)

    class(pf_encap_t), intent(inout) :: L
    class(zndarray), pointer :: L_p
    complex(pfdp),      pointer :: L_array(:,:)
    integer :: Nmat
    
    L_p => cast_as_zndarray(L)
    L_array=>get_array2d(L)
    Nmat = L_p%shape(1) !  Assumes a square matrix

    select case(nprob)
    case (1)
       call init_toda(L_array,Nmat)
    case (2)
       call init_Facke(L_array,Nmat)
    case default
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',nprob)
    end select

!    nullify(L_p)
    nullify(L_array)
  end subroutine initial
  
  subroutine init_toda(L_array,Nmat)
    use probin, only:  toda_periodic

    complex(pfdp),  intent(inout),    pointer :: L_array(:,:)
    integer, intent(in) :: Nmat

    real(pfdp), allocatable :: q(:), p(:)
    real(pfdp) :: alpha_toda
    integer ::  i

    allocate( q(Nmat), p(Nmat))

    L_array = z0
    q = 0.0_pfdp
    p = 0.0_pfdp
    
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
       L_array(i,i) = p(i) ! everyone gets a different initial momentum
    enddo

    do i = 1, Nmat-1
       alpha_toda = exp(-0.5_pfdp * (q(i+1) - q(i)))
       if (i == 1 .and. toda_periodic .eqv. .true.) then
          L_array(i, Nmat) = alpha_toda
          L_array(Nmat, i) = alpha_toda
       endif
       L_array(i, i+1) = alpha_toda
       L_array(i+1, i) = alpha_toda
    enddo

    L_array(:,:) = 0.5_pfdp * L_array(:,:)

    deallocate(q, p)
  end subroutine init_toda

  !>  Initial condition  for Facke
  subroutine init_Facke(L_array,Nmat)
    implicit none
    complex(pfdp),  intent(inout),    pointer :: L_array(:,:)
    integer, intent(in) :: Nmat


    real(pfdp), allocatable :: q(:), x(:)
    real(pfdp) :: dx,ip
    integer ::  i
    
    allocate( q(Nmat), x(Nmat))

    L_array = z0

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
       L_array(:,i) = q*q(i)
    end do
    
    deallocate(q, x)

  end subroutine init_Facke
  
  
  subroutine compute_F_toda(L_array,B_array,Nmat,t,level)
    use probin, only: toda_periodic
    ! RHS for Toda lattice problem
    complex(pfdp), intent(inout),  pointer :: L_array(:,:), B_array(:,:)
    integer,intent(in) :: Nmat
    real(pfdp), intent(in) :: t
    integer, intent(in) :: level
    
    integer :: i


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
  
  
  subroutine compute_Facke(L_array,B_array,Nmat,t,level)
    use probin, only: Znuc,E0,Xmax

    !  RHS for fake Fock matrix example
    complex(pfdp), intent(inout),  pointer :: L_array(:,:), B_array(:,:)
    integer,intent(in) :: Nmat
    real(pfdp), intent(in) :: t
    integer, intent(in) :: level
    
    integer :: i,j,n,m
    real(pfdp) :: xi,xj,xn,xm,cst,dx
    real(pfdp),allocatable :: x(:),x1,x2

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
          x1=abs(B_array(i,j))
          do m = 1, Nmat
             do n = 1, Nmat
                cst = E0*L_array(m,n)*conjg(L_array(n,m))*L_array(i,j)*conjg(L_array(j,i))
                B_array(i,j) = B_array(i,j) + cst*two_electron(x(i),x(j),x(m),x(n))
                cst = E0*conjg(L_array(m,n))*L_array(n,m)*L_array(i,j)*conjg(L_array(j,i))                   
                B_array(i,j) =  B_array(i,j) + cst*(-0.5_pfdp*two_electron(x(i),x(n),x(m),x(j))) 
             enddo
          enddo
          x2=abs(B_array(i,j))
!          if (i .eq. j) print *,x1/x2
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
