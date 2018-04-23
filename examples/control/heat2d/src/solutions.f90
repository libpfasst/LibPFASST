module solutions
  use pf_mod_dtype
  use pf_mod_ndarray_oc
  use probin
  implicit none

contains


  ! Set initial condition.
  subroutine initial(q0, t0, tend)
    type(ndarray_oc), intent(inout) :: q0
    real(pfdp),       intent(in)    :: t0, tend
    real(pfdp), pointer :: q0arr(:,:)
    
    q0%yflatarray = 0.0_pfdp
    q0%pflatarray = 0.0_pfdp 
    
    q0arr => get_array2d_oc(q0, 1)
!     print *, "q0arr", shape(q0arr), t0
    call exact_y(q0arr, shape(q0arr), t0)   
  end subroutine initial

  subroutine y_desired(yd, shape, t)
    integer,    intent(in)    :: shape(2)
    real(pfdp), intent(inout) :: yd(shape(1),shape(2))
    real(pfdp), intent(in)    :: t
    
    integer    :: nx, ny, i, j
    real(pfdp) :: dx, dy, xcoor, ycoor, factor
    
    nx = shape(1)
    ny = shape(2)    
    dx = Lx/dble(nx)
    dy = Ly/dble(ny)
    
!     factor = (ndim*pi*pi/4.+4./(ndim*pi*pi*alpha))*dexp(Tfin) + (1-ndim*pi*pi/4.-4./(4*alpha+ndim*pi*pi*alpha))*dexp(t)
    factor = dexp(t)
!     factor = -dexp(t)*(alpha*(1+8*pi*pi)**2-1)
!     factor = (-1./(8*pi*pi*alpha)-8*pi*pi)*dexp(Tfin) - (-1./((8*pi*pi+1)*alpha) - (8*pi*pi-1))*dexp(t)

!     print *, "ydesired, t = " , t, " factor = ", factor
    
    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [-1,1]^ndim
          ycoor = (j-1)*dy
          yd(i,j) = factor*dsin(2.*pi*xcoor)*dsin(2.*pi*ycoor)
!           print *, "ydesired ", i, j, " value = ", yd(i,j), xcoor, ycoor
       end do
    end do
  end subroutine y_desired
  
  subroutine exact_y(y, shape, t)
    integer,    intent(in)    :: shape(2)
    real(pfdp), intent(inout) :: y(shape(1),shape(2))
    real(pfdp), intent(in)    :: t
    
    integer    :: nx, ny, i, j
    real(pfdp) :: dx, dy, xcoor, ycoor, factor
    
    nx = shape(1)
    ny = shape(2)    
    dx = Lx/dble(nx)
    dy = Ly/dble(ny)
    
!     factor = 4.0/(ndim*pi*pi*alpha)*dexp(Tfin)-4/(4*alpha+ndim*pi*pi*alpha)*dexp(t)
    factor = dexp(t)
!    factor = -1./(8*pi*pi*alpha)*dexp(Tfin)+1./((8*pi*pi+1)*alpha)*dexp(t)
!     print *, "exact_y, t = " , t, " factor = ", factor

    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [-1,1]^ndim
          ycoor = (j-1)*dy
          y(i,j) = factor*dsin(2.*pi*xcoor)*dsin(2.*pi*ycoor)
!           print *, "exact_y", i, j, " values = ", y(i,j)
       end do
    end do
!     print *, "exact_y ", minval(y), maxval(y)
  end subroutine exact_y
  
  subroutine exact_p(p, shape, t)    
    integer,    intent(in)    :: shape(2)
    real(pfdp), intent(inout) :: p(shape(1),shape(2))
    real(pfdp), intent(in)    :: t
    
    integer    :: nx, ny, i, j
    real(pfdp) :: dx, dy, xcoor, ycoor, factor
    
    nx = shape(1)
    ny = shape(2)    
    dx = Lx/dble(nx)
    dy = Ly/dble(ny)
    
    factor = dexp(Tfin)-dexp(t)
!     factor = dexp(t)*alpha*(1+8*pi*pi)
    
    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [-1,1]^ndim
          ycoor = (j-1)*dy
          p(i,j) = factor*dsin(2.*pi*xcoor)*dsin(2.*pi*ycoor)
!           print *, "exact_p, t = " , t, " factor = ", factor, " val", p(i,j)
       end do
    end do
  end subroutine exact_p
  
  subroutine exact_u(u, shape, t)
    integer,    intent(in)    :: shape(2)
    real(pfdp), intent(inout) :: u(shape(1),shape(2))
    real(pfdp), intent(in)    :: t
        
    integer    :: nx, ny, i, j
    real(pfdp) :: dx, dy, xcoor, ycoor, factor
    
    nx = shape(1)
    ny = shape(2)    
    dx = Lx/dble(nx)
    dy = Ly/dble(ny)
    

    factor = dexp(t)*(1+8*pi*pi)

    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [-1,1]^ndim
          ycoor = (j-1)*dy
          u(i,j) = factor*dsin(2.*pi*xcoor)*dsin(2.*pi*ycoor)
!           print *, "exact_y", i, j, " values = ", y(i,j)
       end do
    end do
!     call exact_p(u, shape, t)
!     u = u / alpha
!     print *, "exact_u", t, maxval(u), minval(u(:,:))
    
  end subroutine exact_u

end module solutions
