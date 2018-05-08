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
    
!     factor = (-1./(8*pi*pi*alpha)-8*pi*pi)*dexp(Tfin) + (1./((8*pi*pi+1)*alpha) - (1-8*pi*pi))*dexp(t)
!     factor = -t*exp(Tfin)
    factor = -1._pfdp-1._pfdp/(64._pfdp*alpha*pi**4) + (-8._pfdp*pi*pi-1._pfdp/(8._pfdp*pi*pi*alpha))*Tfin &
            + (8._pfdp*pi*pi+1._pfdp/(8._pfdp*pi*pi*alpha))*t
    
    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [0,1]^ndim
          ycoor = (j-1)*dy
          yd(i,j) = factor*sin(2._pfdp*pi*xcoor)*sin(2._pfdp*pi*ycoor)
!           yd(i,j) = factor
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
    
!    factor = -1./(8*pi*pi*alpha)*dexp(Tfin)+1./((8*pi*pi+1)*alpha)*dexp(t)
!     factor = -t*exp(Tfin)+ exp(t)
    factor = -1._pfdp/(64._pfdp*pi**4*alpha) - 1._pfdp/(8._pfdp*pi*pi*alpha)*Tfin + 1._pfdp/(8._pfdp*pi*pi*alpha)*t
    
    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [-1,1]^ndim
          ycoor = (j-1)*dy
          y(i,j) = factor*sin(2._pfdp*pi*xcoor)*sin(2._pfdp*pi*ycoor)
!           y(i,j) = factor
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
    
!     factor = dexp(Tfin)-dexp(t)
!     factor = exp(Tfin)-exp(t)
    factor = Tfin-t
    
    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx 
          ycoor = (j-1)*dy
          p(i,j) = factor*sin(2._pfdp*pi*xcoor)*sin(2._pfdp*pi*ycoor)
!           p(i,j) = factor
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
    

    call exact_p(u, shape, t)
    u = u / (-alpha)
!     print *, "u", t, maxval(u)
    
  end subroutine exact_u
  
  subroutine p_tilde(p, shape, t, tend)    
    integer,    intent(in)    :: shape(2)
    real(pfdp), intent(inout) :: p(shape(1),shape(2))
    real(pfdp), intent(in)    :: t, tend
    
    integer    :: nx, ny, i, j
    real(pfdp) :: dx, dy, xcoor, ycoor, factor
    
    nx = shape(1)
    ny = shape(2)    
    dx = Lx/dble(nx)
    dy = Ly/dble(ny)
    
!     factor = dexp(Tfin)-dexp(t)+(exp(tend)-exp(Tfin))*exp(-8*pi*pi*(tend-t))
!     factor = -dexp(t)+dexp(tend)
    factor = Tfin-t + (tend-Tfin)*exp(-8.0_pfdp*pi*pi*(tend-t))
    
    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [0,1]^ndim
          ycoor = (j-1)*dy
          p(i,j) = factor*sin(2._pfdp*pi*xcoor)*sin(2._pfdp*pi*ycoor)
!           p(i,j) = factor
       end do
    end do
  end subroutine p_tilde
  
  subroutine lap_p_tilde(p, shape, t, tend)    
    integer,    intent(in)    :: shape(2)
    real(pfdp), intent(inout) :: p(shape(1),shape(2))
    real(pfdp), intent(in)    :: t, tend
    
    integer    :: nx, ny, i, j
    real(pfdp) :: dx, dy, xcoor, ycoor, factor
    
    nx = shape(1)
    ny = shape(2)    
    dx = Lx/dble(nx)
    dy = Ly/dble(ny)
    
    factor = dexp(Tfin)-dexp(t)+(exp(tend)-exp(Tfin))*exp(-8*pi*pi*(tend-t))
    factor = factor * (-8*pi*pi)
!     print *, "exact_p, t = " , t, " factor = ", factor
!     factor = dexp(t)*alpha*(1+8*pi*pi)
    
    do i=1, nx
       do j=1, ny
          xcoor = (i-1)*dx  ! domain is [-1,1]^ndim
          ycoor = (j-1)*dy
          p(i,j) = factor*dsin(2.*pi*xcoor)*dsin(2.*pi*ycoor)
!           print *, "exact_p, t = " , t, " factor = ", factor, " val", p(i,j)
       end do
    end do
  end subroutine lap_p_tilde
  
  subroutine ymyd_exact(y, shape, t)    
    integer,    intent(in)    :: shape(2)
    real(pfdp), intent(inout) :: y(shape(1),shape(2))
    real(pfdp), intent(in)    :: t
    
    integer    :: nx, ny, i, j
    real(pfdp) :: dx, dy, xcoor, ycoor, factor
    real(pfdp) :: yd(shape(1),shape(2))
    
    call exact_y(y, shape, t)
    call y_desired(yd, shape, t)
    
    y = y - yd;
  end subroutine ymyd_exact

end module solutions
