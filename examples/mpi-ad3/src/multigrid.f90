

module mg
  use iso_c_binding
  use pf_mod_ndarray
  implicit none
contains
  subroutine relax(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
    !  Do Nrelax of relaxations
    use probin, only:  Lx, spatial_order
    type(ndarray), intent(inout) :: yptr
    type(ndarray),  intent(in)    :: rhsptr
    real(pfdp),     intent(in)    :: t, nudt
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer      :: n, Nx, irb
    real(pfdp) :: dx, sigi

    
    !  Do one sweep of relaxation  
    Nx = yptr%shape(1)
    dx    = Lx/dble(Nx)

    select case (yptr%dim)
    case (1)  
      call relax_1d(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
!!$    case (2)  
!!$      call relax_2d(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
!!$    case (3)  
!!$      call relax_3d(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
    case default
       write(*,*) 'Bad case in multigrid.f90, relax:  ndim=',yptr%dim
    end select
  end subroutine relax

  subroutine relax_1d(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
    !  Do Nrelax of relaxations in 1d
    use probin, only: Lx, spatial_order
    type(ndarray), intent(inout) :: yptr
    type(ndarray), intent(in)    :: rhsptr
    real(pfdp),     intent(in)    :: t, nudt
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer      :: i,n, Nx, irb
    real(pfdp) :: dx, sigi,sig
    real(pfdp), allocatable :: ybc(:)
    real(pfdp), pointer :: y(:),rhs(:)

    
    !  Do one sweep of relaxation  
    Nx = yptr%shape(1)
    dx    = Lx/dble(Nx)

    call c_f_pointer(yptr%aptr, y, yptr%shape)
    call c_f_pointer(rhsptr%aptr, rhs, rhsptr%shape)
    allocate(ybc(1-spatial_order:Nx+spatial_order))


    select case (spatial_order)
    case (2)  ! Centered diff
       sig = nudt/(dx*dx)
       sigi = (dx*dx)/(nudt)
!!$       do n = 1,Nrelax
!!$          call set_bc_1d(ybc,Nx,spatial_order)
!!$          do i = 1,Nx-1,2
!!$             ybc(i) = (rhs(i)+sig*(ybc(i-1)+ybc(i+1)))/(1.0_pfdp + 2.0_pfdp*sig)
!!$          end do
!!$          call set_bc_1d(ybc,Nx,spatial_order)
!!$          do i = 2,Nx,2
!!$             ybc(i) = (rhs(i)+sig*(ybc(i-1)+ybc(i+1)))/(1.0_pfdp + 2.0_pfdp*sig)
!!$          end do
!!$       end do
!!$       y = ybc(1:Nx) 
! Jacobi
!       do n = 1,Nrelax
!          ybc(1:Nx) = y
!          call set_bc_1d(ybc,Nx,spatial_order)
!          y = (rhs+sig*(ybc(2:Nx+1)+ybc(0:Nx-1)))/(1.0_pfdp + 2.0_pfdp*sig)
!       end do
! Gauss-Seidel
       do n = 1,Nrelax
          call fill_bc_1d(y,ybc,Nx,spatial_order)
          do i = 1,Nx
            ybc(i) = (rhs(i)+sig*(ybc(i-1)+ybc(i+1)))/(1.0_pfdp + 2.0_pfdp*sig)
          end do
!          ybc(Nx+1)=ybc(1)
!          i = Nx
!          ybc(i) = (rhs(i)+sig*(ybc(i-1)+ybc(i+1)))/(1.0_pfdp + 2.0_pfdp*sig)

       end do
       y = ybc(1:Nx) 

!!$    case (4)  ! 4th order Centered diff
!!$       sigi = (12.0_pfdp*dx*dx)/(nudt)
!!$       do n = 1,Nrelax
!!$          !  Red
!!$          do irb = 1,4
!!$             yp  = cshift(y%array,-1)
!!$             ym  = cshift(y%array,1)
!!$             ypp = cshift(y%array,-2)
!!$             ymm = cshift(y%array,2)
!!$             
!!$             y%array(irb::4) = (sigi*rhs%array(irb::4) &
!!$                  + 16.0_pfdp*(yp(irb::4)+ym(irb::4)) - (ypp(irb::4)+ymm(irb::4)))/(sigi+30.0_pfdp)
!!$          end do
!!$       end do

    case default
       write(*,*) 'Bad case in multigrid.f90, relax:  spatial_order=',spatial_order
    end select

    deallocate(ybc)
  end subroutine relax_1d

  subroutine relax_2d(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
    !  Do Nrelax of relaxations in 2d
    use probin, only:  Lx, spatial_order
    type(ndarray), intent(inout) :: yptr
    type(ndarray), intent(in)    :: rhsptr
    real(pfdp),     intent(in)    :: t, nudt
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer      :: n, Nx,Ny, irb
    real(pfdp) :: dx, sigi


    
    !  Do one sweep of relaxation  
    Nx = yptr%shape(1)
    dx    = Lx/dble(Nx)

    select case (spatial_order)
    case (2)  ! Centered diff
!!$       sigi = (dx*dx)/(nudt)
!!$       do n = 1,Nrelax
!!$          do irb = 1,2
!!$             yp = cshift(y%array,-1)
!!$             ym = cshift(y%array,1)
!!$             y%array(irb::2) = (sigi*rhs%array(irb::2) +yp(irb::2)+ym(irb::2))/(sigi+TWO)
!!$          end do
!!$       end do
       print *,'no 2d here'
!!$    case (4)  ! 4th order Centered diff
!!$       sigi = (12.0_pfdp*dx*dx)/(nudt)
!!$       do n = 1,Nrelax
!!$          !  Red
!!$          do irb = 1,4
!!$             yp  = cshift(y%array,-1)
!!$             ym  = cshift(y%array,1)
!!$             ypp = cshift(y%array,-2)
!!$             ymm = cshift(y%array,2)
!!$             
!!$             y%array(irb::4) = (sigi*rhs%array(irb::4) &
!!$                  + 16.0_pfdp*(yp(irb::4)+ym(irb::4)) - (ypp(irb::4)+ymm(irb::4)))/(sigi+30.0_pfdp)
!!$          end do
!!$       end do

    case default
       write(*,*) 'Bad case in multigrid.f90, relax:  spatial_order=',spatial_order
    end select
  end subroutine relax_2d
  
  recursive subroutine Vcycle(yptr, t, nudt, rhsptr, level, ctx, Nrelax,maxresid)
    !  Do one Vcycle
    use probin, only:  Lx, spatial_order
!    use transfer, only: interpolate, restrict
 
    type(ndarray), pointer,intent(in) :: yptr
    type(ndarray), pointer, intent(in)    :: rhsptr
    real(pfdp),     intent(in)    :: t, nudt
    real(pfdp),     intent(inout)    :: maxresid
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer          :: i, Nx,ikind
    real(pfdp)     :: dx,maxres0
    integer,    allocatable :: shape_c(:)
!    type(c_ptr) :: p_res,p_resc,p_corrf,p_corr
    type(ndarray), target :: res,resc,corrf,corr

    Nx = yptr%shape(1)
    dx    = Lx/dble(Nx)

    maxres0=maxresid
    if (Nx .lt. 8)  then
       !  Do bottom solve
       do i = 1,4
          call relax(yptr, t, nudt, rhsptr, level, ctx, Nrelax)
       end do
    else
       ikind = 1
       allocate(shape_c(yptr%dim))
       shape_c=yptr%shape/2

       call ndarray_create_simple(res, yptr%shape)
       call ndarray_create_simple(resc,  shape_c)
       call ndarray_create_simple(corrf,  yptr%shape)
       call ndarray_create_simple(corr,  shape_c)
       call resid(yptr, t, nudt, rhsptr, level, res)
!       print *, 'Before all Nx=',Nx, ' max resid= ',maxval(abs(res%flatarray)), Nrelax,  yptr%shape
!           print *, 'res',res%flatarray
!           print *, 'y',yptr%flatarray
       !  Do one sweep of relaxation  
       call relax(yptr, t, nudt, rhsptr, level, ctx, Nrelax)

       !  Compute residual
       call resid(yptr, t, nudt, rhsptr, level, res)
!       print *, 'Nx=',Nx, ' max resid= ',maxval(abs(res%flatarray)), Nrelax
!           print *, 'res',res%flatarray
!           print *, 'y',yptr%flatarray
       if (1 .eq. 0) then
           print *,'=================================================================='
 
           print *,'----------'
           print *, 'rhs',rhsptr%flatarray
           print *,'----------'
           print *, 'y',yptr%flatarray
           print *,'----------'
           print *, 'res',res%flatarray
        end if

       !  Coarsen residual
       call coarsen(res,resc)
!       print *,'----------'
!       print *, 'resc',resc%flatarray

       !  Call Vcyle recursively
       corr%flatarray= 0.0_pfdp
       call Vcycle(corr, t, nudt, resc, level-1, ctx,Nrelax,maxresid)
       

       !  Interpolate correction
       call interp(corrf,corr)

       !  Add correction
       yptr%flatarray = yptr%flatarray+corrf%flatarray
       call resid(yptr, t, nudt, rhsptr, level, res)
       ! if (level .eq. q) then
!           print *, 'After V Nx=',Nx, ' max resid= ',maxval(abs(res%flatarray)), maxval(abs(corr%flatarray)) 
!           print *,  res%flatarray
       ! end if

       !  Do one sweep of relaxation  
       call relax(yptr, t, nudt, rhsptr, level, ctx, Nrelax)
       call resid(yptr, t, nudt, rhsptr, level, res)
       ! if (level .eq. q) then
!           print *, 'After final relax Nx=',Nx, ' max resid= ',maxval(abs(res%flatarray)), Nrelax 
!           print *,  res%flatarray
       ! end if
       maxresid=maxval(abs(res%flatarray))


       deallocate(res%flatarray)
       deallocate(corr%flatarray)
       deallocate(resc%flatarray)
       deallocate(corrf%flatarray)
       deallocate(res%shape)
       deallocate(corr%shape)
       deallocate(resc%shape)
       deallocate(corrf%shape)
       deallocate(shape_c)

    end if
  end subroutine Vcycle

  subroutine resid(yptr, t, nudt, rhsptr, level, resptr)
    !  Compute the residual for MG
    use probin, only:  Lx, spatial_order
    type(ndarray), intent(inout) :: yptr
    type(ndarray), intent(in)    :: rhsptr
    real(pfdp),     intent(in)    :: t, nudt
    integer,          intent(in)    :: level
    type(ndarray), intent(inout) :: resptr

    real(pfdp) :: Ndim

    !  Compute the residual
    
    Ndim = yptr%dim

    select case (yptr%dim)

    case (1)  ! Centered diff
       call resid_1d(yptr, t, nudt, rhsptr, level, resptr)
!!$    case (4)  ! 4th order Centered diff
!!$       yp  = cshift(y%array,-1)
!!$       ym  = cshift(y%array,1)
!!$       ypp = cshift(y%array,-2)
!!$       ymm = cshift(y%array,2)
!!$       Lap = (-30.0_pfdp*y%array + 16.0_pfdp*(yp + ym)-(ypp+ymm))/(12.0_pfdp*dx*dx)
!!$
    case default
       write(*,*) 'Bad case in feval.f90, feval_init:  spatial_order=', spatial_order
    end select
  end subroutine resid

  subroutine coarsen(yptr, ycptr)
    !  coarsen solution y to yc
    use probin, only:  Lx, spatial_order
    type(ndarray), intent(in) :: yptr
    type(ndarray),  intent(in)    :: ycptr
    real(pfdp), pointer :: y(:),yc(:)
    integer      :: n, Nx, irb
    real(pfdp), allocatable :: ybc(:)
    
    call c_f_pointer(yptr%aptr, y, yptr%shape)
    call c_f_pointer(ycptr%aptr, yc, ycptr%shape)

    Nx = yptr%shape(1)


    select case (yptr%dim)
    case (1)  
    !  call relax_1d(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
       allocate(ybc(-spatial_order+1:Nx+spatial_order))
       call fill_bc_1d(y,ybc, Nx, spatial_order)
       yc=0.5_pfdp*ybc(1:Nx-1:2) + 0.25_pfdp*(ybc(0:Nx-2:2)+ybc(2:Nx:2))
       deallocate(ybc)
    case default
       write(*,*) 'Bad case in multigrid.f90, coarsen:  ndim=',yptr%dim
    end select
  end subroutine coarsen

  subroutine interp(yptr, ycptr)
    !  interp  solution  yc to y
    use probin, only:  Lx, spatial_order
    type(ndarray), intent(inout) :: yptr
    type(ndarray),  intent(in)    :: ycptr
    real(pfdp), pointer :: y(:),yc(:)
    integer      :: n, Nx, irb
    real(pfdp), allocatable :: ybc(:)
    
    call c_f_pointer(yptr%aptr, y, yptr%shape)
    call c_f_pointer(ycptr%aptr, yc, ycptr%shape)

    Nx = yptr%shape(1)


    select case (yptr%dim)
    case (1)  
    !  call relax_1d(yptr, t, nudt, rhsptr, level, ctx,Nrelax)
       allocate(ybc(-spatial_order+1:Nx/2+spatial_order))
       call fill_bc_1d(yc,ybc, Nx/2, spatial_order)
       y(1:Nx-1:2) = ybc(1:Nx/2)
!       y(2:Nx:2) = 0.5_pfdp*(ybc(1:Nx/2)+ybc(2:Nx/2+1))
       y(2:Nx:2) = (-ybc(0:Nx/2-1)+9.0_pfdp*(ybc(1:Nx/2)+ybc(2:Nx/2+1)) -ybc(3:Nx/2+2) )/16.0_pfdp
       deallocate(ybc)
    case default
       write(*,*) 'Bad case in multigrid.f90, interp:  ndim=',yptr%dim
    end select
  end subroutine interp


  subroutine resid_1d(yptr, t, nudt, rhsptr, level, resptr)
    !  Compute the residual for MG
    use probin, only: Lx, spatial_order
    type(ndarray), intent(inout) :: yptr
    type(ndarray), intent(in)    :: rhsptr
    real(pfdp),     intent(in)    :: t, nudt
    integer,          intent(in)    :: level
    type(ndarray), intent(inout) :: resptr
    real(pfdp), pointer :: y(:),rhs(:),res(:)
    real(pfdp) :: dx

    real(pfdp), allocatable :: ybc(:)
    real(pfdp), allocatable :: Lap(:)
    integer :: i,Nx

    Nx = yptr%shape(1)
    dx = Lx/dble(Nx)
    call c_f_pointer(yptr%aptr, y, yptr%shape)
    call c_f_pointer(rhsptr%aptr, rhs, rhsptr%shape)
    call c_f_pointer(resptr%aptr, res, resptr%shape)

    allocate(ybc(-spatial_order+1:Nx+spatial_order))
    allocate(Lap(1:Nx))

    select case (yptr%dim)
    case (1)  
       call fill_bc_1d(y,ybc, Nx, spatial_order)
       Lap = (-2.0_pfdp*ybc(1:Nx) + ybc(0:Nx-1) + ybc(2:Nx+1))/(dx*dx)
!!$    case (4)  ! 4th order Centered diff
!!$       yp  = cshift(y%array,-1)
!!$       ym  = cshift(y%array,1)
!!$       ypp = cshift(y%array,-2)
!!$       ymm = cshift(y%array,2)
!!$       Lap = (-30.0_pfdp*y%array + 16.0_pfdp*(yp + ym)-(ypp+ymm))/(12.0_pfdp*dx*dx)
!!$
    case default
       write(*,*) 'Bad case in feval.f90, feval_init:  spatial_order=', spatial_order
    end select
    res = rhs - (y - nudt*Lap)
!!$    print *,'----------------------------------------------------'
!!$    print *,'in resid Nx=',Nx,size(resptr%flatarray )
!!$    print *, 'y',yptr%flatarray 
!!$    print *, 'ybc',ybc
!!$    print *, 'Lap',Lap
!!$    print *, 'rhs',rhsptr%flatarray 
!!$    print *, 'res',resptr%flatarray 

    deallocate(ybc)

  end subroutine resid_1d


  subroutine fill_bc_1d(y,ybc, Nx, Nbc)
    integer, intent(in)      ::  Nx, Nbc
    real(pfdp), intent(inout) :: y(1:Nx)
    real(pfdp), intent(inout) :: ybc(1-Nbc:Nx+Nbc)
    integer   :: i

    ybc(1:Nx)=y
    do i = 1,Nbc
       ybc(Nx+i)=ybc(i)
       ybc(1-i)=ybc(Nx+1-i)
    end do
    
  end subroutine fill_bc_1d
  subroutine set_bc_1d(ybc, Nx, Nbc)
    integer, intent(in)      ::  Nx, Nbc
    real(pfdp), intent(inout) :: ybc(1-Nbc:Nx+Nbc)
    integer   :: i
    do i = 1,Nbc
       ybc(Nx+i)=ybc(i)
       ybc(1-i)=ybc(Nx+1-i)
    end do

!!$    ybc(Nx+1)=ybc(1)
!!$    ybc(Nx+2)=-ybc(Nx)
!!$    ybc(0)=-ybc(1)
!!$    ybc(-1)=-ybc(2)

    
  end subroutine set_bc_1d
  subroutine set_bc_2d(ybc, Nx, Nbc)
    integer, intent(in)      ::  Nx, Nbc
    real(pfdp),  intent(inout) :: ybc(:,:)
    integer   :: i,j
    
  end subroutine set_bc_2d
  subroutine set_bc_3d(ybc, Nx, Nbc)
    integer, intent(in)      ::  Nx, Nbc
    real(pfdp),  intent(inout) :: ybc(:,:,:)
    integer   :: i,j,k


    
  end subroutine set_bc_3d
end module mg
