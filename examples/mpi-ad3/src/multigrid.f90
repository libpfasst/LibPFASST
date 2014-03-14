module mg
  use iso_c_binding
  use encap
  implicit none
contains
  subroutine relax(y, t, dt, rhs, level, ctx,Nrelax)
    !  Do Nrelax of relaxations
    use probin_mod, only: nu, Lx, spatial_order
    type(pf_encap_t), intent(inout) :: y
    type(pf_encap_t), intent(in)    :: rhs
    real(kind=8),     intent(in)    :: t, dt
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer      :: n, Nx, irb
    real(kind=8) :: dx, sigi
    real(kind=8) :: ym(y%nx), ymm(y%nx), yp(y%nx), ypp(y%nx)

    
    !  Do one sweep of relaxation  
    Nx = y%shape(1)
    dx    = Lx/dble(nvars)

    select case (y%ndim)
    case (1)  
      call relax_1d(y, t, dt, rhs, level, ctx,Nrelax)
    case (2)  
      call relax_2d(y, t, dt, rhs, level, ctx,Nrelax)
    case (3)  
      call relax_3d(y, t, dt, rhs, level, ctx,Nrelax)
    case default
       write(*,*) 'Bad case in multigrid.f90, relax:  ndim=',y%ndim
    end select
  end subroutine relax

  subroutine relax_1d(yptr, t, dt, rhsptr, level, ctx,Nrelax)
    !  Do Nrelax of relaxations in 1d
    use probin_mod, only: nu, Lx, spatial_order
    type(pf_encap_t), intent(inout) :: yptr
    type(pf_encap_t), intent(in)    :: rhs
    real(kind=8),     intent(in)    :: t, dt
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer      :: n, Nx, irb
    real(kind=8) :: dx, sigi
    real(kind=8), allocatable :: ybc(:)
    real(kind=8), pointer :: y(:),rhs(:)

    
    !  Do one sweep of relaxation  
    Nx = yptr%shape(1)
    dx    = Lx/dble(Nx)
    y => yptr%array

    call c_f_pointer(yptr%aptr, y, yptr%shape)
    call c_f_pointer(rhsptr%aptr, rhs, rhsptr%shape)
    allocate(ybc(-spatial_order+1:Nx+spatial_order))

    !  Set values of ybc
    ybc(1:Nx) = y

    select case (spatial_order)
    case (2)  ! Centered diff
       sigi = (dx*dx)/(nu*dt)
       do n = 1,Nrelax
          call set_bc_1d(ybc,Nx,spatial_order)
          do i = 1,Nx-1,2
             ybc(i) = (sigi*rhs(i)+ybc(i-1)+ybc(i+1))/(sigi+TWO)
          end do
          call set_bc_1d(ybc,Nx,spatial_order)
          do i = 2,Nx-1,2
             ybc(i) = (sigi*rhs(i)+ybc(i-1)+ybc(i+1))/(sigi+TWO)
          end do
       end do
!!$    case (4)  ! 4th order Centered diff
!!$       sigi = (12.0_pfdp*dx*dx)/(nu*dt)
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

  subroutine relax_2d(y, t, dt, rhs, level, ctx,Nrelax)
    !  Do Nrelax of relaxations in 2d
    use probin_mod, only: nu, Lx, spatial_order
    type(pf_encap_t), intent(inout) :: y
    type(pf_encap_t), intent(in)    :: rhs
    real(kind=8),     intent(in)    :: t, dt
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer      :: n, Nx,Ny, irb
    real(kind=8) :: dx, sigi
    real(kind=8) :: ym(y%nx), ymm(y%nx), yp(y%nx), ypp(y%nx)

    
    !  Do one sweep of relaxation  
    nvars = size(y%array)
    dx    = Lx/dble(nvars)

    select case (spatial_order)
    case (2)  ! Centered diff
       sigi = (dx*dx)/(nu*dt)
       do n = 1,Nrelax
          do irb = 1,2
             yp = cshift(y%array,-1)
             ym = cshift(y%array,1)
             y%array(irb::2) = (sigi*rhs%array(irb::2) +yp(irb::2)+ym(irb::2))/(sigi+TWO)
          end do
       end do

    case (4)  ! 4th order Centered diff
       sigi = (12.0_pfdp*dx*dx)/(nu*dt)
       do n = 1,Nrelax
          !  Red
          do irb = 1,4
             yp  = cshift(y%array,-1)
             ym  = cshift(y%array,1)
             ypp = cshift(y%array,-2)
             ymm = cshift(y%array,2)
             
             y%array(irb::4) = (sigi*rhs%array(irb::4) &
                  + 16.0_pfdp*(yp(irb::4)+ym(irb::4)) - (ypp(irb::4)+ymm(irb::4)))/(sigi+30.0_pfdp)
          end do
       end do

    case default
       write(*,*) 'Bad case in multigrid.f90, relax:  spatial_order=',spatial_order
    end select
  end subroutine relax_2d
  
  recursive subroutine Vcycle(yptr, t, dt, rhsptr, level, ctx, Nrelax)
    !  Do one Vcycle
    use probin_mod, only: nu, Lx, spatial_order
    use transfer, only: interpolate, restrict
 
    type(pf_encap_t), intent(inout) :: yptr
    type(pf_encap_t), intent(in)    :: rhsptr
    real(kind=8),     intent(in)    :: t, dt
    integer,          intent(in)    :: level
    type(c_ptr),      intent(in)    :: ctx
    integer,          intent(in)    :: Nrelax

    integer          :: i, Nx
    real(kind=8)     :: dx
    integer,    allocatable :: shape_c(:)
    type(pf_encap_t) :: res,resc,cf,c

    Nx = size(yptr%shape(1))
    dx    = Lx/dble(nvars)

    allocate(shape_c(yptr%dim)))

    if (nvars .lt. 8)  then
       !  Do bottom solve
       do i = 1,3
          call relax(yptr, t, dt, rhsptr, level, ctx, Nrelax)
       end do
    else
       kind = 1
       shape_c=yptr%shape/2

       call create(res, level,kind , 0, yptr%shape, c_null_ptr)
       call create(resc, level-1, kind, 0, shape_c, c_null_ptr)
       call create(cf, level, kind, 0, yptr%shape, c_null_ptr)
       call create(c, level-1, kind, 0, shape_c, c_null_ptr)

       !  Do one sweep of relaxation  
       call relax(yptr, t, dt, rhsptr, level, ctx, Nrelax)

       !  Compute residual
       call resid(yptr, t, dt, rhsptr, level, res)
       ! if (level .eq. q) then
       !    print *, 'max resid', maxval(abs(res%array))
       ! end if

       !  Coarsen residual
       call restrict(res, resc, level, ctx, level-1, ctx)       

       !  Call Vcyle recursively
       call setval(c, 0.0d0)
       call Vcycle(c, t, dt, resc, level-1, ctx,Nrelax)

       !  Interpolate correction
       call interpolate(cf, c, level, ctx,level-1 , ctx, do_mg=.true.)

       !  Add correction
       yptr%array = yptr%array+cf%array

       !  Do one sweep of relaxation  
       call relax(y, t, dt, rhs, level, ctx, Nrelax)
       call destroy(res)
       call destroy(resc)
       call destroy(c)
       call destroy(cf)
       deallocate(shape_c)

    end if
  end subroutine Vcycle

  subroutine resid(yptr, t, dt, rhsptr, level, resptr)
    !  Compute the residual for MG
    use probin_mod, only: nu, Lx, spatial_order
    type(pf_encap_t), intent(inout) :: yptr
    type(pf_encap_t), intent(in)    :: rhsprt
    real(kind=8),     intent(in)    :: t, dt
    integer,          intent(in)    :: level
    type(pf_encap_t), intent(inout) :: resptr

    real(kind=8) :: Ndim

    !  Do one sweep of relaxation  
    
    Ndim = yptr%dim

    select case (spatial_order)

    case (2)  ! Centered diff
       call resid_1d(yptr, t, dt, rhsptr, level, resptr)
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
    
    res%array = rhs%array - (y%array - nu*dt*Lap)
  end subroutine resid

  subroutine resid_1d(yptr, t, dt, rhsptr, level, resptr)
    !  Compute the residual for MG
    use probin_mod, only: nu, Lx, spatial_order
    type(pf_encap_t), intent(inout) :: yptr
    type(pf_encap_t), intent(in)    :: rhsptr
    real(kind=8),     intent(in)    :: t, dt
    integer,          intent(in)    :: level
    type(pf_encap_t), intent(inout) :: resptr

    real(kind=8) :: dx,Nx
    real(kind=8) ::  Lap(y%nx)
    real(kind=8), allocatable :: ybc(:)
    real(kind=8), allocatable :: Lap(:)
    integer :: i

    !  Do one sweep of relaxation  
    Nx = yptr%shape(1)
    dx = Lx/dble(y%nx)

    allocate(ybc(-spatial_order+1:Nx+spatial_order))
    allocate(Lap(1:Nx))

    select case (spatial_order)

    case (2)  ! Centered diff
       call set_bc_1d(ubc, Nx, spatial_order)
       Lap = (-TWO*ybc + ybc(0:Nx-1) + ybc(2:Nx+1)/(dx*dx)
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
    
    resptr%array = rhsptr%array - (yptr%array - nu*dt*Lap)
    deallocate(ybc)

  end subroutine resid_1d


  subroutine set_bc_1d(ybc, Nx, Nbc)
    integer, intent(in)      ::  Nx, Nbc
    real(kind=8), , intent(inout) :: ybc(-Nbc+1:Nx+Nbc)


    integer   :: i

    do i = 1,Nbc
       ybc(Nx+i)=ybc(i)
       ybc(1-i)=ybc(Nx+1-i)
    end do
    
  end subroutine set_bc_1d
  subroutine set_bc_2d(ybc, Nx, Nbc)
    integer, intent(in)      ::  Nx, Nbc
    real(kind=8), , intent(inout) :: y(:,:)
    integer   :: i,j
    
  end subroutine set_bc_2d
  subroutine set_bc_3d(ybc, Nx, Nbc)
    integer, intent(in)      ::  Nx, Nbc
    real(kind=8), , intent(inout) :: ybc(:,:,:)
    integer   :: i,j,k


    
  end subroutine set_bc_3d
end module mg
