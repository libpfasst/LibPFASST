!! Useful subroutines to return exact solutions for various problems
!
! This file is part of LIBPFASST.
!
!> Module defining exact solutions for various PDEs
module pf_mod_solutions
  use pf_mod_dtype
  use pf_mod_stop
  implicit none

  interface exact_ad_cos
     module procedure exact_ad_cos_1d
     module procedure exact_ad_cos_1ds
     module procedure exact_ad_cos_1dz
     module procedure exact_ad_cos_2d
     module procedure exact_ad_cos_2dz
     module procedure exact_ad_cos_3d
     module procedure exact_ad_cos_3dz
  end interface exact_ad_cos
  interface exact_ad_exp
     module procedure exact_ad_exp_1d
     module procedure exact_ad_exp_1dz
     module procedure exact_ad_exp_2d
     module procedure exact_ad_exp_2dz
     module procedure exact_ad_exp_3d
     module procedure exact_ad_exp_3dz
  end interface exact_ad_exp
  interface exact_burg_sin
     module procedure exact_burg_sin_1d
     module procedure exact_burg_sin_1dz
     module procedure exact_burg_sin_2d
     module procedure exact_burg_sin_2dz
     module procedure exact_burg_sin_3d
     module procedure exact_burg_sin_3dz
  end interface exact_burg_sin
  interface exact_nls
     module procedure exact_nls_1dz
     module procedure exact_nls_2dz
     module procedure exact_nls_3dz
  end interface exact_nls
  interface exact_nls_sg
     module procedure exact_nls_sg_1dz
     module procedure exact_nls_sg_2dz
     module procedure exact_nls_sg_3dz
     module procedure exact_nls_sg_1d
     module procedure exact_nls_sg_2d
     module procedure exact_nls_sg_3d
  end interface exact_nls_sg
  interface exact_nls_pert
     module procedure exact_nls_pert_1dz
     module procedure exact_nls_pert_2dz
     module procedure exact_nls_pert_3dz
     module procedure exact_nls_pert_1d
     module procedure exact_nls_pert_2d
     module procedure exact_nls_pert_3d
  end interface exact_nls_pert
  interface exact_kdv
     module procedure exact_kdv_1d
     module procedure exact_kdv_1dz
!!$     module procedure exact_kdv_2d
!!$     module procedure exact_kdv_3d
!!$     module procedure exact_kdv_2dz
!!$     module procedure exact_kdv_3dz
  end interface exact_kdv

  
contains
  !> Routine to return the exact solution for advection diffusion
  
  !> Routine to return the exact solution for advection diffusion
  function ad_cos_ex(t, x,nu,v,kfreq,Lx) result(u)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(in)  :: x
    real(pfdp), intent(in)  :: v    
    real(pfdp), intent(in)  :: nu
    real(pfdp), intent(in)  :: kfreq
    real(pfdp), intent(in)  :: Lx        
    real(pfdp)  :: u

    real(pfdp) ::  omega
    omega = kfreq*two_pi/Lx
    u=cos(omega*(x-t*v))*exp(-omega*omega*nu*t)
  end function ad_cos_ex
  
  subroutine exact_ad_cos_1d(t, uex,nu,v,kfreq,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: nu,v,kfreq(1),Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = ad_cos_ex(t, x,nu,v,kfreq(1),Lx)
    end do

  end subroutine exact_ad_cos_1d
  subroutine exact_ad_cos_1ds(t, uex,nu,v,kfreq,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: nu,v,kfreq,Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = ad_cos_ex(t, x,nu,v,kfreq,Lx)
    end do

  end subroutine exact_ad_cos_1ds
  
  subroutine exact_ad_cos_1dz(t, uex,nu,v,kfreq,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: nu,v,kfreq(1),Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = ad_cos_ex(t, x,nu,v,kfreq(1),Lx)       
    end do

  end subroutine exact_ad_cos_1dz
  
  subroutine exact_ad_cos_2d(t, uex,nu,v,kfreq,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: nu,v(2),kfreq(2),Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,uy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    

    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       uy=ad_cos_ex(t, y,nu,v(2),kfreq(2),Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = ad_cos_ex(t, x,nu,v(1),kfreq(1),Lx(1))*uy
       end do
    end do
    
  end subroutine exact_ad_cos_2d
  
  subroutine exact_ad_cos_2dz(t, uex,nu,v,kfreq,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: nu,v(2),kfreq(2),Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,uy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       uy=ad_cos_ex(t, y,nu,v(2),kfreq(2),Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = ad_cos_ex(t, x,nu,v(1),kfreq(1),Lx(1))*uy
       end do
    end do
       
  end subroutine exact_ad_cos_2dz
  subroutine exact_ad_cos_3d(t, uex,nu,v,kfreq,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: nu,v(3),kfreq(3),Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,uy,uz
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    

    do k = 1, nz
       z = Lx(3)*REAL(k-1,pfdp)/REAL(nz,pfdp) 
       uz=ad_cos_ex(t, z,nu,v(3),kfreq(3),Lx(3))              
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          uy=ad_cos_ex(t, y,nu,v(2),kfreq(2),Lx(2))              
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = ad_cos_ex(t, x,nu,v(1),kfreq(1),Lx(1))*uy*uz
          end do
       end do
    end do
    
  end subroutine exact_ad_cos_3d
  
  subroutine exact_ad_cos_3dz(t, uex,nu,v,kfreq,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: nu,v(3),kfreq(3),Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,uy,uz
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    

    do k = 1, nz
       z = Lx(3)*REAL(k-1,pfdp)/REAL(nz,pfdp) 
       uz=ad_cos_ex(t, z,nu,v(3),kfreq(3),Lx(3))              
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          uy=ad_cos_ex(t, y,nu,v(2),kfreq(2),Lx(2))              
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = ad_cos_ex(t, x,nu,v(1),kfreq(1),Lx(1))*uy*uz
          end do
       end do
    end do

  end subroutine exact_ad_cos_3dz

  function ad_exp_ex(t, x,nu,v,Lx) result(u)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(in)  :: x
    real(pfdp), intent(in)  :: v    
    real(pfdp), intent(in)  :: nu
    real(pfdp), intent(in)  :: Lx        
    real(pfdp)  :: u

    integer    :: nx, i, ii, k,nbox
        
    real(pfdp) ::  xx,c,t0

    u=0.0_pfdp
    if (nu .gt. 0.0) then
       t0=0.0025_pfdp/nu
       nbox = ceiling(sqrt(4.0_pfdp*nu*(t0+t)*37.0_pfdp))  !  Decide how many periodic images
       do k = -nbox,nbox
          xx = x- 0.5_pfdp*Lx - t*v + REAL(k,pfdp)*Lx
          u = u + sqrt(t0)/sqrt(t0+t)*exp(-xx*xx/(4.0_pfdp*nu*(t0+t)))
       end do
    else
       nbox = ceiling(sqrt(37.0d0))  !  Decide how many periodic images
       do k = -nbox,nbox
          xx = x - 0.5_pfdp*Lx- t*v  + REAL(k,pfdp)*Lx
          u = u + exp(-xx*xx/(4.0_pfdp*0.0025_pfdp))
       end do
       
    end if  

  end function ad_exp_ex
  
  subroutine exact_ad_exp_1d(t, uex,nu,v,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: nu,v,Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = ad_exp_ex(t, x,nu,v,Lx)
    end do

  end subroutine exact_ad_exp_1d
  
  subroutine exact_ad_exp_1dz(t, uex,nu,v,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: nu,v,Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = ad_exp_ex(t, x,nu,v,Lx)       
    end do

  end subroutine exact_ad_exp_1dz
  
  subroutine exact_ad_exp_2d(t, uex,nu,v,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: nu,v(2),Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,uy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    

    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       uy=ad_exp_ex(t, y,nu,v(2),Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = ad_exp_ex(t, x,nu,v(1),Lx(1))*uy
       end do
    end do
    
  end subroutine exact_ad_exp_2d
  
  subroutine exact_ad_exp_2dz(t, uex,nu,v,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: nu,v(2),Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,uy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       uy=ad_exp_ex(t, y,nu,v(2),Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = ad_exp_ex(t, x,nu,v(1),Lx(1))*uy
       end do
    end do
       
  end subroutine exact_ad_exp_2dz
  subroutine exact_ad_exp_3d(t, uex,nu,v,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: nu,v(3),Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,uy,uz
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    

    do k = 1, nz
       z = Lx(3)*REAL(k-1,pfdp)/REAL(nz,pfdp) 
       uz=ad_exp_ex(t, z,nu,v(3),Lx(3))              
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          uy=ad_exp_ex(t, y,nu,v(2),Lx(2))              
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = ad_exp_ex(t, x,nu,v(1),Lx(1))*uy*uz
          end do
       end do
    end do
    
  end subroutine exact_ad_exp_3d
  
  subroutine exact_ad_exp_3dz(t, uex,nu,v,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: nu,v(3),Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,uy,uz
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    

    do k = 1, nz
       z = Lx(3)*REAL(k-1,pfdp)/REAL(nz,pfdp) 
       uz=ad_exp_ex(t, z,nu,v(3),Lx(3))              
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          uy=ad_exp_ex(t, y,nu,v(2),Lx(2))              
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = ad_exp_ex(t, x,nu,v(1),Lx(1))*uy*uz
          end do
       end do
    end do

  end subroutine exact_ad_exp_3dz
  
  !> Routine to return the exact solution for inviscid Burgers based on Platzman
  !>See "A Simple Illustration of a Weak Spectral Cascade", Muraki D,SIAM J Appl Math,2007
  function burg_sin_ex(t, x,nu,Lx) result(u)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(in)  :: x
    real(pfdp), intent(in)  :: nu
    real(pfdp), intent(in)  :: Lx        
    real(pfdp)  :: u

    real(pfdp) :: rn,ts,xs,a,b,c,phi,dphi,o1,o2
    integer :: n,nterms

    nterms = 100
    u=0.0_pfdp
    xs=x
    ts=t

    if (t .gt. 0.0_pfdp) then
       do n =  nterms,1,-1
          rn = REAL(n,pfdp)
          u=u - 2.0_pfdp*bessel_jn(n,nu*rn*ts)/(rn*ts)*sin(rn*xs)
       end do
    else
       u=-sin(xs)    
    end if

  end function burg_sin_ex
  
  
  subroutine exact_burg_sin_1d(t, uex,nu,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: nu    
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = burg_sin_ex(t, x,nu,Lx)
    end do

  end subroutine exact_burg_sin_1d
  
  subroutine exact_burg_sin_1dz(t, uex,nu,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: nu    
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = burg_sin_ex(t, x,nu,Lx)       
    end do

  end subroutine exact_burg_sin_1dz
  
  subroutine exact_burg_sin_2d(t, uex,nu,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: nu    
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,uy,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    L=0.5_pfdp*(Lx(1)+Lx(2))
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = burg_sin_ex(2.0_pfdp*t, x+y,nu,L)
       end do
    end do
    
  end subroutine exact_burg_sin_2d
  
  subroutine exact_burg_sin_2dz(t, uex,nu,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: nu    
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    L=0.5_pfdp*(Lx(1)+Lx(2))
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = burg_sin_ex(2.0_pfdp*t,x+y,nu,L)
       end do
    end do
       
  end subroutine exact_burg_sin_2dz
  subroutine exact_burg_sin_3d(t, uex,nu,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: nu
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    L=(Lx(1)+Lx(2)+Lx(3))/3.0_pfdp
    do k = 1, nz
       z = Lx(3)*REAL(k-1,pfdp)/REAL(nz,pfdp) 
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = burg_sin_ex(3.0_pfdp*t, x+y+z,nu,L)
          end do
       end do
    end do
    
  end subroutine exact_burg_sin_3d
  
  subroutine exact_burg_sin_3dz(t, uex,nu,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: nu
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    L=(Lx(1)+Lx(2)+Lx(3))/3.0_pfdp
    do k = 1, nz
       z = Lx(3)*real(k-1,pfdp)/REAL(nz,pfdp) 
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = burg_sin_ex(3.0_pfdp*t,x+y+z,nu,L)
          end do
       end do
    end do

  end subroutine exact_burg_sin_3dz
  
  
  !> Routine to return the exact solution for the nonlinear Schoedinger Eq  u_t=i2}u|^2u +i u_xx
  !>See Tuncay Aktosun et al 2007 Inverse Problems 23 217
  !>  Note the domain size should be 2*pi
  function nls_ex(t, x,Lx) result(u)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(in)  :: x
    real(pfdp), intent(in)  :: Lx        
    complex(pfdp)  :: u

    real(pfdp) :: a,b,taa,c
    complex(pfdp) :: ae,dn

    a=1.0_pfdp/sqrt(2.0_pfdp)
    b=1.0_pfdp
    taa=1.0_pfdp
    ae=a*exp(zi*t)
    dn=cosh(t) + zi*sinh(t)
    c = 4.0_pfdp
       
    u=exp(zi*c*(x-c*t))*ae*(dn/(cosh(t)-sqrt(2.0_pfdp)/2.0_pfdp*cos(x-2.0_pfdp*c*t))-1.0_pfdp)

  end function nls_ex
  
  subroutine exact_nls_1dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = nls_ex(t, x,Lx)       
    end do

  end subroutine exact_nls_1dz
  
  subroutine exact_nls_2dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    L=0.5_pfdp*(Lx(1)+Lx(2))
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = nls_ex(2.0_pfdp*t,x+y,L)
       end do
    end do
       
  end subroutine exact_nls_2dz
  
  subroutine exact_nls_3dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    L=(Lx(1)+Lx(2)+Lx(3))/3.0_pfdp
    do k = 1, nz
       z = Lx(3)*real(k-1,pfdp)/REAL(nz,pfdp) 
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = nls_ex(3.0_pfdp*t,x+y+z,L)
          end do
       end do
    end do

  end subroutine exact_nls_3dz
  function nls_sg_ex(t, x,Lx) result(u)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(in)  :: x
    real(pfdp), intent(in)  :: Lx        
    real(pfdp)  :: u

    real(pfdp) :: sigma,xx
    complex(pfdp) :: ae,dn
    sigma=0.02_pfdp
    xx=(x-Lx*0.5_pfdp)/(sigma*Lx)   
    u=exp(-(xx*xx*xx*xx))

  end function nls_sg_ex
  subroutine exact_nls_sg_1dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = nls_sg_ex(t, x,Lx)       
    end do

  end subroutine exact_nls_sg_1dz
  
  subroutine exact_nls_sg_2dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       sy=nls_sg_ex(t, y,Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = nls_sg_ex(t, x,Lx(1))*sy
       end do
    end do
       
  end subroutine exact_nls_sg_2dz
  
  subroutine exact_nls_sg_3dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    real(pfdp) :: sz, sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    do k = 1, nz
       z = Lx(3)*real(k-1,pfdp)/REAL(nz,pfdp) 
       sz = nls_sg_ex(t, z,Lx(3))
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          sy = nls_sg_ex(t, y,Lx(2))
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = nls_sg_ex(t, x,Lx(1))*sy*sz
          end do
       end do
    end do

  end subroutine exact_nls_sg_3dz

  subroutine exact_nls_sg_1d(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = nls_sg_ex(t, x,Lx)       
    end do

  end subroutine exact_nls_sg_1d
  
  subroutine exact_nls_sg_2d(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       sy=nls_sg_ex(t, y,Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = nls_sg_ex(t, x,Lx(1))*sy
       end do
    end do
       
  end subroutine exact_nls_sg_2d
  
  subroutine exact_nls_sg_3d(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    real(pfdp) :: sz, sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    do k = 1, nz
       z = Lx(3)*real(k-1,pfdp)/REAL(nz,pfdp) 
       sz = nls_sg_ex(t, z,Lx(3))
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          sy = nls_sg_ex(t, y,Lx(2))
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = nls_sg_ex(t, x,Lx(1))*sy*sz
          end do
       end do
    end do

  end subroutine exact_nls_sg_3d
  function nls_pert_ex(t, x,Lx) result(u)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(in)  :: x
    real(pfdp), intent(in)  :: Lx        
    complex(pfdp)  :: u

    real(pfdp) :: eps
    eps=0.01_pfdp

    !    u=1.0_pfdp+eps*exp(zi*x*0.25_pfdp-2.0_pfdp*two_pi)
    !    u=1.0_pfdp+eps*exp(zi*(x-2.0_pfdp*two_pi)*0.25_pfdp)
    u=1.0_pfdp+eps*exp(zi*(x)*0.25_pfdp)        

  end function nls_pert_ex
  subroutine exact_nls_pert_1dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp)-0.5_pfdp*Lx 
       uex(i) = nls_pert_ex(t, x,Lx)
    end do

  end subroutine exact_nls_pert_1dz
  
  subroutine exact_nls_pert_2dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       sy=nls_pert_ex(t, y,Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = nls_pert_ex(t, x,Lx(1))*sy
       end do
    end do
       
  end subroutine exact_nls_pert_2dz
  
  subroutine exact_nls_pert_3dz(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    real(pfdp) :: sz, sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    do k = 1, nz
       z = Lx(3)*real(k-1,pfdp)/REAL(nz,pfdp) 
       sz = nls_pert_ex(t, z,Lx(3))
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          sy = nls_pert_ex(t, y,Lx(2))
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = nls_pert_ex(t, x,Lx(1))*sy*sz
          end do
       end do
    end do

  end subroutine exact_nls_pert_3dz
  subroutine exact_nls_pert_1d(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = nls_pert_ex(t, x,Lx)       
    end do

  end subroutine exact_nls_pert_1d
  
  subroutine exact_nls_pert_2d(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       sy=nls_pert_ex(t, y,Lx(2))              
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = nls_pert_ex(t, x,Lx(1))*sy
       end do
    end do
       
  end subroutine exact_nls_pert_2d
  
  subroutine exact_nls_pert_3d(t, uex,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    real(pfdp) :: sz, sy
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    do k = 1, nz
       z = Lx(3)*real(k-1,pfdp)/REAL(nz,pfdp) 
       sz = nls_pert_ex(t, z,Lx(3))
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          sy = nls_pert_ex(t, y,Lx(2))
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = nls_pert_ex(t, x,Lx(1))*sy*sz
          end do
       end do
    end do

  end subroutine exact_nls_pert_3d

  !  Soliton exact solution to kdv  moving at speed one (must scale equation properly)
  function kdv_ex(t, x,beta,Lx) result(u)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(in)  :: x
    real(pfdp), intent(in)  :: beta
    real(pfdp), intent(in)  :: Lx        
    real(pfdp)  :: u
    
    real(pfdp) :: a,s
    
    a=beta*(x-0.375_pfdp*Lx-t)
    s = 2.0_pfdp/(exp(a)+exp(-a))
    u=s*s
    
  end function kdv_ex
  
  subroutine exact_kdv_1dz(t, uex,beta,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in)  :: beta
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = kdv_ex(t, x,beta,Lx)       
    end do

  end subroutine exact_kdv_1dz
  subroutine exact_kdv_1d(t, uex,beta,Lx)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(inout) :: uex(:)
    real(pfdp), intent(in)  :: beta
    real(pfdp), intent(in) :: Lx
    
    integer    :: nx, i
    real(pfdp) :: x
    
    nx = SIZE(uex)
    do i = 1, nx
       x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) 
       uex(i) = kdv_ex(t, x,beta,Lx)       
    end do

  end subroutine exact_kdv_1d
  
  subroutine exact_kdv_2dz(t, uex,beta,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:)
    real(pfdp), intent(in)  :: beta
    real(pfdp), intent(in) :: Lx(2)
    
    integer    :: nx,ny, i,j
    real(pfdp) :: x, y,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)
    L=0.5_pfdp*(Lx(1)+Lx(2))
    do j = 1, ny
       y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
       do i = 1, nx
          x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
          uex(i,j) = kdv_ex(2.0_pfdp*t,x,beta,L)
       end do
    end do
       
  end subroutine exact_kdv_2dz
  
  subroutine exact_kdv_3dz(t, uex,beta,Lx)
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: uex(:,:,:)
    real(pfdp), intent(in)  :: beta
    real(pfdp), intent(in) :: Lx(3)
    
    integer    :: nx,ny,nz, i,j,k
    real(pfdp) :: x, y,z,L
    
    nx = SIZE(uex,1)
    ny = SIZE(uex,2)    
    nz = SIZE(uex,3)    
    L=(Lx(1)+Lx(2)+Lx(3))/3.0_pfdp
    do k = 1, nz
       z = Lx(3)*REAL(k-1,pfdp)/REAL(nz,pfdp) 
       do j = 1, ny
          y = Lx(2)*REAL(j-1,pfdp)/REAL(ny,pfdp)
          do i = 1, nx
             x = Lx(1)*REAL(i-1,pfdp)/REAL(nx,pfdp) 
             uex(i,j,k) = kdv_ex(3.0_pfdp*t,x,beta,L)
          end do
       end do
    end do

  end subroutine exact_kdv_3dz

  !  Initial conditions for ML games
  subroutine init_ML_1dz(u,Lx,ic_bar,ic_type)
    complex(pfdp), intent(inout) :: u(:)
    real(pfdp), intent(in) :: Lx
    integer, intent(in) :: ic_type
    real(pfdp), intent(in) :: ic_bar(:)
    
    integer    :: nx, i,j
    real(pfdp) :: x, xrand(8,2),ixrand(2)

    nx = SIZE(u)
    select case (ic_type)
    case (0)
       !  eight random fourier modes
       call random_number(xrand)
       do i = 1, nx
          x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp)
          u(i)=0.0_pfdp
          do j = 1,8
             u(i) = u(i)+ xrand(j,1)*sin(real(j,pfdp)*x)+xrand(j,2)*cos(real(j,pfdp)*x)       
          end do
       end do
    case (1)
       ! weird exponential thing
       call random_number(ixrand)
       ixrand=real(ceiling(ixrand*7))
       do i = 1, nx
          x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp)
          u(i)=sin(ixrand(1)*x)*exp(cos(ixrand(2)*x))
       end do
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',ic_type)
    end select

  end subroutine init_ML_1dz
  
  
end module pf_mod_solutions
