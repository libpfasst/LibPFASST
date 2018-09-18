!
! This file is part of LIBPFASST.
!
!
! RHS routines for 2-D advection/diffusion example.
!     u_t + v*u_x = nu*(u_xx + u_yy)
module feval
  use pf_mod_dtype
  use pf_mod_ndsysarray
  use pf_mod_imexQ
  implicit none

  include 'fftw3.f03'
  real(pfdp), parameter ::  Lx     = 1.0_pfdp    ! domain size
  real(pfdp), parameter ::  Ly     = 1.0_pfdp    ! domain size  

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_imexQ_t) :: ad_sweeper_t
     type(c_ptr) :: ffft, ifft
     integer ::  ierror, nx,ny
     complex(pfdp), pointer :: wk(:,:)              ! work space     
     complex(pfdp), allocatable :: lap(:,:)         ! Lapclacian operators
     complex(pfdp), allocatable :: ddy(:,:), ddx(:,:) ! First derivative operators
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves 

  end type ad_sweeper_t

contains

  function as_ad_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(ad_sweeper_t), pointer :: r
    select type(sweeper)
    type is (ad_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_ad_sweeper

  subroutine sweeper_setup(sweeper, arr_shape)
    use probin, only:  imex_stat 
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: arr_shape(3)

    class(ad_sweeper_t), pointer :: this
    integer     :: i,j,ierror,nx,ny
    type(c_ptr) :: wk
    real(pfdp)  :: kxj,kxi

    this => as_ad_sweeper(sweeper)

    if (imex_stat .eq. 0 ) then
       this%explicit=.TRUE.
       this%implicit=.FALSE.
    elseif (imex_stat .eq. 1 ) then
       this%implicit=.TRUE.
       this%explicit=.FALSE.
    else
       this%implicit=.TRUE.
       this%explicit=.TRUE.
    end if  

    nx = arr_shape(1)
    ny = arr_shape(2)
    this%nx = nx
    this%ny = ny

    ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nx*ny, c_size_t))
    call c_f_pointer(wk, this%wk, [nx,ny])
    
    this%ffft = fftw_plan_dft_2d(nx,ny, &
         this%wk, this%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    this%ifft = fftw_plan_dft_2d(nx,ny, &
         this%wk, this%wk, FFTW_BACKWARD, FFTW_ESTIMATE)


    allocate(this%lap(nx,ny))
    allocate(this%ddx(nx,ny))
    allocate(this%ddy(nx,ny))        
    do j = 1, ny
       if (j <= ny/2+1) then
          kxj = two_pi / Ly * dble(j-1)
       else
          kxj = two_pi / Ly * dble(-ny + j - 1)
       end if
       do i = 1, nx
       if (i <= nx/2+1) then
          kxi = two_pi / Lx * dble(i-1)
       else
          kxi = two_pi / Lx * dble(-nx + i - 1)
       end if

       if (kxi**2+kxj**2 < 1e-13) then
          this%lap(i,j) = 0.0_pfdp
       else
          this%lap(i,j) = -(kxi**2+kxj**2)
       end if
       this%ddx(i,j) = (0.0_pfdp,1.0_pfdp)*kxi
       this%ddy(i,j) = (0.0_pfdp,1.0_pfdp)*kxj       
    end do
    end do
  end subroutine sweeper_setup
    !>  destroy the sweeper type

  subroutine sweeper_destroy(sweeper)
    class(pf_sweeper_t), intent(inout) :: sweeper
    
    class(ad_sweeper_t), pointer :: this
    this => as_ad_sweeper(sweeper)    
    deallocate(this%wk)
    deallocate(this%lap)
    deallocate(this%ddx)
    deallocate(this%ddy)
    call fftw_destroy_plan(this%ffft)
    call fftw_destroy_plan(this%ifft)
    
    
!    call this%imexQ_destroy(lev)

  end subroutine sweeper_destroy

  subroutine initial(y_0)
    type(ndsysarray), intent(inout) :: y_0
    call exact(0.0_pfdp, y_0)
  end subroutine initial
    
  !> Routine to return the exact solution
  subroutine exact(t, y_ex)
    use probin, only: nprob,nu, a,b, t00, kfreq
    real(pfdp), intent(in)  :: t
    type(ndsysarray), intent(inout) :: y_ex

    integer    :: ny,nx
    integer    :: i,j, ii, k,nbox
    real(pfdp) :: tol, x,y, t0,Dx, omega
    real(pfdp),pointer :: u_ex(:,:),v_ex(:,:)
    
    nx=y_ex%arr_shape(1)
    ny=y_ex%arr_shape(2)

    u_ex=>get_array2d(y_ex,1)
    v_ex=>get_array2d(y_ex,2)
    

!    if (nprob .eq. 1) then
       !  Using sin wave initial condition
       omega = 2*pi*kfreq
       do j=1, ny
          y = Ly*dble(j-1-ny/2)/dble(ny) - t*b 
          do i = 1, nx
             x = Lx*dble(i-1-nx/2)/dble(nx) - t*a 
             u_ex(i,j) = sin(omega*x)*sin(omega*y)*exp(-2.0_pfdp*omega*omega*nu*t)
             v_ex(i,j) = cos(omega*x)*cos(omega*y)*exp(-2.0_pfdp*omega*omega*nu*t)             
          end do
       end do
       
!!$    else  !  Use periodic image of Gaussians
!!$       yex=0
!!$       if (nu .gt. 0) then
!!$          nbox = ceiling(sqrt(4.0*nu*(t00+t)*37.0d0/(Lx*Lx)))  !  Decide how many periodic images
!!$          do k = -nbox,nbox
!!$             do i = 1, nx
!!$                x = (i-1)*Dx-0.5d0 - t*v + dble(k)*Lx
!!$                yex(i,:) = yex(i,:) + dsqrt(t00)/dsqrt(t00+t)*dexp(-x*x/(4.0*nu*(t00+t)))
!!$             end do
!!$          end do
!!$       else
!!$          nbox = ceiling(sqrt(37.0d0/(Lx*Lx)))  !  Decide how many periodic images
!!$          do k = -nbox,nbox
!!$             do i = 1, nx
!!$                x = i*Dx-0.5d0 - t*v + dble(k)*Lx
!!$                yex(i) = yex(i) + dexp(-x*x)
!!$             end do
!!$          end do
!!$          
!!$       end if  ! nbox
!       print *,yex
!    end if
  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    use probin, only:  imex_stat ,nu, a, b
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:,:), fvec(:,:)
    complex(pfdp), pointer :: wk(:,:)
    integer ::  ierror,nx,ny,n


    do n=1,2
    yvec  => get_array2d(y,n)
    fvec => get_array2d(f,n)
    wk => this%wk
    wk = yvec
    nx=size(yvec,1)
    ny=size(yvec,2)

    call fftw_execute_dft(this%ffft, wk, wk)
    
    select case (piece)
    case (1)  ! Explicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          wk = (-a * this%ddx -b * this%ddy + nu * this%lap) *  wk
       case (1)  ! Fully Implicit
          print *,'Should not be in this case in feval'
       case (2)  ! IMEX
          wk = (-a * this%ddx -b * this%ddy) *  wk
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case (2)  ! Implicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          print *,'Should not be in this case in feval'
       case (1)  ! Fully Implicit
          wk = (-a * this%ddx -b * this%ddy + nu * this%lap) *  wk
       case (2)  ! IMEX
          wk = (nu * this%lap) *  wk
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select
    wk =  wk / dble(nx*ny)
    call fftw_execute_dft(this%ifft, wk, wk)
    fvec = real(wk)
    end do
  end subroutine f_eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  imex_stat ,nu,a,b
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),      pointer :: yvec(:,:), rhsvec(:,:), fvec(:,:)
    complex(pfdp), pointer :: wk(:,:)
    integer ::  ierror,nx,ny,n
    
    wk => this%wk
    do n=1,2
       yvec=>get_array2d(y,n)
       fvec=>get_array2d(f,n)
       rhsvec=>get_array2d(rhs,n)                
           
       wk = rhsvec
       nx = size(yvec,1)
       ny = size(yvec,2)

       if (piece == 2) then
          
          call fftw_execute_dft(this%ffft, wk, wk)
!          if (ierror .ne. 0) &
!               stop "error  calling cfft1f in f_comp"
          if (imex_stat .eq. 2) then
             wk = wk / (1.0_pfdp - dtq*nu*this%lap) 
          else  ! fully implicit
             wk = wk / (1.0_pfdp - dtq*(-a * this%ddx -b * this%ddy + nu * this%lap))               
          end if
       else
          print *,'Bad piece in f_comp ',piece
          call exit(0)
       end if
       wk = wk /dble(nx*ny)
       call fftw_execute_dft(this%ifft, wk, wk)
       yvec = real(wk)
       fvec = (yvec - rhsvec) / dtq
    end do
    
  end subroutine f_comp
  subroutine interp2(qF, qG, adF, adG)
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),  intent(inout) :: qF(:,:), qG(:,:)

    complex(pfdp), pointer :: wkF(:,:), wkG(:,:)
    integer      :: nxF, nxG, nyF, nyG, xrat,yrat,i,j,ii,jj
    nxF=size(qF,1)
    nyF=size(qF,2)
    nxG=size(qG,1)
    nyG=size(qG,2)
    xrat  = nxF/nxG
    yrat  = nyF/nyG
    

    if (xrat == 1 .and. yrat==1) then
       qF = qG
       return
    endif

    
    wkF => adF%wk
    wkG => adG%wk
       
    wkG = qG
    call fftw_execute_dft(adG%ffft, wkG, wkG)
    wkG = wkG / (nxG*nyG)
       
    wkF = 0.0d0
    do j = 1, nyG
       if (j <= nyG/2) then
          jj = j
       else if (j > nyG/2+1) then
          jj = nyF - nyG + j
       else
          cycle
       end if
       
       do i = 1, nxG
          if (i <= nxG/2) then
             ii = i
          else if (i > nxG/2+1) then
             ii = nxF - nxG + i
          else
             cycle
          end if
          
          wkF(ii, jj) = wkG(i, j)
       end do
    end do

    call fftw_execute_dft(adF%ifft, wkF, wkF)
       
    qF = real(wkF)    
  end subroutine interp2

  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags


    integer :: nvarF, nvarG, xrat,n
    class(ad_sweeper_t), pointer :: adF, adG
    type(ndsysarray), pointer :: ndsysarray_F
    real(pfdp),         pointer :: f(:,:), g(:,:)

    integer ::  ierror
    adG => as_ad_sweeper(levelG%ulevel%sweeper)
    adF => as_ad_sweeper(levelF%ulevel%sweeper)

    ndsysarray_F => cast_as_ndsysarray(qF)
    
    do n=1,ndsysarray_F%ncomp
       f => get_array2d(qF,n); 
       g => get_array2d(qG,n)
       
       call interp2(f, g, adF, adG)  
    end do
    

  end subroutine interpolate

  subroutine restrict(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t   
    integer, intent(in), optional :: flags

    integer :: NxF, NxG,NyF, NyG, xrat, yrat,n
    real(pfdp), pointer :: f(:,:), g(:,:)

    type(ndsysarray), pointer :: ndsysarray_F
    ndsysarray_F => cast_as_ndsysarray(qF)


    
    do n=1,ndsysarray_F%ncomp
       f => get_array2d(qF,n)
       g => get_array2d(qG,n)
       
       NxF=size(f,1)
       NyF=size(f,2)
       NxG=size(g,1)
       NyG=size(g,2)
       xrat  = NxF/NxG
       yrat  = NyF/NyG
       g = f(::xrat,::yrat)
    end do
 
  end subroutine restrict

end module feval
