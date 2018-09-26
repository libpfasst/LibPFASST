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
  use pf_mod_fftpackage
  implicit none

  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_imexQ_t) :: ad_sweeper_t
     type(pf_fft_t), pointer :: fft_tool     
     integer ::   nx,ny
     complex(pfdp), allocatable :: lap(:,:)         ! Lapclacian operators
     complex(pfdp), allocatable :: ddx(:,:) ! First derivative operators
     complex(pfdp), allocatable :: ddy(:,:) ! First derivative operators     
     
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
    integer     :: i,j,nx,ny
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
    allocate(this%lap(nx,ny))
    allocate(this%ddx(nx,ny))
    allocate(this%ddy(nx,ny))

    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nx,ny],2)

    call this%fft_tool%make_lap_2d(this%lap)
    call this%fft_tool%make_deriv_2d(this%ddx,1)
    call this%fft_tool%make_deriv_2d(this%ddy,2)        
  end subroutine sweeper_setup

  !>  destroy the sweeper type
  subroutine sweeper_destroy(sweeper)
    class(pf_sweeper_t), intent(inout) :: sweeper
    
    class(ad_sweeper_t), pointer :: this
    this => as_ad_sweeper(sweeper)    

    deallocate(this%lap)
    deallocate(this%ddx)
    deallocate(this%ddy)    
    
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)    

  end subroutine sweeper_destroy

  subroutine initial(y_0)
    type(ndsysarray), intent(inout) :: y_0
    call exact(0.0_pfdp, y_0)
  end subroutine initial
    
  !> Routine to return the exact solution
  subroutine exact(t, y_ex)
    use probin, only: nu, a,b, kfreq
    real(pfdp), intent(in)  :: t
    type(ndsysarray), intent(inout) :: y_ex

    integer    :: ny,nx
    integer    :: i,j, ii, k,nbox
    real(pfdp) :: tol, x,y, omega
    real(pfdp),pointer :: u_ex(:,:),v_ex(:,:)
    
    nx=y_ex%arr_shape(1)
    ny=y_ex%arr_shape(2)

    u_ex=>get_array2d(y_ex,1)
    v_ex=>get_array2d(y_ex,2)
    
    omega = two_pi*kfreq
    do j=1, ny
       y = Ly*dble(j-1-ny/2)/dble(ny) - t*b 
       do i = 1, nx
          x = Lx*dble(i-1-nx/2)/dble(nx) - t*a 
          u_ex(i,j) = sin(omega*x)*sin(omega*y)*exp(-2.0_pfdp*omega*omega*nu*t)
          v_ex(i,j) = cos(omega*x)*cos(omega*y)*exp(-2.0_pfdp*omega*omega*nu*t)             
       end do
    end do
    
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
    type(pf_fft_t),     pointer :: fft    
    integer ::  n

    fft => this%fft_tool
    do n=1,2
       yvec  => get_array2d(y,n)
       fvec => get_array2d(f,n)
       wk => fft%get_wk_ptr_2d()
       
       wk = yvec
       call fft%fftf()          
    
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
       call fft%fftb()
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
    complex(pfdp),   pointer :: wk(:,:)
    type(pf_fft_t),     pointer :: fft        
    integer ::  n
    
    fft => this%fft_tool
    wk => fft%get_wk_ptr_2d()
    
    do n=1,2
       yvec=>get_array2d(y,n)
       fvec=>get_array2d(f,n)
       rhsvec=>get_array2d(rhs,n)                

       if (piece == 2) then
          wk = rhsvec
          call fft%fftf()          
          if (imex_stat .eq. 2) then
             wk = wk / (1.0_pfdp - dtq*nu*this%lap) 
          else  ! fully implicit
             wk = wk / (1.0_pfdp - dtq*(-a * this%ddx -b * this%ddy + nu * this%lap))
          end if
       else
          print *,'Bad piece in f_comp ',piece
          call exit(0)
       end if
       call fft%fftb()
       yvec=real(wk)
       fvec = (yvec - rhsvec) / dtq
    end do
    
  end subroutine f_comp
  subroutine interp2(q_f, q_c, sweeper_f, sweeper_c)
    class(ad_sweeper_t), pointer :: sweeper_f, sweeper_c
    real(pfdp), pointer, intent(inout) :: q_f(:,:), q_c(:,:)

    integer      :: nx_f, nx_c, ny_f, ny_c, xrat,yrat,i,j,ii,jj
    complex(pfdp),         pointer ::  wk_f(:,:),wk_c(:,:)    
    type(pf_fft_t),     pointer :: fft_f,fft_c
    
    nx_f=size(q_f,1)
    ny_f=size(q_f,2)
    nx_c=size(q_c,1)
    ny_c=size(q_c,2)
    xrat  = nx_f/nx_c
    yrat  = ny_f/ny_c
    
    if (xrat == 1 .and. yrat==1) then
       q_f = q_c
       return
    endif

    
    fft_c => sweeper_c%fft_tool
    fft_f => sweeper_f%fft_tool    
    wk_f => fft_f%get_wk_ptr_2d()
    wk_c => fft_c%get_wk_ptr_2d()
       
    wk_c = q_c
    call fft_c%fftf()              
    wk_c = wk_c *4.0_pfdp

    !  fill fine spectral reprentations of y with coarse spectrum
    wk_f = 0.0d0
    do j = 1, ny_c
       if (j <= ny_c/2) then
          jj = j
       else if (j > ny_c/2+1) then
          jj = ny_f - ny_c + j
       else
          cycle
       end if
       
       do i = 1, nx_c
          if (i <= nx_c/2) then
             ii = i
          else if (i > nx_c/2+1) then
             ii = nx_f - nx_c + i
          else
             cycle
          end if
          
          wk_f(ii, jj) = wk_c(i, j)
       end do
    end do

    call fft_f%fftb()
    q_f = real(wk_f)    
    
  end subroutine interp2

  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags


    integer :: nvarF, nvarG, xrat,n
    class(ad_sweeper_t), pointer :: sweeper_f, sweeper_c
    type(ndsysarray), pointer :: ndsysarray_F
    real(pfdp),         pointer :: f(:,:), g(:,:)

    sweeper_c => as_ad_sweeper(levelG%ulevel%sweeper)
    sweeper_f => as_ad_sweeper(levelF%ulevel%sweeper)

    ndsysarray_F => cast_as_ndsysarray(qF)
    
    do n=1,ndsysarray_F%ncomp
       f => get_array2d(qF,n) 
       g => get_array2d(qG,n)
       call interp2(f, g, sweeper_f, sweeper_c)  
    end do
    

  end subroutine interpolate

  subroutine restrict(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t   
    integer, intent(in), optional :: flags

    integer :: xrat, yrat,n
    real(pfdp), pointer :: q_f(:,:), q_c(:,:)

    type(ndsysarray), pointer :: ndsysarray_F
    ndsysarray_F => cast_as_ndsysarray(qF)

    
    do n=1,ndsysarray_F%ncomp
       q_f => get_array2d(qF,n)
       q_c => get_array2d(qG,n)
       
       xrat  = size(q_f,1)/size(q_c,1)
       yrat  = size(q_f,2)/size(q_c,2)
       q_c = q_f(::xrat,::yrat)
    end do
 
  end subroutine restrict

end module feval
