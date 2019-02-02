!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 3-D advection/diffusion example.
!>     u_t + a*u_x+b*u_y+c*u_z = nu*u_xx
module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imexQ
  use pf_mod_fftpackage
  implicit none



  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imexQ_t) :: ad_sweeper_t
     integer ::     nx,ny
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: lap(:,:)         ! Lapclacian operators
     complex(pfdp), allocatable :: ddx(:,:) ! First derivative operators
     complex(pfdp), allocatable :: ddy(:,:) ! First derivative operators
     complex(pfdp), allocatable :: opE(:,:) ! Explicit part operator
     complex(pfdp), allocatable :: opI(:,:) ! Implicit part operator
     
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves 

  end type ad_sweeper_t

contains

  !>  Helper function to return sweeper pointer
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

  !>  Routine to set up sweeper variables and operators
  subroutine sweeper_setup(sweeper, grid_shape)
    use probin, only:  imex_stat,nu, a,b

    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: grid_shape(3)

    class(ad_sweeper_t), pointer :: this
    integer     :: i,nx,ny

    nx=grid_shape(1)
    ny=grid_shape(2)
    this => as_ad_sweeper(sweeper)

    !>  Set variables for explicit and implicit parts
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

    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nx,ny],2)

    allocate(this%lap(nx,ny))
    allocate(this%ddx(nx,ny))
    allocate(this%ddy(nx,ny))
    allocate(this%opE(nx,ny))
    allocate(this%opI(nx,ny))

    call this%fft_tool%make_lap_2d(this%lap)
    call this%fft_tool%make_deriv_2d(this%ddx,1)
    call this%fft_tool%make_deriv_2d(this%ddy,2)
    select case (imex_stat)
       case (0)  ! Fully Explicit        
          this%opE = -a*this%ddx-b*this%ddy + nu * this%lap
          this%opI = 0.0_pfdp          
       case (1)  ! Fully Implicit
          this%opE = 0.0_pfdp
          this%opI = -a*this%ddx-b*this%ddy + nu * this%lap 
       case (2)  ! IMEX
          this%opE = -a*this%ddx-b*this%ddy
          this%opI =  nu * this%lap           
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select

       
  end subroutine sweeper_setup

  !>  destroy the sweeper type
  subroutine sweeper_destroy(sweeper)
    class(pf_sweeper_t), intent(inout) :: sweeper
    
    class(ad_sweeper_t), pointer :: this
    this => as_ad_sweeper(sweeper)    

    deallocate(this%opE)
    deallocate(this%opI)
       deallocate(this%lap)
       deallocate(this%ddx)
       deallocate(this%ddy)    
    !    call this%imexQ_destroy(lev)
    
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)    
  end subroutine sweeper_destroy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    use probin, only:  imex_stat ,nu, a,b
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:,:), fvec(:,:)
    type(pf_fft_t),     pointer :: fft
    complex(pfdp),      pointer :: wk(:,:)

    yvec  => get_array2d(y)
    fvec => get_array2d(f)
    fft => this%fft_tool
    wk => fft%get_wk_ptr_2d()
    
    ! Take the fft of y (stored in fft%workhat)
    wk=yvec

    ! Apply spectral operators
    select case (piece)
    case (1)  ! Explicit piece
       call fft%conv(this%opE)            
    case (2)  ! Implicit piece
       call fft%conv(this%opI)            
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select


    fvec=real(wk)
  end subroutine f_eval

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
    complex(pfdp),      pointer :: wk(:,:)
    type(pf_fft_t),     pointer :: fft

    fft => this%fft_tool
    wk => fft%get_wk_ptr_2d()
    if (piece == 2) then
       yvec  => get_array2d(y)
       rhsvec => get_array2d(rhs)
       fvec => get_array2d(f)

       ! Take the fft of y (stored in fft%workhat)
       wk=rhsvec
       call fft%conv(1.0_pfdp/(1.0_pfdp - dtq*this%opI))
       yvec=real(wk)
       !  The function is easy to derive
       fvec = (yvec - rhsvec) / dtq

    else
       print *,'Bad piece in f_comp ',piece
       call exit(0)
    end if
  end subroutine f_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  These are the transfer functions that must be  provided for the level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF,qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    integer      :: nx_f, nx_c, ny_f, ny_c
    integer      :: irat,jrat,i,j,ii,jj
    class(ad_sweeper_t), pointer :: sweeper_f, sweeper_c
    real(pfdp),         pointer :: yvec_f(:,:), yvec_c(:,:)
    complex(pfdp),         pointer ::  wk_f(:,:),wk_c(:,:)    
    type(pf_fft_t),     pointer :: fft_f,fft_c

    sweeper_c => as_ad_sweeper(levelG%ulevel%sweeper)
    sweeper_f => as_ad_sweeper(levelF%ulevel%sweeper)

    yvec_f => get_array2d(qF) 
    yvec_c => get_array2d(qG)

    nx_f=size(yvec_f,1)
    ny_f=size(yvec_f,2)
    nx_c=size(yvec_c,1)
    ny_c=size(yvec_c,2)
    irat  = nx_f/nx_c
    jrat  = ny_f/ny_c

    if (irat == 1 .and. jrat==1) then
       yvec_f = yvec_c
       return
    endif

    fft_c => sweeper_c%fft_tool
    fft_f => sweeper_f%fft_tool    

    wk_f => fft_f%get_wk_ptr_2d()
    wk_c => fft_c%get_wk_ptr_2d()
    wk_c=yvec_c
    call fft_c%fftf()    

    wk_f = 0.0_pfdp

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
    wk_f=wk_f*4.0_pfdp
    call fft_f%fftb()
    yvec_f=real(wk_f)
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, levelf, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelf  !<  fine level
    class(pf_level_t), intent(inout) :: levelG  !<  coarse level
    class(pf_encap_t), intent(inout) :: qF    !<  fine solution
    class(pf_encap_t), intent(inout) :: qG    !<  coarse solution
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags


    real(pfdp), pointer :: yvec_f(:,:), yvec_c(:,:)  
    integer      :: irat,jrat
    yvec_f => get_array2d(qF)
    yvec_c => get_array2d(qG)

    irat  = size(yvec_f,1)/size(yvec_c,1)
    jrat  = size(yvec_f,2)/size(yvec_c,2)

    if (irat == 1 .and. jrat==1) then
       yvec_c = yvec_f
    else
       yvec_c = yvec_f(::irat,::jrat)       
    endif


  end subroutine restrict


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial(y_0)
    type(ndarray), intent(inout) :: y_0
    real(pfdp), pointer :: yvec(:,:)
    yvec => get_array2d(y_0)    
    call exact(0.0_pfdp, yvec)
  end subroutine initial

  !> Routine to return the exact solution
  subroutine exact(t, yex)
    use probin, only: nu, a,b, kfreq
    real(pfdp), intent(in)  :: t
    
    real(pfdp), intent(out) :: yex(:,:)

    integer    :: nx,ny, i, j
    real(pfdp) :: tol, x,y,  omega

    nx = size(yex,1)
    ny = size(yex,2)

    !  Using sin wave initial condition
    omega = two_pi*kfreq
    do j = 1, ny
       y = dble(j-1-ny/2)/dble(ny) - t*b 
       do i = 1, nx
          x = dble(i-1-nx/2)/dble(nx) - t*a 
          yex(i,j) = sin(omega*x)*sin(omega*y)*exp(-2.0_pfdp*omega*omega*nu*t)
       end do
    end do
  end subroutine exact


end module feval
