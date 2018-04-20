!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for 1-D advection/diffusion example.
!     u_t + v*u_x = nu*u_xx
module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imexQ
  use pf_mod_rkstepper

  implicit none

  real(pfdp), parameter :: &
       Lx     = 1.0_pfdp, &    ! domain size
       kfreq  = 8.0_pfdp, &    ! Frequency of initial conditions
       t00    = 0.15_pfdp      ! initial time for exact solution starting from delta function (Gaussian)

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type, extends(pf_user_level_t) :: ad_level_t
   contains

     procedure :: restrict    => restrict
     procedure :: interpolate => interpolate

  end type ad_level_t


  ! Define the derived sweeper type
  type, extends(pf_imexQ_t) :: ad_sweeper_t
!     type(c_ptr) :: ffft, ifft
!     complex(pfdp), pointer :: wk(:)              ! work space
     real(pfdp), allocatable :: wsave(:)           ! work space
     real(pfdp), allocatable :: work(:)             ! work space
     complex(pfdp), allocatable :: workhat(:)      ! work space
     integer ::  lenwrk, lensav, ierror, nx
     complex(pfdp), allocatable :: ddx(:), lap(:) ! operators
   contains

     procedure :: f_eval => ad_sweeper_f_eval
     procedure :: f_comp => ad_sweeper_f_comp

  end type ad_sweeper_t

  ! Define the derived stepper type
  type, extends(pf_ark_t) :: ad_stepper_t
     
     real(pfdp), allocatable :: wsave(:)           ! work space
     real(pfdp), allocatable :: work(:)             ! work space
     complex(pfdp), allocatable :: workhat(:)      ! work space
     integer ::  lenwrk, lensav, ierror, nx
     complex(pfdp), allocatable :: ddx(:), lap(:) ! operators
   contains
       
     procedure :: f_eval => ad_stepper_f_eval
     procedure :: f_comp => ad_stepper_f_comp
     
  end type ad_stepper_t

contains

  ! Helper functions

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

  function as_ad_stepper(stepper) result(r)
    class(pf_stepper_t), intent(inout), target :: stepper
    class(ad_stepper_t), pointer :: r
    select type(stepper)
    type is (ad_stepper_t)
       r => stepper
    class default
       stop
    end select
  end function as_ad_stepper

  subroutine sweeper_setup(sweeper, nx)
    use probin, only: use_LUq, imex_stat 
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: nx

    class(ad_sweeper_t), pointer :: this
    integer     :: i,ierror
    type(c_ptr) :: wk
    real(pfdp)  :: kx

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
              
    this%use_LUq =use_LUq
    this%nx = nx 
    this%lenwrk = 2*nx 
    this%lensav = 2*nx + int(log(real(nx,kind=8))/log(2.0d+00))+4

    ! create complex fft plans
    allocate(this%workhat(nx))   !  complex transform
    allocate(this%work(this%lenwrk))
    allocate(this%wsave(this%lensav))
    call cfft1i ( nx, this%wsave, this%lensav, ierror )
    if (ierror .ne. 0) then
       stop "error  initializing fft"
    end if

    ! create spectral operators
    allocate(this%ddx(nx))
    allocate(this%lap(nx))
    do i = 1, nx
       if (i <= nx/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-nx + i - 1)
       end if

       this%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          this%lap(i) = 0.0_pfdp
       else
          this%lap(i) = -kx**2
       end if
    end do
  end subroutine sweeper_setup
  
  subroutine stepper_setup(stepper, nx)
    use probin, only:  imex_stat 
    class(pf_stepper_t), intent(inout) :: stepper
    integer,             intent(in   ) :: nx

    class(ad_stepper_t), pointer :: this
    integer     :: i,ierror
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    this => as_ad_stepper(stepper)


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
              
    this%nx = nx 
    this%lenwrk = 2*nx 
    this%lensav = 2*nx + int(log(real(nx,kind=8))/log(2.0d+00))+4

    ! create complex fft plans
    allocate(this%workhat(nx))   !  complex transform
    allocate(this%work(this%lenwrk))
    allocate(this%wsave(this%lensav))
    call cfft1i ( nx, this%wsave, this%lensav, ierror )
    if (ierror .ne. 0) then
       stop "error  initializing fft"
    end if

    ! create spectral operators
    allocate(this%ddx(nx))
    allocate(this%lap(nx))
    do i = 1, nx
       if (i <= nx/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-nx + i - 1)
       end if

       this%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          this%lap(i) = 0.0_pfdp
       else
          this%lap(i) = -kx**2
       end if
    end do
  end subroutine stepper_setup
  

  subroutine sweeper_destroy(this, lev)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev

    deallocate(this%work)
    deallocate(this%workhat)
    deallocate(this%wsave)
    deallocate(this%ddx)
    deallocate(this%lap)
!    call fftw_destroy_plan(this%ffft)
!    call fftw_destroy_plan(this%ifft)
    
    call this%imexQ_destroy(lev)

  end subroutine sweeper_destroy
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0
    call exact(0.0_pfdp, q0%flatarray)
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, yex)
    use probin, only: nu, v
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    integer    :: nx, i, ii, nbox
    real(pfdp) :: tol, x

    nx = size(yex)
    yex   = 0

    if (nu .ge. 0) then
       do i = 1, nx
          !          x = Lx*dble(i-nx/2-1)/dble(nx) - t*v
          x = Lx*dble(i)/dble(nx) - t*v
          yex(i) = yex(i) + dsin(2.0_pfdp*pi*x)*dexp(-4.0_pfdp*pi*pi*nu*t)
       end do
    else
       ! decide how many images so that contribution is neglible
       tol  = 1e-16
       !nbox = 1 + ceiling( sqrt( -(4.0*t00)*log((4.0*pi*(t00))**(0.5)*tol) ))
       nbox = 0

       do ii = -nbox, nbox
          do i = 1, nx
             x = Lx*dble(i-nx/2-1)/dble(nx) + ii*Lx - t*v
             yex(i) = yex(i) + 1.0/(4.0*pi*t00)**(0.5)*dexp(-x**2/(4.0*t00))
             x = Lx*dble(i)/dble(nx)  - t*v
             yex(i) =  dsin(two_pi*kfreq*x)
          end do
       end do
    end if
  end subroutine exact


  ! Evaluate the explicit function at y, t.
  subroutine ad_sweeper_f_eval(this, y, t, level_index, f, piece)
    use probin, only:  imex_stat ,nu, v
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:), fvec(:)
    integer ::  ierror

    yvec  => array1(y)
    fvec => array1(f)

    this%workhat = yvec

    call cfft1f (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in f_eval"
    end if

    select case (piece)
    case (1)  ! Explicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          this%workhat = (-v * this%ddx + nu * this%lap) *  this%workhat
       case (1)  ! Fully Implicit
          print *,'Should not be in this case in feval'
       case (2)  ! IMEX
          this%workhat = -v * this%ddx *  this%workhat
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case (2)  ! Implicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          print *,'Should not be in this case in feval'
       case (1)  ! Fully Implicit
          this%workhat = (-v * this%ddx + nu * this%lap) *  this%workhat
       case (2)  ! IMEX
          this%workhat = (nu * this%lap) *  this%workhat
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select

    call cfft1b (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in f_eval"
    end if

    fvec = real(this%workhat)

  end subroutine ad_sweeper_f_eval
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Solve for y and return f2 also.
  subroutine ad_sweeper_f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  imex_stat ,nu,v
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),      pointer :: yvec(:), rhsvec(:), fvec(:)
    integer ::  ierror
    
    if (piece == 2) then
       yvec  => array1(y)
       rhsvec => array1(rhs)
       fvec => array1(f)
       this%workhat = rhsvec
       
       call cfft1f (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
       if (ierror .ne. 0) &
          stop "error  calling cfft1f in f_comp"
       if (imex_stat .eq. 2) then
          this%workhat =  this%workhat/ (1.0_pfdp - nu*dtq*this%lap)
       else  ! fully implicit
          this%workhat =  this%workhat/ (1.0_pfdp - dtq*(-v * this%ddx +nu*this%lap))
       end if
       call cfft1b (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )      
       if (ierror .ne. 0) &
                      stop "error  calling cfft1f in f_comp"
      yvec  = real(this%workhat)
      fvec = (yvec - rhsvec) / dtq
    else
      print *,'Bad piece in f_comp ',piece
      call exit(0)
    end if
  end subroutine ad_sweeper_f_comp
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! Evaluate the explicit function at y, t.
  subroutine ad_stepper_f_eval(this, y, t, level_index, f, piece)
    use probin, only:  imex_stat ,nu, v
    class(ad_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:), fvec(:)
    integer ::  ierror

    yvec  => array1(y)
    fvec => array1(f)

    this%workhat = yvec

    call cfft1f (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in f_eval"
    end if

    select case (piece)
    case (1)  ! Explicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          this%workhat = (-v * this%ddx + nu * this%lap) *  this%workhat
       case (1)  ! Fully Implicit
          print *,'Should not be in this case in feval'
       case (2)  ! IMEX
          this%workhat = -v * this%ddx *  this%workhat
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case (2)  ! Implicit piece
       select case (imex_stat)
       case (0)  ! Fully Explicit        
          print *,'Should not be in this case in feval'
       case (1)  ! Fully Implicit
          this%workhat = (-v * this%ddx + nu * this%lap) *  this%workhat
       case (2)  ! IMEX
          this%workhat = (nu * this%lap) *  this%workhat
       case DEFAULT
          print *,'Bad case for imex_stat in f_eval ', imex_stat
          call exit(0)
       end select
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select

    call cfft1b (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in f_eval"
    end if

    fvec = real(this%workhat)

  end subroutine ad_stepper_f_eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Solve for y and return f2 also.
  subroutine ad_stepper_f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  imex_stat ,nu,v
    class(ad_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),      pointer :: yvec(:), rhsvec(:), fvec(:)
    integer ::  ierror
    
    if (piece == 2) then
       yvec  => array1(y)
       rhsvec => array1(rhs)
       fvec => array1(f)
       this%workhat = rhsvec
       
       call cfft1f (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )
       if (ierror .ne. 0) &
          stop "error  calling cfft1f in f_comp"
       if (imex_stat .eq. 2) then
          this%workhat =  this%workhat/ (1.0_pfdp - nu*dtq*this%lap)
       else  ! fully implicit
          this%workhat =  this%workhat/ (1.0_pfdp - dtq*(-v * this%ddx +nu*this%lap))
       end if
       call cfft1b (this%nx, 1, this%workhat, this%nx, this%wsave, this%lensav, this%work, this%lenwrk, ierror )      
       if (ierror .ne. 0) &
                      stop "error  calling cfft1f in f_comp"
      yvec  = real(this%workhat)
      fvec = (yvec - rhsvec) / dtq
    else
      print *,'Bad piece in f_comp ',piece
      call exit(0)
    end if
  end subroutine ad_stepper_f_comp


  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t    
    integer, intent(in), optional :: flags


    integer :: nvarF, nvarG, xrat
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),         pointer :: f(:), g(:)

    integer ::  ierror
    adG => as_ad_sweeper(levelG%ulevel%sweeper)
    adF => as_ad_sweeper(levelF%ulevel%sweeper)

    f => array1(qF); 
    g => array1(qG)

    nvarF = size(f)
    nvarG = size(g)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       f = g
       return
    endif

    adG%workhat=g
    call cfft1f (adG%nx, 1, adG%workhat, adG%nx, adG%wsave, adG%lensav, adG%work, adG%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1f in interpolate"
    end if

!    call fftw_execute_dft(adG%ffft, wkG, wkG)
!    adG%workhat = adG%workhat / real(nvarG, kind=8)


    adF%workhat = 0.0d0
    adF%workhat(1:nvarG/2) = adG%workhat(1:nvarG/2)
    adF%workhat(nvarF-nvarG/2+2:nvarF) = adG%workhat(nvarG/2+2:nvarG)


!    call fftw_execute_dft(adF%ifft, adF%wk, adF%wk)
    call cfft1b (adF%nx, 1, adF%workhat, adF%nx, adF%wsave, adF%lensav, adF%work, adF%lenwrk, ierror )
    if (ierror .ne. 0) then
       stop "error  calling cfft1b in interpolate"
    end if

    f = real(adF%workhat)

  end subroutine interpolate

  subroutine restrict(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    real(pfdp), pointer :: f(:), g(:)

    integer :: nvarF, nvarG, xrat

    f => array1(qF)
    g => array1(qG)

    nvarF = size(f)
    nvarG = size(g)
    xrat  = nvarF / nvarG

    g = f(::xrat)
  end subroutine restrict

end module feval
