module feval
  use pf_mod_dtype
  use pf_mod_ndarray_oc
  use pf_mod_restrict
  use pf_mod_imexQ_oc
  use pf_mod_misdcQ_oc
  use pf_mod_fftpackage
  use probin
  use solutions
  implicit none
!   include 'fftw3.f03'

!   real(pfdp), parameter :: pi = 3.141592653589793_pfdp    !defined in probin already
!   real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

type, extends(pf_imexQ_oc_t) :: ad_sweeper_t  !generalize so we can use misdc as well?
!  type, extends(pf_misdcQ_oc_t) :: ad_sweeper_t  !generalize so we can use misdc as well?
!      type(c_ptr) :: ffft, ifft
     type(pf_fft_t), pointer :: fft_tool
!      complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), allocatable :: lap(:) ! operators
     real(pfdp), pointer    :: u(:,:,:)
     real(pfdp), pointer    :: ydesired(:,:,:)
     real(pfdp)             :: alpha
     integer                :: nsteps_per_rank, nproc, myrank
     real(pfdp), allocatable:: newton_f(:), newton_fprime(:), delta_s(:), delta_sbar(:)
   contains

     procedure :: f_eval
     procedure :: f_comp
!     final :: destroy0, destroy1

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

      ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece, flags, idx, step)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    integer, intent(in)              :: flags
    integer, intent(in), optional    :: idx       ! index of quadrature node
    integer, intent(in), optional    :: step   ! time step for sequential version

    real(pfdp),    pointer :: yvec(:), fvec(:), p(:)
    complex(pfdp), pointer :: wk(:)
    integer :: l, loopstart, loopend, thisstep, mystep

    thisstep = 1
    if(present(step)) thisstep = step
    mystep = (thisstep - 1 - this%myrank)/this%nproc + 1
    
!     print *, this%myrank, thisstep, mystep, size(this%u)
    
    select case (piece)
    case (1)  ! Explicit piece
      select case (flags)
      case (0)
        ! first component: y
        fvec => get_array1d_oc(f, 1)
        yvec  => get_array1d_oc(y, 1)
        if (do_imex .eq. 1) then
          fvec = this%u(mystep, idx, :) - (gamma*1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
        else
          fvec = this%u(mystep, idx, :)
       endif
       ! second component: p
       fvec => get_array1d_oc(f, 2)
       if (do_imex .eq. 1) then
         p  => get_array1d_oc(y, 2)
         fvec = (yvec - this%ydesired(mystep, idx,:)) - (gamma*yvec*yvec-1.0_pfdp)*p
       else
         fvec = (yvec - this%ydesired(mystep, idx,:))
       end if
      case (1)
         fvec => get_array1d_oc(f, 1)
         if (do_imex .eq. 1) then
            yvec  => get_array1d_oc(y, 1)
!             print *, yvec
            fvec = this%u(mystep, idx, :) - (gamma*1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
         else
            fvec = this%u(mystep, idx, :)
         endif
       case (2)
         ! evaluate y-y_d
         fvec => get_array1d_oc(f, 2)
         yvec  => get_array1d_oc(y, 1)
         if (do_imex .eq. 1) then
           p  => get_array1d_oc(y, 2)
           fvec = (yvec -this%ydesired(mystep, idx,:)) - (gamma*yvec*yvec-1.0_pfdp)*p
         else
           fvec = (yvec -this%ydesired(mystep, idx,:))
         end if
       case default
         stop "ERROR in f_eval: only 0, 1, 2 allowed as flags"
       end select

    case (2)  ! Implicit piece
       select case (flags)
       case (0)
         loopstart = 1
         loopend = 2
       case (1)
         loopstart = 1
         loopend = 1
       case (2)
         loopstart = 2
         loopend = 2
       case default
         stop "ERROR in f_eval: only 0, 1, 2 allowed as flags"
       end select
       do l = loopstart, loopend
        yvec  => get_array1d_oc(y, l)
        fvec  => get_array1d_oc(f, l)
        call this%fft_tool%conv(yvec,this%lap,fvec)            

!         wk => this%fft_tool%get_wk_ptr_1d() 
!         wk = yvec
! !         call fftw_execute_dft(this%ffft, wk, wk)
! !         wk = nu * this%lap * wk / size(wk)
! !         call fftw_execute_dft(this%ifft, wk, wk)
!         call this%fft_tool%fftf() 
!         wk = nu * this%lap * wk
!         call this%fft_tool%fftb()        
!         fvec = real(wk)
       end do

    case (3) ! third piece (misdc)
    select case (flags)
      case (0)
        ! first component: y
        fvec => get_array1d_oc(f, 1)
        yvec  => get_array1d_oc(y, 1)
        fvec = - (gamma*1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
        ! second component: p
        fvec => get_array1d_oc(f, 2)
        p  => get_array1d_oc(y, 2)
        fvec = - (gamma*yvec*yvec-1.0_pfdp)*p
      case (1)
        fvec => get_array1d_oc(f, 1)
        yvec  => get_array1d_oc(y, 1)
        fvec = - (gamma*1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
      case (2)
        fvec => get_array1d_oc(f, 2)
        yvec  => get_array1d_oc(y, 1)
        p  => get_array1d_oc(y, 2)
        fvec = - (gamma*yvec*yvec-1.0_pfdp)*p
       case default
         stop "ERROR in f_eval: only 0, 1, 2 allowed as flags"
       end select

    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select
  end subroutine f_eval



  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f, piece, flags)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    integer,             intent(in   ) :: flags

    real(pfdp),      pointer :: s(:), rhsvec(:), fvec(:), p(:)
!     complex(pfdp),   pointer :: wk(:)
    type(pf_fft_t),     pointer :: fft

    integer :: l, loopstart, loopend, newtoniter
    real(pfdp) :: delta_s_norm2, damping
    
    fft => this%fft_tool
!     wk => fft%get_wk_ptr_1d()

    
    
    if (piece == 2) then
      select case (flags)
      case (0)
        loopstart = 1
        loopend = 2
      case (1)
        loopstart = 1
        loopend = 1
      case (2)
        loopstart = 2
        loopend = 2
      case default
        stop "ERROR in f2comp: only 0, 1, 2 allowed as flags"
      end select

      do l = loopstart, loopend
        s  => get_array1d_oc(y, l)
        rhsvec => get_array1d_oc(rhs, l)
        fvec => get_array1d_oc(f, l)
        call fft%conv(rhsvec,1.0_pfdp/(1.0_pfdp - dtq*this%lap),s)
        fvec = (s - rhsvec) / dtq

! !         wk => this%wk
!         wk = rhsvec
! !         call fftw_execute_dft(this%ffft, wk, wk)
! !         wk = wk / (1.0_pfdp - nu*dtq*this%lap) / size(wk)
! !         call fftw_execute_dft(this%ifft, wk, wk)
!         call fft%fftf()            
!         wk = wk / (1.0_pfdp - nu*dtq*this%lap) 
!         call fft%fftb()
!         s = real(wk)
!         fvec = (s - rhsvec) / dtq
! !         wk = s
! !         call fftw_execute_dft(this%ffft, wk, wk)
! !         wk = nu * this%lap * wk / size(wk)
! !         call fftw_execute_dft(this%ifft, wk, wk)
! !         fvec = real(wk)
      end do
    else if (piece == 3) then
      select case (flags)
      case (0)
       fvec => get_array1d_oc(f, 1)
       rhsvec => get_array1d_oc(rhs, 1)
       s  => get_array1d_oc(y, 1)
!        allocate(yold(size(s)))
!        yold=s    
       s = rhsvec/(1.0_pfdp+dtq*(gamma*1.0_pfdp/3.0_pfdp*s*s-1.0_pfdp))
       fvec = - (gamma*1.0_pfdp/3.0_pfdp*s*s-1.0_pfdp)*s
!        deallocate(yold)
       fvec => get_array1d_oc(f, 2)
       p  => get_array1d_oc(y, 2)
       s  => get_array1d_oc(y, 1)
       rhsvec => get_array1d_oc(rhs, 2)
       p = rhsvec/(1.0_pfdp + dtq*(gamma*s*s-1.0_pfdp))
       fvec = - (gamma*s*s-1.0_pfdp)*p
      case (1)
       ! solving y+dt*y*(1/3*yold*yold-1.0) = rhs
       fvec => get_array1d_oc(f, 1)
       rhsvec => get_array1d_oc(rhs, 1)
       s  => get_array1d_oc(y, 1)
       
       ! with lagging: avoid nonlinear solve
       if(lagging .eqv. .true.) then
          s = rhsvec/(1.0_pfdp+dtq*(gamma*1.0_pfdp/3.0_pfdp*s*s-1.0_pfdp))    
       else
        ! without lagging: use plain Newton to solve, damping with natural monotonicity test
         newtoniter = 0
         this%newton_f = gamma*1.0_pfdp/3.0_pfdp*dtq*s*s*s+ (1.0_pfdp-dtq)*s- rhsvec
         do           
!            if (maxval(abs(this%newton_f)) < 1e-6) exit
           this%newton_fprime = gamma*dtq*s*s+ (1.0_pfdp-dtq)
           if (maxval(abs(this%newton_fprime)) < 1e-12) then
             print *, "Newton: derivative is close to zero", &
                       maxval(abs(this%newton_f)), maxval(this%newton_fprime), minval(this%newton_fprime)
             exit
           end if
           this%delta_s = this%newton_f/this%newton_fprime
           delta_s_norm2 = sum(this%delta_s*this%delta_s)
           
           if (sqrt(delta_s_norm2) < 1e-12 ) exit
           damping = 1
           do
            s = s - damping * this%delta_s
            this%newton_f = gamma*1.0_pfdp/3.0_pfdp*dtq*s*s*s+ (1.0_pfdp-dtq)*s- rhsvec
            this%delta_sbar = this%newton_f/this%newton_fprime
!             print *, newtoniter, damping, delta_s_norm2, sum(this%delta_sbar*this%delta_sbar)
            if ( sum(this%delta_sbar*this%delta_sbar) .le. (1-damping/2.0)*delta_s_norm2 ) exit
            s = s + damping * this%delta_s
            damping = damping/2
            if(damping < 1e-7) exit
           end do
           if (damping < 1e-7) then 
            print *, "Newton: no damping parameter found", &
                      maxval(abs(this%newton_f)), maxval(this%newton_fprime), minval(this%newton_fprime)
            exit
           end if
           newtoniter = newtoniter + 1
           if (newtoniter > max_newton_iter) then
             print *, "Newton did not converge in maximum iterations", &
                      maxval(abs(this%newton_f)), maxval(this%newton_fprime), minval(this%newton_fprime)
             exit
           end if
!            s = s - this%newton_f/this%newton_fprime
         end do     
       end if
       
        fvec = - (gamma*1.0_pfdp/3.0_pfdp*s*s-1.0_pfdp)*s 
      case (2)
       ! solving p+dt*(y*y-1.0)p=rhs
       fvec => get_array1d_oc(f, 2)
       p  => get_array1d_oc(y, 2)
       s  => get_array1d_oc(y, 1)
       rhsvec => get_array1d_oc(rhs, 2)
       p = rhsvec/(1.0_pfdp + dtq*(gamma*s*s-1.0_pfdp))
       fvec = - (gamma*s*s-1.0_pfdp)*p
    case default
       stop "ERROR in f3eval: only 0, 1, 2 allowed as flags"
    end select    
    else
      print *,'Bad piece in f_comp ',piece
      call exit(0)
    end if
  end subroutine f_comp



  subroutine setup(sweeper, nvars, nnodes, nsteps, nproc)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,     intent(in)  :: nvars, nnodes, nsteps, nproc

    class(ad_sweeper_t), pointer :: this
    integer     :: i
!     type(c_ptr) :: wk
    real(pfdp)  :: kx, dx

    this => as_ad_sweeper(sweeper)
    this%use_LUq = .true.

 ! create in-place, complex fft plans
!     wk = fftw_alloc_complex(int(nvars, c_size_t))
!     call c_f_pointer(wk, this%wk, [nvars])
    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nvars],1,[Lx])

!     this%ffft = fftw_plan_dft_1d(nvars, &
!          this%wk, this%wk, FFTW_FORWARD, FFTW_ESTIMATE)
!     this%ifft = fftw_plan_dft_1d(nvars, &
!          this%wk, this%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
!     allocate(this%ddx(nvars))
    allocate(this%lap(nvars))
!     do i = 1, nvars
!        if (i <= nvars/2+1) then
!           kx = two_pi / Lx * dble(i-1)
!        else
!           kx = two_pi / Lx * dble(-nvars + i - 1)
!        end if
! 
!        this%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx
! 
!        if (kx**2 < 1e-13) then
!           this%lap(i) = 0.0_pfdp
!        else
!           this%lap(i) = -kx**2
!        end if
!     end do
    call this%fft_tool%make_lap(this%lap)
!     call this%fft_tool%make_deriv_1d(this%ddx) 

    this%nproc = nproc
    this%nsteps_per_rank = nsteps
    this%myrank = -1
    
    if(lagging .eqv. .false.) then
      allocate(this%newton_f(nvars))
      allocate(this%newton_fprime(nvars))
      allocate(this%delta_s(nvars))
      allocate(this%delta_sbar(nvars))
    end if
    
    ! allocate control and desired state
    allocate(this%u(nsteps,nnodes,nvars))
    allocate(this%ydesired(nsteps,nnodes,nvars))
    !work%u = 0
    !work%ydesired = 0
  end subroutine setup
  
  subroutine setrank(sweeper, rank)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer, intent(in) :: rank
    class(ad_sweeper_t), pointer :: this
    this => as_ad_sweeper(sweeper)
    this%myrank = rank
  end subroutine setrank


  subroutine destroy(sweeper)
    class(pf_sweeper_t), intent(inout) :: sweeper
    class(ad_sweeper_t), pointer :: this
    this => as_ad_sweeper(sweeper) 
!     deallocate(this%wk)
!     deallocate(this%ddx)
    deallocate(this%lap)
!     call fftw_destroy_plan(this%ffft)
!     call fftw_destroy_plan(this%ifft)
! 
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)
    
    if(allocated(this%newton_f))      deallocate(this%newton_f)
    if(allocated(this%newton_fprime)) deallocate(this%newton_fprime)
    if(allocated(this%delta_s))       deallocate(this%delta_s)
    if(allocated(this%delta_sbar))    deallocate(this%delta_sbar)
    
    deallocate(this%u)
    deallocate(this%ydesired)

!     call this%destroy(lev)

  end subroutine destroy


  subroutine initialize_oc(s, alpha)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp), intent(in)    :: alpha

    integer :: nnodes, i, nsteps, step
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nsteps = size(sweeper%u,1)
    nnodes = size(sweeper%u,2)
    do step=1, nsteps
      do i=1, nnodes
        sweeper%u(step,i,:) = 0.0_pfdp
      end do
    end do
    sweeper%alpha = alpha
  end subroutine initialize_oc



  subroutine fill_ydesired_nagumo(s, pf, step)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: step
    integer :: nnodes, i
    character(len=256)     :: fname
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = pf%levels(pf%nlevels)%nnodes
    do i = 1, nnodes
      call pf%levels(pf%nlevels)%Q(i)%pack(sweeper%ydesired(step,i,:), 1)
!       write(fname, "('y_d','r',i0.2,'s',i0.2,'m',i0.2,'.npy')") pf%rank, step, i
      write(fname, "('y_d','r',i0.2,'m',i0.2,'.npy')") (step-1)*sweeper%nproc + sweeper%myrank, i

!       call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
!            1, pf%levels(pf%nlevels)%lev_shape, size(sweeper%ydesired, 3), sweeper%ydesired(step,i,:))
    end do
  end subroutine fill_ydesired_nagumo



  subroutine initialize_control(s, solAt25, t0, dt, nodes)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),        intent(in   ) :: solAt25(:), nodes(:), t0, dt
    real(pfdp), pointer :: wk(:)
    integer                 :: nnodes, m, nsteps, n, thisstep
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nsteps = size(sweeper%u, 1)
    nnodes = size(nodes)
    allocate(wk(size(solAt25)))
!     wk => sweeper%wk
!     wk => sweeper%fft_tool%get_wk_ptr_1d()
     wk = solAt25
! !     call fftw_execute_dft(sweeper%ffft, wk, wk)
! !     wk = nu * sweeper%lap * wk / size(wk)
! !     call fftw_execute_dft(sweeper%ifft, wk, wk)
!     call sweeper%fft_tool%fftf()            
!     wk = nu*sweeper%lap*wk 
!     call sweeper%fft_tool%fftb()     
    call sweeper%fft_tool%conv(wk,nu*sweeper%lap,wk)
    
    do n=1, nsteps
      do m = 1, nnodes
        thisstep = (n-1)*sweeper%nproc + sweeper%myrank
        sweeper%u(n, m,:) = 0.0_pfdp
        if ( t0+thisstep*dt + dt*nodes(m) > 2.5 ) then
          sweeper%u(n,m,:)  = (gamma*1.0_pfdp/3.0_pfdp*solAt25(:)**3 - solAt25(:) - wk(:)) *0.5 
        end if
      end do
    end do
    deallocate(wk)
  end subroutine

  subroutine dump_stuff(pf,data, fbase)    
    type(pf_pfasst_t),  intent(inout) :: pf
    real(pfdp), intent(in) :: data(:,:,:)
    character(len = *), intent(in   ) :: fbase

    integer :: nnodes, i, nsteps, n
    character(len=256)     :: fname

    nsteps = size(data, 1)
    nnodes = size(data, 2)
    do n = 1, nsteps
      do i=1, nnodes
        write(fname, "(A,'s',i0.4,'m',i0.2,'.npy')") trim(fbase), n, i!(n-1)*stepper%nproc + stepper%myrank, i

       !call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
       !    1, pf%levels(pf%nlevels)%lev_shape, size(data, 3), data(n,i,:))
      end do
    end do
  end subroutine dump_stuff

  subroutine dump_control(s, pf, fbase)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t),  intent(inout) :: pf
    character(len = *), intent(in   ) :: fbase

    integer :: nnodes, i, nsteps, n
    character(len=256)     :: fname
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nsteps = size(sweeper%u, 1)
    nnodes = pf%levels(pf%nlevels)%nnodes
    do n = 1, nsteps
      do i = 1, nnodes
        write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase), (n-1)*sweeper%nproc + sweeper%myrank, i

        !call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
        !    1, pf%levels(pf%nlevels)%lev_shape, size(sweeper%u, 3), sweeper%u(n, i,:))
      end do
    end do
  end subroutine dump_control



  subroutine dump_exact_control(s, pf, solAt25, t0, dt, nodes, L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl)
    type(pf_pfasst_t),    intent(inout) :: pf
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),           intent(in   ) :: nodes(:), t0, dt, solAt25(:)
    real(pfdp),           intent(  out) :: L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl

    integer                :: nnodes, i, nvars, n, nsteps, step, thisstep
    character(len=256)     :: fname
    real(pfdp), pointer    :: uexact(:), udiff(:), ex(:), diff(:)
    real(pfdp), pointer :: wk(:)

    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nvars = product(pf%levels(pf%nlevels)%lev_shape)
    nnodes = size(nodes,1)
    nsteps = size(sweeper%u,1)
    
    
    allocate(uexact(nvars))
    allocate(udiff(nvars))
    allocate(ex(nnodes))
    allocate(diff(nnodes))

    LinfErrorCtrl = 0.0_pfdp
    LinfExactCtrl = 0.0_pfdp
    L2errorCtrl = 0.0
    L2exactCtrl = 0.0

!     wk => sweeper%wk
!     wk => sweeper%fft_tool%get_wk_ptr_1d()
    allocate(wk(size(solAt25)))
    wk = solAt25
!     call sweeper%fft_tool%fftf()            
!     wk = nu*sweeper%lap*wk 
!     call sweeper%fft_tool%fftb()   
    call sweeper%fft_tool%conv(wk,nu*sweeper%lap,wk)
    
    do step = 1, nsteps
      do i = 1, nnodes
        thisstep = (step-1)*sweeper%nproc + sweeper%myrank
        uexact = 0.0_pfdp
        if ( t0 + thisstep*dt + dt*nodes(i) > 2.5 ) then
          uexact(:)  = (1.0_pfdp/3.0_pfdp*solAt25(:)**3 - solAt25(:) - wk(:))
        end if

        write(fname, "('uexact','r',i0.2,'m',i0.2,'.npy')") (step-1)*sweeper%nproc + sweeper%myrank, i
!         call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
!             1, pf%levels(pf%nlevels)%lev_shape, nvars, uexact)

        udiff(:) = uexact(:) - sweeper%u(step,i,:)
        write(fname, "('udiff','r',i0.2,'m',i0.2,'.npy')") (step-1)*sweeper%nproc + sweeper%myrank, i
!         call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
!             1, pf%levels(pf%nlevels)%lev_shape, nvars, udiff(:))

        ! compute norm
        ex(i) = 0.0_pfdp
        diff(i) = 0.0_pfdp
        do n = 1, nvars
          ex(i) = ex(i) + uexact(n)**2.0
          diff(i) = diff(i) + udiff(n)**2.0
        end do
        ex(i) = ex(i)*Lx/(nvars)
        diff(i) = diff(i)*Lx/(nvars)

        LinfExactCtrl =  max(maxval(abs(uexact(:))), LinfExactCtrl)
        LinfErrorCtrl =  max(maxval(abs(udiff(:))), LinfErrorCtrl)
      end do

      do i=1, nnodes-1
        L2errorCtrl = L2errorCtrl + (diff(i)+diff(i+1))*(nodes(i+1)-nodes(i))*dt
        L2exactCtrl = L2exactCtrl + (ex(i)+ex(i+1))*(nodes(i+1)-nodes(i))*dt
      end do
    
    end do

    L2errorCtrl = 0.5*L2errorCtrl
    L2exactCtrl = 0.5*L2exactCtrl

    deallocate(ex)
    deallocate(diff)
    deallocate(uexact)
    deallocate(udiff)
    deallocate(wk)
  end subroutine dump_exact_control


  subroutine construct_gradient(s, grad, nodes, LinftyNormGrad, L2NormGradSq)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),          intent(inout) :: grad(:,:,:)
    real(pfdp),          intent(in)    :: nodes(:)
    real(pfdp),          intent(out)   :: LinftyNormGrad, L2NormGradSq

    integer              :: m, n, nnodes, nvars, step, nsteps
    real(pfdp),  pointer :: obj(:)
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nsteps = size(grad,1)
    nnodes = size(grad,2)
    nvars  = size(grad,3)
    allocate(obj(nnodes))

    LinftyNormGrad =  0
    L2NormGradSq = 0
    do step = 1, nsteps
      do m = 1, nnodes
        grad(step,m,:) = grad(step,m,:) + sweeper%alpha * sweeper%u(step,m,:)
        LinftyNormGrad = max(maxval(abs(grad(step,m,:))), LinftyNormGrad)
        !obj(m) = sum((grad(m,:)**2.0))*Lx/nvars
        obj(m) = 0.0_pfdp
        do n = 1, nvars
          obj(m) = obj(m) + grad(step,m,n)**2.0
        end do
        obj(m) = obj(m)*Lx/(nvars)
      end do

      do m=1, nnodes-1
        L2normGradSq = L2normGradSq + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
      end do
    end do
    
    L2NormGradSq = 0.5*L2NormGradSq !0.5 for trapezoidal rule

    deallocate(obj)
  end subroutine construct_gradient


  subroutine update_control(s, delta, stepsize)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),          intent(in)    :: delta(:,:,:), stepsize

    integer              :: m, nnodes, step, nsteps
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = size(delta,2)
    nsteps = size(delta,1)
  
    do step = 1, nsteps
      do m = 1, nnodes
        sweeper%u(step,m,:) = sweeper%u(step,m,:) + stepsize * delta(step,m,:)
      end do
    end do
 end subroutine update_control




  subroutine control_L2Q(s, dt, nodes, nvars, L2normSq)
    class(pf_sweeper_t), intent(inout) :: s
    integer,             intent(in)    :: nvars
    real(pfdp),          intent(in)    :: nodes(:), dt
    real(pfdp),          intent(out)   :: L2normSq

    real(pfdp),  pointer   :: obj(:)
    integer                :: m, n, nnodes, step, nsteps

    class(ad_sweeper_t), pointer :: sweeper;
    sweeper => as_ad_sweeper(s)


    nsteps = size(sweeper%u, 1)
    nnodes = size(nodes)
    allocate(obj(nnodes))

    L2normSq = 0.0
    
    do step=1, nsteps
      do m=1, nnodes
        obj(m) = 0.0_pfdp
        do n = 1, nvars
          obj(m) = obj(m) + sweeper%u(step,m,n)**2.0
        end do
        obj(m) = obj(m)*Lx/(nvars)
      end do

      do m=1, nnodes-1
        L2normSq = L2normSq + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
      end do
    end do

    L2normSq = 0.5*L2normSq !0.5 for trapezoidal rule

    deallocate(obj)
  end subroutine control_L2Q


  subroutine objective_function(s, sol, nvars, m, objective, step)
  ! actually just the tracking part of the objective
    class(pf_sweeper_t), intent(inout) :: s
    class(pf_encap_t), intent(in   )   :: sol
    integer,             intent(in)    :: nvars, m, step
    real(pfdp),          intent(out)   :: objective

    real(pfdp),  pointer   :: y(:), f(:) !, obj(:)
    integer                :: n !, nnodes
    complex(pfdp), pointer :: wk(:)

    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    allocate(f(nvars))

       y => get_array1d_oc(sol, 1)
       f = (y -sweeper%ydesired(step,m,:))
       objective = 0.0_pfdp
       do n = 1, nvars
         objective = objective + f(n)**2
       end do
       objective = objective * Lx / (nvars)

    deallocate(f)
  end subroutine objective_function


  ! this shouldbe in a separate module
  function compute_scalar_prod(f, g, nodes, dt) result(r)
    real(pfdp), intent(in)  :: f(:,:,:), g(:,:,:), nodes(:), dt
    real(pfdp)              :: r
    real(pfdp), pointer     :: obj(:)
    integer                 :: nnodes, nvars, m, i, step, nsteps
    r = 0.0_pfdp
    nsteps = size(f,1)
    nnodes = size(f,2)
    nvars  = size(f,3)
    allocate(obj(nnodes))
    do step = 1, nsteps
      do m = 1, nnodes
        obj(m) = 0.0_pfdp
        do i = 1, nvars
          obj(m) = obj(m) + f(step,m,i)*g(step,m,i)
        end do
        obj(m) = obj(m)*Lx/(nvars)
      end do
      do m=1, nnodes-1
        r = r + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
      end do
    end do
    r = r * 0.5
    deallocate(obj)
  end function compute_scalar_prod


  subroutine restrict_control(sG, sF)
    class(pf_sweeper_t), intent(inout) :: sG, sF

    real(pfdp), pointer  :: uF(:,:,:), uG(:,:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m, nsteps, step

    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    uF => sweeperF%u
    uG => sweeperG%u

    nsteps = size(uF,1)
    nnodesF = size(uF,2)
    nnodesG = size(uG,2)
    nvarF = size(uF,3)
    nvarG = size(uG,3)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict u', xrat, trat

    do step=1,nsteps
       uG(step,:,:) = uF(step,::trat,::xrat)
    end do
  end subroutine restrict_control


  subroutine restrict_ydesired(sG, sF)
    class(pf_sweeper_t), intent(inout) :: sG, sF

    real(pfdp), pointer  :: ydesiredF(:,:,:), ydesiredG(:,:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m, nsteps, step

    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    ydesiredF => sweeperF%ydesired
    ydesiredG => sweeperG%ydesired

    nsteps = size(ydesiredF, 1)
    nnodesF = size(ydesiredF,2)
    nnodesG = size(ydesiredG,2)
    nvarF = size(ydesiredF,3)
    nvarG = size(ydesiredG,3)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict ydesired', xrat, trat

    do step=1,nsteps
       ydesiredG(step,:,:) = ydesiredF(step,::trat,::xrat)
    end do
  end subroutine restrict_ydesired

  
    subroutine restrict_for_adjoint(pf, t0, dt, which, step)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp), intent(in) :: t0, dt
    integer, intent(in) :: which, step
    real(pfdp), pointer :: tF(:), tG(:) !zF(:), zG(:)
    integer :: l, m !, nnodesF, nnodesG, nvarsF, nvarsG

!     call pf%levels(pf%nlevels)%ulevel%sweeper%evaluate_all(pf%levels(pf%nlevels), t0+dt*pf%levels(pf%nlevels)%nodes, which, step)
    do l = pf%nlevels, 2, -1
!       allocate(tF(pf%levels(l)%nnodes))
!       allocate(tG(pf%levels(l-1)%nnodes))
!       tF = pf%state%t0 + pf%state%dt*pf%levels(l)%nodes
!       tG = pf%state%t0 + pf%state%dt*pf%levels(l-1)%nodes
      call restrict_ts(pf%levels(l), pf%levels(l-1), pf%levels(l)%Q, pf%levels(l-1)%Q, t0+dt*pf%levels(l)%nodes, which) !tF
!       call pf%levels(l-1)%ulevel%sweeper%evaluate_all(pf%levels(l-1), t0+dt*pf%levels(l-1)%nodes, which, step)
!       deallocate(tF)
!       deallocate(tG)
    end do
  end subroutine restrict_for_adjoint

  subroutine restrict_and_evaluate(pf, t0, dt, which, step)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp), intent(in) :: t0, dt
    integer, intent(in) :: which, step
    real(pfdp), pointer :: tF(:), tG(:) !zF(:), zG(:)
    integer :: l, m !, nnodesF, nnodesG, nvarsF, nvarsG

!     call pf%levels(pf%nlevels)%ulevel%sweeper%evaluate_all(pf%levels(pf%nlevels), t0+dt*pf%levels(pf%nlevels)%nodes, which, step)
    do l = pf%nlevels, 2, -1
!       allocate(tF(pf%levels(l)%nnodes))
!       allocate(tG(pf%levels(l-1)%nnodes))
!       tF = pf%state%t0 + pf%state%dt*pf%levels(l)%nodes
!       tG = pf%state%t0 + pf%state%dt*pf%levels(l-1)%nodes
!       call restrict_sdc(pf%levels(l), pf%levels(l-1), pf%levels(l)%Q, pf%levels(l-1)%Q, .false., t0+dt*pf%levels(l)%nodes, which) !tF
!       call pf%levels(l-1)%ulevel%sweeper%evaluate_all(pf%levels(l-1), t0+dt*pf%levels(l-1)%nodes, which, step)
        call pf_residual(pf, l, dt, which)  
        call restrict_time_space_fas(pf, t0, dt, l, flags=which, mystep=step)
!         call save(pf%levels(l), which)
!       deallocate(tF)
!       deallocate(tG)
    end do
  end subroutine restrict_and_evaluate


  subroutine interp1(qF, qC, adF, adC)
    class(ad_sweeper_t), pointer :: adF, adC
    real(pfdp), pointer, intent(inout) :: qF(:), qC(:)

    class(ad_sweeper_t), pointer :: sweeperF, sweeperC
!     complex(pfdp),       pointer :: wkF(:), wkC(:)
    type(pf_fft_t),      pointer :: fftF,fftC

    integer      :: NxF, NxC, xrat,i

    sweeperC => as_ad_sweeper(adC)
    sweeperF => as_ad_sweeper(adF)
    fftC => sweeperC%fft_tool
    fftF => sweeperF%fft_tool
    
    NxF = size(qF)
    NxC = size(qC)
    xrat  = NxF / NxC

    if (xrat == 1) then
       qF = qC
       return
    endif

    call fftC%interp(qC,fftF,qF)

!     wkF => fftF%get_wk_ptr_1d()
!     wkC => fftC%get_wk_ptr_1d()
! 
!     wkC = qC
!     call fftC%fftf()    
! 
!     wkF = 0.0d0
!     wkF(1:NxC/2) = wkC(1:NxC/2)
!     wkF(NxF-NxC/2+2:NxF) = wkC(NxC/2+2:NxC)
! !     wkF = wkF*2.0_pfdp
! 
!     call fftF%fftb()
!     qF = real(wkF)
  end subroutine interp1


  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: f_lev
    class(pf_level_t), intent(inout) :: c_lev
    class(pf_encap_t), intent(inout) :: f_vec, c_vec
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    integer :: which
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),          pointer :: f(:), g(:)

    which = 0
    if(present(flags)) which = flags

    adG => as_ad_sweeper(c_lev%ulevel%sweeper)
    adF => as_ad_sweeper(f_lev%ulevel%sweeper)

    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array1d_oc(f_vec,1)
      g => get_array1d_oc(c_vec,1)
      call interp1(f, g, adF, adG)
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array1d_oc(f_vec,2)
      g => get_array1d_oc(c_vec,2)
      call interp1(f, g, adF, adG)
    end if
  end subroutine interpolate


  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: f_lev
    class(pf_level_t), intent(inout) :: c_lev
    class(pf_encap_t), intent(inout) :: f_vec, c_vec
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    real(pfdp), pointer :: f(:), g(:)

    integer :: nvarF, nvarG, xrat, which
    which = 0
    if(present(flags)) which = flags

    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array1d_oc(f_vec,1)
      g => get_array1d_oc(c_vec,1)
      nvarF = size(f)
      nvarG = size(g)
      xrat  = nvarF / nvarG
      g = f(::xrat)
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array1d_oc(f_vec,2)
      g => get_array1d_oc(c_vec,2)
      nvarF = size(f)
      nvarG = size(g)
      xrat  = nvarF / nvarG
      g = f(::xrat)
    end if
  end subroutine restrict

end module feval
