module feval
  use pf_mod_dtype
  use pf_mod_ndarray_oc 
  use pf_mod_restrict
  use pf_mod_imexQ_oc
  use probin
  use solutions
  implicit none
  include 'fftw3.f03'

!   real(pfdp), parameter :: pi = 3.141592653589793_pfdp    !defined in probin already
!   real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp
  
  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_imexQ_oc_t) :: ad_sweeper_t  !generalize so we can use misdc as well?
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), allocatable :: ddx(:), lap(:) ! operators
     real(pfdp), pointer    :: u(:,:)
     real(pfdp), pointer    :: ydesired(:,:)
     real(pfdp)             :: alpha
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
    integer :: l, loopstart, loopend

      
    select case (piece)
    case (1)  ! Explicit piece
      select case (flags)
      case (0)
        ! first component: y
        fvec => get_array1d_oc(f, 1)
        yvec  => get_array1d_oc(y, 1)
        if (do_imex .eq. 1) then
          fvec = this%u(idx, :) - (1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
        else
          fvec = this%u(idx, :)
       endif
       ! second component: p
       fvec => get_array1d_oc(f, 2)      
       if (do_imex .eq. 1) then
         p  => get_array1d_oc(y, 2)
         fvec = (yvec - this%ydesired(idx,:)) - (yvec*yvec-1.0_pfdp)*p
       else
         fvec = (yvec - this%ydesired(idx,:))
       end if 
      case (1)
         fvec => get_array1d_oc(f, 1)
         if (do_imex .eq. 1) then
            yvec  => get_array1d_oc(y, 1)
            fvec = this%u(idx, :) - (1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
         else
            fvec = this%u(idx, :)
         endif
       case (2)
         ! evaluate y-y_d
         fvec => get_array1d_oc(f, 2)
         yvec  => get_array1d_oc(y, 1)
         if (do_imex .eq. 1) then
           p  => get_array1d_oc(y, 2)
           fvec = (yvec -this%ydesired(idx,:)) - (yvec*yvec-1.0_pfdp)*p
         else
           fvec = (yvec -this%ydesired(idx,:)) 
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
        wk => this%wk
        wk = yvec
        call fftw_execute_dft(this%ffft, wk, wk)
        wk = nu * this%lap * wk / size(wk)
        call fftw_execute_dft(this%ifft, wk, wk)
        fvec = real(wk)
       end do
      
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select
  end subroutine f_eval

  
!   subroutine f1eval1(this, sol, t, level, f, flags, idx, step)
!     class(ad_sweeper_t), intent(inout) :: this
!     class(pf_encap_t),   intent(in   ) :: sol
!     class(pf_encap_t),   intent(inout) :: f
!     real(pfdp),  intent(in)           :: t
!     integer,     intent(in)           :: level
!     integer,     intent(in)           :: flags
!     integer,     intent(in)           :: idx      ! index of the current quadrature node
!     integer,     intent(in), optional :: step
! 
!     real(pfdp),    pointer :: y(:), f1(:), p(:)
! 
!     select case (flags)
!     case (0)
!        ! first component: y
!        f1 => get_array1d_oc(f, 1)
!        y  => get_array1d_oc(sol, 1)
!        
!        if (do_imex .eq. 1) then
!          f1 = this%u(idx, :) - (1.0_pfdp/3.0_pfdp*y*y-1.0_pfdp)*y 
!        else
!          f1 = this%u(idx, :)
!        endif
! 
!        ! second component: p
!        f1 => get_array1d_oc(f, 2)      
!        if (do_imex .eq. 1) then
!          p  => get_array1d_oc(sol, 2)
!          f1 = (y - this%ydesired(idx,:)) - (y*y-1.0_pfdp)*p
!        else
!          f1 = (y - this%ydesired(idx,:))
!        end if 
! 
!     case (1)
!        f1 => get_array1d_oc(f, 1)
!        if (do_imex .eq. 1) then
!          y  => get_array1d_oc(sol, 1)
!          f1 = this%u(idx, :) - (1.0_pfdp/3.0_pfdp*y*y-1.0_pfdp)*y 
!        else
!          f1 = this%u(idx, :)
!        endif
! 
!     case (2)
!        ! evaluate y-y_d
!        f1 => get_array1d_oc(f, 2)
!        y  => get_array1d_oc(sol, 1)
!        if (do_imex .eq. 1) then
!            p  => get_array1d_oc(sol, 2)
!            f1 = (y -this%ydesired(idx,:)) - (y*y-1.0_pfdp)*p
!          else
!            f1 = (y -this%ydesired(idx,:)) 
!          end if
!     case default
!        stop "ERROR in f1eval1: only 0, 1, 2 allowed as flags"
!     end select
!   end subroutine f1eval1
! 
!   
!   subroutine f2eval1(this, sol, t, level, f, flags)
!     class(ad_sweeper_t), intent(inout) :: this
!     class(pf_encap_t),   intent(in   ) :: sol
!     class(pf_encap_t),   intent(inout) :: f
!     real(pfdp),          intent(in)    :: t
!     integer,             intent(in)    :: level
!     integer,             intent(in)    :: flags
! 
!     real(pfdp),    pointer :: s(:), f2(:)
!     complex(pfdp), pointer :: wk(:)
!     integer ::  l, loopstart, loopend
! 
!     select case (flags)
!     case (0)
!        loopstart = 1
!        loopend = 2
!     case (1)
!        loopstart = 1
!        loopend = 1
!     case (2)
!        loopstart = 2
!        loopend = 2
!     case default
!        stop "ERROR in f2eval1: only 0, 1, 2 allowed as flags"
!     end select
!        
!     do l = loopstart, loopend	
!        s   => get_array1d_oc(sol, l)
!        f2  => get_array1d_oc(f, l)
!        wk => this%wk
!        wk = s
!        call fftw_execute_dft(this%ffft, wk, wk)
!        wk = nu * work%lap * wk / size(wk)
!        call fftw_execute_dft(this%ifft, wk, wk)
!        f2 = real(wk)
!     end do
!   end subroutine f2eval1

    
    
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

    real(pfdp),      pointer :: s(:), rhsvec(:), fvec(:)
    complex(pfdp),   pointer :: wk(:)
    integer :: l, loopstart, loopend
    
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
        wk => this%wk
      
        wk = rhsvec
        call fftw_execute_dft(this%ffft, wk, wk)
        wk = wk / (1.0_pfdp - nu*dtq*this%lap) / size(wk)
        call fftw_execute_dft(this%ifft, wk, wk)
        s = real(wk)
        fvec = (s - rhsvec) / dtq
!         wk = s
!         call fftw_execute_dft(this%ffft, wk, wk)
!         wk = nu * this%lap * wk / size(wk)
!         call fftw_execute_dft(this%ifft, wk, wk)
!         fvec = real(wk)
      end do
    else
      print *,'Bad piece in f_comp ',piece
      call exit(0)
    end if
  end subroutine f_comp


  
  subroutine setup(sweeper, nvars, nnodes)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,     intent(in)  :: nvars, nnodes

    class(ad_sweeper_t), pointer :: this
    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx, dx

    this => as_ad_sweeper(sweeper)

 ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nvars, c_size_t))
    call c_f_pointer(wk, this%wk, [nvars])

    this%ffft = fftw_plan_dft_1d(nvars, &
         this%wk, this%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    this%ifft = fftw_plan_dft_1d(nvars, &
         this%wk, this%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
    allocate(this%ddx(nvars))
    allocate(this%lap(nvars))
    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-nvars + i - 1)
       end if

       this%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          this%lap(i) = 0.0_pfdp
       else
          this%lap(i) = -kx**2
       end if
    end do
    
    ! allocate control and desired state
    allocate(this%u(nnodes,nvars))
    allocate(this%ydesired(nnodes,nvars))
    !work%u = 0
    !work%ydesired = 0  
  end subroutine setup

  
  subroutine destroy(this, lev)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev

    deallocate(this%wk)
    deallocate(this%ddx)
    deallocate(this%lap)
    call fftw_destroy_plan(this%ffft)
    call fftw_destroy_plan(this%ifft)
    
    deallocate(this%u)
    deallocate(this%ydesired)
    
    call this%imexQ_oc_destroy(lev)

  end subroutine destroy


  subroutine initialize_oc(s, t0, dt, nodes, nvars, alpha)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp), intent(in)    :: nodes(:), t0, dt, alpha
    integer,             intent(in)    :: nvars

    integer :: nnodes, i
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nnodes = size(nodes)

    do i = 1, nnodes
       sweeper%u(i,:) = 0.0_pfdp
!        call exact_rhs(t0+dt*nodes(i), nvars, work%u(i,:))
!        work%u(i,:) = work%u(i,:)*0.5
    end do
    sweeper%alpha = alpha
  end subroutine initialize_oc
  
  
  
  subroutine fill_ydesired_nagumo(s, pf)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t), intent(inout) :: pf
    integer :: nnodes, i
    character(len=256)     :: fname
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = pf%levels(pf%nlevels)%nnodes
    do i = 1, nnodes
      call pf%levels(pf%nlevels)%Q(i)%pack(sweeper%ydesired(i,:), 1)
      write(fname, "('y_d','r',i0.2,'m',i0.2,'.npy')") pf%rank, i

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           1, pf%levels(pf%nlevels)%shape, size(sweeper%ydesired, 2), sweeper%ydesired(i,:))
    end do
  end subroutine fill_ydesired_nagumo
  
  
  
  subroutine initialize_control(s, solAt25, t0, dt, nodes)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),        intent(in   ) :: solAt25(:), nodes(:), t0, dt
    complex(pfdp), pointer :: wk(:)
    integer                :: nnodes, m
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nnodes = size(nodes)
       
    do m = 1, nnodes
       sweeper%u(m,:) = 0.0_pfdp
       if ( t0+dt*nodes(m) > 2.5 ) then
         wk => sweeper%wk
         wk = solAt25
         call fftw_execute_dft(sweeper%ffft, wk, wk)
         wk = nu * sweeper%lap * wk / size(wk)
         call fftw_execute_dft(sweeper%ifft, wk, wk)
         sweeper%u(m,:)  = (1.0_pfdp/3.0_pfdp*solAt25(:)**3 - solAt25(:) -real(wk)) *0.5
       end if
    end do
  end subroutine


  
  subroutine dump_control(s, pf, fbase)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t),  intent(inout) :: pf
    character(len = *), intent(in   ) :: fbase

    integer :: nnodes, i
    character(len=256)     :: fname
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = pf%levels(pf%nlevels)%nnodes
    do i = 1, nnodes
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase), pf%rank, i

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           1, pf%levels(pf%nlevels)%shape, size(sweeper%u, 2), sweeper%u(i,:))
    end do
  end subroutine dump_control
  
  

  subroutine dump_exact_control(s, pf, solAt25, t0, dt, nodes, L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl)
    type(pf_pfasst_t),    intent(inout) :: pf
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),           intent(in   ) :: nodes(:), t0, dt, solAt25(:)
    real(pfdp),           intent(  out) :: L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl

    integer                :: nnodes, i, nvars, n
    character(len=256)     :: fname
    real(pfdp), pointer    :: uexact(:), udiff(:), ex(:), diff(:)
    complex(pfdp), pointer :: wk(:)

    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nvars = product(pf%levels(pf%nlevels)%shape)
    nnodes = size(nodes,1)
    
    allocate(uexact(nvars))   
    allocate(udiff(nvars))   
    allocate(ex(nnodes))
    allocate(diff(nnodes))   
    
    LinfErrorCtrl = 0.0_pfdp
    LinfExactCtrl = 0.0_pfdp
    
    do i = 1, nnodes
!       call exact_rhs(t0+dt*nodes(i), nvars, uexact)
      uexact = 0.0_pfdp
      if ( t0+dt*nodes(i) > 2.5 ) then
        wk => sweeper%wk
        wk = solAt25
        call fftw_execute_dft(sweeper%ffft, wk, wk)
        wk = nu * sweeper%lap * wk / size(wk)
        call fftw_execute_dft(sweeper%ifft, wk, wk)
        uexact(:)  = (1.0_pfdp/3.0_pfdp*solAt25(:)**3 - solAt25(:) -real(wk))
      end if

      write(fname, "('uexact','r',i0.2,'m',i0.2,'.npy')") pf%rank, i
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           1, pf%levels(pf%nlevels)%shape, nvars, uexact)
           
      udiff(:) = uexact(:) - sweeper%u(i,:)
      write(fname, "('udiff','r',i0.2,'m',i0.2,'.npy')") pf%rank, i
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           1, pf%levels(pf%nlevels)%shape, nvars, udiff(:))
           
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

    L2errorCtrl = 0.0
    L2exactCtrl = 0.0
    do i=1, nnodes-1
       L2errorCtrl = L2errorCtrl + (diff(i)+diff(i+1))*(nodes(i+1)-nodes(i))*dt
       L2exactCtrl = L2exactCtrl + (ex(i)+ex(i+1))*(nodes(i+1)-nodes(i))*dt
    end do

    L2errorCtrl = 0.5*L2errorCtrl
    L2exactCtrl = 0.5*L2exactCtrl

    deallocate(ex)
    deallocate(diff)          
    deallocate(uexact)
    deallocate(udiff)          
  end subroutine dump_exact_control
  
  
  subroutine construct_gradient(s, grad, nodes, LinftyNormGrad, L2NormGradSq)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),          intent(inout) :: grad(:,:)
    real(pfdp),          intent(in)    :: nodes(:)
    real(pfdp),          intent(out)   :: LinftyNormGrad, L2NormGradSq
    
    integer              :: m, n, nnodes, nvars
    real(pfdp),  pointer :: obj(:)
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = size(grad,1)
    nvars  = size(grad,2)
    allocate(obj(nnodes))

    LinftyNormGrad =  0
    L2NormGradSq = 0
    do m = 1, nnodes
       grad(m,:) = grad(m,:) + sweeper%alpha * sweeper%u(m,:)
       LinftyNormGrad = max(maxval(abs(grad(m,:))), LinftyNormGrad)  
       !obj(m) = sum((grad(m,:)**2.0))*Lx/nvars
       obj(m) = 0.0_pfdp
       do n = 1, nvars
         obj(m) = obj(m) + grad(m,n)**2.0
       end do
       obj(m) = obj(m)*Lx/(nvars)
     end do

    L2NormGradSq = 0.0
    do m=1, nnodes-1
       L2normGradSq = L2normGradSq + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
    end do

    L2NormGradSq = 0.5*L2NormGradSq !0.5 for trapezoidal rule

    deallocate(obj)
  end subroutine construct_gradient
  

  subroutine update_control(s, delta, step)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),          intent(in)    :: delta(:,:), step
    
    integer              :: m, nnodes   
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nnodes = size(delta,1)

    do m = 1, nnodes
       sweeper%u(m,:) = sweeper%u(m,:) + step * delta(m,:)
    end do
 end subroutine update_control


 
  subroutine write_control(s, k, fbase)
    class(pf_sweeper_t), intent(inout) :: s
    integer,            intent(in) :: k
    character(len = *), intent(in) :: fbase

    character(len = 64)  :: fname
    integer              :: nvars, i
    real(pfdp)           :: dx
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nvars = size(sweeper%u, 2)
    dx = 1.0_pfdp / nvars

    write(fname, "(A,I0.6,A)") trim(fbase), k, ".txt"
    fname = trim(fname)
    open(unit=106, file = fname, status = 'unknown', action = 'write') 
    do i=1,nvars
       write(106,*) (i-1)*dx, sweeper%u(:, i)
    end do
    close(106)
  end subroutine write_control



  subroutine control_L2Q(s, dt, nodes, nvars, L2normSq)
    class(pf_sweeper_t), intent(inout) :: s
    integer,             intent(in)    :: nvars
    real(pfdp),          intent(in)    :: nodes(:), dt
    real(pfdp),          intent(out)   :: L2normSq
    
    real(pfdp),  pointer   :: obj(:)
    integer                :: m, n, nnodes
 
    class(ad_sweeper_t), pointer :: sweeper;
    sweeper => as_ad_sweeper(s)

    
    nnodes = size(nodes)
    allocate(obj(nnodes))

    do m=1, nnodes
       !print *, 'sum u/nvars = ', sum(work%u(m,:))*Lx/nvars, 'sum u2/nvars = ', sum(work%u(m,:)**2.0)*Lx/nvars, &
       !          maxval(work%u(m,:)), minval(work%u(m,:))
       obj(m) = 0.0_pfdp
       !obj(m) = sum(work%u(m,:)**2.0)*Lx/nvars !rectangle rule
       do n = 1, nvars
         obj(m) = obj(m) + sweeper%u(m,n)**2.0
       end do
       obj(m) = obj(m)*Lx/(nvars)
    end do

    L2normSq = 0.0
    do m=1, nnodes-1
       L2normSq = L2normSq + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
    end do

    L2normSq = 0.5*L2normSq !0.5 for trapezoidal rule

    deallocate(obj)
  end subroutine control_L2Q

  
  subroutine objective_function(s, sol, nvars, m, objective) 
  ! actually just the tracking part of the objective
    class(pf_sweeper_t), intent(inout) :: s
    class(pf_encap_t), intent(in   )   :: sol
    integer,             intent(in)    :: nvars, m
    real(pfdp),          intent(out)   :: objective
    
    real(pfdp),  pointer   :: y(:), f(:) !, obj(:)
    integer                :: n !, nnodes
    complex(pfdp), pointer :: wk(:)

    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    allocate(f(nvars))

       y => get_array1d_oc(sol, 1)
       f = (y -sweeper%ydesired(m,:))   
       objective = 0.0_pfdp
       do n = 1, nvars
         objective = objective + f(n)**2
       end do
       objective = objective * Lx / (nvars)

    deallocate(f)
  end subroutine objective_function
  
  
  ! this shouldbe in a separate module
  function compute_scalar_prod(f, g, nodes, dt) result(r)
    real(pfdp), intent(in)  :: f(:,:), g(:,:), nodes(:), dt
    real(pfdp)              :: r
    real(pfdp), pointer     :: obj(:)
    integer                 :: nnodes, nvars, m, i
    r = 0.0_pfdp
    nnodes = size(f,1)
    nvars  = size(f,2)
    allocate(obj(nnodes))
    do m = 1, nnodes
      obj(m) = 0.0_pfdp
      do i = 1, nvars
        obj(m) = obj(m) + f(m,i)*g(m,i) 
      end do
      obj(m) = obj(m)*Lx/(nvars)
    end do
    do m=1, nnodes-1
      r = r + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
    end do
    r = r * 0.5
    deallocate(obj)
  end function compute_scalar_prod

  
  subroutine restrict_control(sG, sF)
    class(pf_sweeper_t), intent(inout) :: sG, sF
        
    real(pfdp), pointer  :: uF(:,:), uG(:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m
    
    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    uF => sweeperF%u
    uG => sweeperG%u

    nnodesF = size(uF,1)
    nnodesG = size(uG,1)
    nvarF = size(uF,2)
    nvarG = size(uG,2)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict u', xrat, trat

    !do m=1,nnodesG
       uG(:,:) = uF(::trat,::xrat)
    !end do
  end subroutine restrict_control

  
  subroutine restrict_ydesired(sG, sF)
    class(pf_sweeper_t), intent(inout) :: sG, sF
    
    real(pfdp), pointer  :: ydesiredF(:,:), ydesiredG(:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m
    
    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    ydesiredF => sweeperF%ydesired
    ydesiredG => sweeperG%ydesired

    nnodesF = size(ydesiredF,1)
    nnodesG = size(ydesiredG,1)
    nvarF = size(ydesiredF,2)
    nvarG = size(ydesiredG,2)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict ydesired', xrat, trat

    !do m=1,nnodesG
       ydesiredG(:,:) = ydesiredF(::trat,::xrat)
    !end do
  end subroutine restrict_ydesired
  
  
  subroutine restrict_for_adjoint(pf, which)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: which
    real(pfdp), pointer :: tF(:), tG(:) !zF(:), zG(:)
    integer :: l, m !, nnodesF, nnodesG, nvarsF, nvarsG
    
    do l = pf%nlevels, 2, -1
      allocate(tF(pf%levels(l)%nnodes))
      allocate(tG(pf%levels(l-1)%nnodes))
      tF = pf%state%t0 + pf%state%dt*pf%levels(l)%nodes
      tG = pf%state%t0 + pf%state%dt*pf%levels(l-1)%nodes
      call restrict_sdc(pf%levels(l), pf%levels(l-1), pf%levels(l)%Q, pf%levels(l-1)%Q, .false. ,tF, which)

        call pf%levels(l-1)%ulevel%sweeper%evaluate_all(pf%levels(l-1), tG, which)
      deallocate(tF)
      deallocate(tG)
    end do
  end subroutine restrict_for_adjoint
  
    
  subroutine interp1(qF, qG, adF, adG)
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),  intent(inout) :: qF(:), qG(:)

    complex(pfdp), pointer :: wkF(:), wkG(:)
    integer      :: NxF, NxG, xrat,i
    
    NxF = size(qF)
    NxG = size(qG)
    xrat  = NxF / NxG

    if (xrat == 1) then
       qF = qG
       return
    endif
      
    wkF => adF%wk
    wkG => adG%wk
       
    wkG = qG
    call fftw_execute_dft(adG%ffft, wkG, wkG)
    wkG = wkG / NxG
       
    wkF = 0.0d0
    wkF(1:NxG/2) = wkG(1:NxG/2)
    wkF(NxF-NxG/2+2:NxF) = wkG(NxG/2+2:NxG)
       
    call fftw_execute_dft(adF%ifft, wkF, wkF)
       
    qF = real(wkF)    
  end subroutine interp1
    
    
  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags
    
    integer :: which
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),          pointer :: f(:), g(:)
             
    which = 0
    if(present(flags)) which = flags

    adG => as_ad_sweeper(levelG%ulevel%sweeper)
    adF => as_ad_sweeper(levelF%ulevel%sweeper)
  
    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array1d_oc(qF,1)
      g => get_array1d_oc(qG,1)
      call interp1(f, g, adF, adG)  
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array1d_oc(qF,2)
      g => get_array1d_oc(qG,2)
      call interp1(f, g, adF, adG)
    end if
  end subroutine interpolate
  

  subroutine restrict(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    real(pfdp), pointer :: f(:), g(:)

    integer :: nvarF, nvarG, xrat, which
    which = 0
    if(present(flags)) which = flags
        
    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array1d_oc(qF,1)
      g => get_array1d_oc(qG,1)
      nvarF = size(f)
      nvarG = size(g)
      xrat  = nvarF / nvarG
      g = f(::xrat)
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array1d_oc(qF,2)
      g => get_array1d_oc(qG,2)
      nvarF = size(f)
      nvarG = size(g)
      xrat  = nvarF / nvarG
      g = f(::xrat)
    end if    
  end subroutine restrict
 
end module feval
