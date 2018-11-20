module feval
  use pf_mod_dtype
  use pf_mod_ndarray_oc
  use pf_mod_restrict
  use pf_mod_imexQ_oc
  use pf_mod_misdcQ_oc
  use probin
  use solutions
  use pf_mod_fftpackage
  
  implicit none


  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

 type, extends(pf_imexQ_oc_t) :: ad_sweeper_t  !generalize so we can use misdc as well?
!  type, extends(pf_misdcQ_oc_t) :: ad_sweeper_t  !generalize so we can use misdc as well?
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: ddx(:), lap(:) ! operators
     real(pfdp), pointer    :: u(:,:)
     real(pfdp), pointer    :: ydesired(:,:)
     real(pfdp)             :: alpha
   contains

     procedure :: f_eval
     procedure :: f_comp
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
!             print *, yvec
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
        wk => this%fft_tool%get_wk_ptr_1d()
        wk = yvec
        call this%fft_tool%fftf()
        wk = nu * this%lap * wk 
        call this%fft_tool%fftb()
        fvec=real(wk)
       end do

    case (3) ! third piece (misdc)
    select case (flags)
      case (0)
        ! first component: y
        fvec => get_array1d_oc(f, 1)
        yvec  => get_array1d_oc(y, 1)
        fvec = - (1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
        ! second component: p
        fvec => get_array1d_oc(f, 2)
        p  => get_array1d_oc(y, 2)
        fvec = - (yvec*yvec-1.0_pfdp)*p
      case (1)
        fvec => get_array1d_oc(f, 1)
        yvec  => get_array1d_oc(y, 1)
        fvec = - (1.0_pfdp/3.0_pfdp*yvec*yvec-1.0_pfdp)*yvec
      case (2)
        fvec => get_array1d_oc(f, 2)
        yvec  => get_array1d_oc(y, 1)
        p  => get_array1d_oc(y, 2)
        fvec = - (yvec*yvec-1.0_pfdp)*p
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

    real(pfdp),      pointer :: s(:), rhsvec(:), fvec(:), yold(:), p(:)
    complex(pfdp),   pointer :: wk(:)
    type(pf_fft_t),     pointer :: fft
    integer :: l, loopstart, loopend

    fft => this%fft_tool
    wk => fft%get_wk_ptr_1d()
    
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

        wk = rhsvec
        call fft%fftf()            
        wk = wk / (1.0_pfdp - nu*dtq*this%lap) 
        call fft%fftb()
        s = real(wk)
        fvec = (s - rhsvec) / dtq
      end do
    else if (piece == 3) then
      select case (flags)
      case (0)
       fvec => get_array1d_oc(f, 1)
       rhsvec => get_array1d_oc(rhs, 1)
       s  => get_array1d_oc(y, 1)
       allocate(yold(size(s)))
       yold=s    
       s = rhsvec/(1.0_pfdp+dtq*(1.0_pfdp/3.0_pfdp*yold*yold-1.0_pfdp))
       fvec = - (1.0_pfdp/3.0_pfdp*s*s-1.0_pfdp)*s
       deallocate(yold)
       fvec => get_array1d_oc(f, 2)
       p  => get_array1d_oc(y, 2)
       s  => get_array1d_oc(y, 1)
       rhsvec => get_array1d_oc(rhs, 2)
       p = rhsvec/(1.0_pfdp + dtq*(s*s-1.0_pfdp))
       fvec = - (s*s-1.0_pfdp)*p
    case (1)
       ! solving y+dt*y*(1/3*yold*yold-1.0) = rhs
       fvec => get_array1d_oc(f, 1)
       rhsvec => get_array1d_oc(rhs, 1)
       s  => get_array1d_oc(y, 1)
       allocate(yold(size(s)))
       yold=s
       s = rhsvec/(1.0_pfdp+dtq*(1.0_pfdp/3.0_pfdp*yold*yold-1.0_pfdp))
       fvec = - (1.0_pfdp/3.0_pfdp*s*s-1.0_pfdp)*s
       deallocate(yold)
   case (2)
       ! solving p+dt*(y*y-1.0)p=rhs
       fvec => get_array1d_oc(f, 2)
       p  => get_array1d_oc(y, 2)
       s  => get_array1d_oc(y, 1)
       rhsvec => get_array1d_oc(rhs, 2)
       p = rhsvec/(1.0_pfdp + dtq*(s*s-1.0_pfdp))
       fvec = - (s*s-1.0_pfdp)*p
    case default
       stop "ERROR in f3eval: only 0, 1, 2 allowed as flags"
    end select    
    else
      print *,'Bad piece in f_comp ',piece
      call exit(0)
    end if
  end subroutine f_comp


  subroutine sweeper_setup(sweeper, nvars, nnodes)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,     intent(in)  :: nvars, nnodes

    class(ad_sweeper_t), pointer :: this
    integer     :: i
    real(pfdp)  :: kx, dx

    this => as_ad_sweeper(sweeper)
    this%use_LUq = .true.

    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nvars],1,[sizex])

    ! create operators
    allocate(this%ddx(nvars))
    allocate(this%lap(nvars))
    call this%fft_tool%make_lap_1d(this%lap)
    call this%fft_tool%make_deriv_1d(this%ddx) 

    ! allocate control and desired state
    allocate(this%u(nnodes,nvars))
    allocate(this%ydesired(nnodes,nvars))
    !work%u = 0
    !work%ydesired = 0
  end subroutine sweeper_setup


  subroutine sweeper_destroy(sweeper)
    class(pf_sweeper_t), intent(inout) :: sweeper

    class(ad_sweeper_t), pointer :: this
    this => as_ad_sweeper(sweeper) 
    
    deallocate(this%ddx)
    deallocate(this%lap)


    deallocate(this%u)
    deallocate(this%ydesired)

    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)  
    
  end subroutine sweeper_destroy


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

    wk => sweeper%fft_tool%get_wk_ptr_1d()
    do m = 1, nnodes
       sweeper%u(m,:) = 0.0_pfdp
       if ( t0+dt*nodes(m) > 2.5 ) then
       wk = solAt25
       call sweeper%fft_tool%fftf()            
       wk = nu*sweeper%lap*wk 
       call sweeper%fft_tool%fftb()        
       sweeper%u(m,:)  = (1.0_pfdp/3.0_pfdp*solAt25(:)**3 - solAt25(:) -real(wk)) *0.5
       end if
    end do
  end subroutine

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

    wk => sweeper%fft_tool%get_wk_ptr_1d()
    do i = 1, nnodes
!       call exact_rhs(t0+dt*nodes(i), nvars, uexact)
      uexact = 0.0_pfdp
      if ( t0+dt*nodes(i) > 2.5 ) then
        wk = solAt25
        call sweeper%fft_tool%fftf()            
        wk = nu*sweeper%lap*wk 
        call sweeper%fft_tool%fftb()     
        uexact(:)  = (1.0_pfdp/3.0_pfdp*solAt25(:)**3 - solAt25(:) -real(wk))
      end if

      udiff(:) = uexact(:) - sweeper%u(i,:)

      ! compute norm
      ex(i) = 0.0_pfdp
      diff(i) = 0.0_pfdp
      do n = 1, nvars
         ex(i) = ex(i) + uexact(n)**2.0
         diff(i) = diff(i) + udiff(n)**2.0
      end do
      ex(i) = ex(i)*sizex/(nvars)
      diff(i) = diff(i)*sizex/(nvars)

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
       !obj(m) = sum((grad(m,:)**2.0))*sizex/nvars
       obj(m) = 0.0_pfdp
       do n = 1, nvars
         obj(m) = obj(m) + grad(m,n)**2.0
       end do
       obj(m) = obj(m)*sizex/(nvars)
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
       !print *, 'sum u/nvars = ', sum(work%u(m,:))*sizex/nvars, 'sum u2/nvars = ', sum(work%u(m,:)**2.0)*sizex/nvars, &
       !          maxval(work%u(m,:)), minval(work%u(m,:))
       obj(m) = 0.0_pfdp
       !obj(m) = sum(work%u(m,:)**2.0)*sizex/nvars !rectangle rule
       do n = 1, nvars
         obj(m) = obj(m) + sweeper%u(m,n)**2.0
       end do
       obj(m) = obj(m)*sizex/(nvars)
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
       objective = objective * sizex / (nvars)

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
      obj(m) = obj(m)*sizex/(nvars)
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


  subroutine interp1(qF, qC, adF, adC)
    class(ad_sweeper_t), pointer :: adF, adC
    real(pfdp),  intent(inout) :: qF(:), qC(:)

    class(ad_sweeper_t), pointer :: sweeperF, sweeperC
    complex(pfdp),       pointer :: wkF(:), wkC(:)
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

    wkF => fftF%get_wk_ptr_1d()
    wkC => fftC%get_wk_ptr_1d()

    wkC = qC
    call fftC%fftf()    

    wkF = 0.0d0
    wkF(1:NxC/2) = wkC(1:NxC/2)
    wkF(NxF-NxC/2+2:NxF) = wkC(NxC/2+2:NxC)
    wkF = wkF*2.0_pfdp

    call fftF%fftb()
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
