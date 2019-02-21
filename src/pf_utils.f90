!! Useful subroutines that don't  fit in other modules
!
! This file is part of LIBPFASST.
!
!> Module with useful subroutines that don't  fit in other modules
module pf_mod_utils
  use pf_mod_dtype
  use pf_mod_timer
  implicit none
  
contains
  !
  !> Compute full residual at each node and measure its size
  subroutine pf_residual(pf, levind, dt, flag)
    type(pf_pfasst_t), target,  intent(inout) :: pf
    integer,            intent(in) :: levind
    real(pfdp),         intent(in)    :: dt
    integer, optional,  intent(in)    :: flag

    real(pfdp) :: res_norms(pf%levels(levind)%nnodes-1)    !!  Holds norms of residual
    real(pfdp) :: sol_norms(pf%levels(levind)%nnodes)      !!  Holds norms of solution ! for adjoint: need sol at t0 as well, not only t0+dt
    integer :: m
    type(pf_level_t), pointer   :: lev
    
    call start_timer(pf, TRESIDUAL)

    lev => pf%levels(levind)
    call lev%ulevel%sweeper%residual(pf,levind, dt, flag)

    ! compute max residual norm
    sol_norms(1) = lev%Q(1)%norm(flag) ! for adjoint
    do m = 1, lev%nnodes-1
       res_norms(m) = lev%R(m)%norm(flag)
       sol_norms(m+1) = lev%Q(m+1)%norm(flag) ! only the value at lev%nnodes is needed for forward integration, right? 
    end do

    !    lev%residual = res_norms(lev%nnodes-1)
    m = lev%nnodes  ! for usual forward integration
    if(present(flag)) then
      if(flag==2) m = 1
    end if
    lev%residual = maxval(res_norms)    
    if (sol_norms(m) > 0.0d0) then
       lev%residual_rel = lev%residual/sol_norms(m)
    else
       lev%residual_rel = 0.0d0
    end if

    if (pf%save_residuals .and. pf%state%iter>0)  then
       pf%results(lev%index)%residuals(pf%state%iter, pf%state%pfblock, pf%state%sweep) = lev%residual
    end if

    call end_timer(pf, TRESIDUAL)

  end subroutine pf_residual

  !
  !> Generic residual
  !! Each sweeper can define its own residual, or use this generic one
  !! This routine is in the "Q" form, so the residual approximates
  !! R(m)=y(t_n) + \int_{t_n}^t_m f(y,s) ds - y(t_m)
  subroutine pf_generic_residual(this, pf,levind, dt, flags)
    class(pf_sweeper_t), intent(in)  :: this
    type(pf_pfasst_t), target,  intent(inout) :: pf
    integer,              intent(in)    :: levind
    real(pfdp),        intent(in)    :: dt
    integer,  intent(in), optional  :: flags

    integer :: m
    type(pf_level_t), pointer :: lev

    lev=>pf%levels(levind)
    
    !>  Compute the integral of F from t_n to t_m at each node
    call lev%ulevel%sweeper%integrate(lev, lev%Q, lev%F, dt, lev%I, flags)

    !> add tau if it exists
    if (lev%index < pf%state%finest_level) then    
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), flags)
       end do
    end if

    !> subtract out the solution value
    if (present(flags)) then
      do m = 1, lev%nnodes-1      
        if( (flags .eq. 0) .or. (flags .eq. 1) ) then
          call lev%R(m)%copy(lev%I(m), 1)
          call lev%R(m)%axpy(1.0_pfdp, lev%Q(1), 1)
          call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1), 1)
        end if
        if( (flags .eq. 0) .or. (flags .eq. 2) ) then
          call lev%R(m)%copy(lev%I(m), 2)
          call lev%R(m)%axpy(1.0_pfdp, lev%Q(lev%nnodes), 2)
          call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m), 2)
        end if
      end do
    else
      do m = 1, lev%nnodes-1      
        call lev%R(m)%copy(lev%I(m))
        call lev%R(m)%axpy(1.0_pfdp, lev%Q(1))
        call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
      end do
    end if

  end subroutine pf_generic_residual


  !
  !> Generic evaluate all
  !! Each sweeper can define its own evaluate_all or use this generic one
  subroutine pf_generic_evaluate_all(this, lev, t, flags, step)
    class(pf_sweeper_t), intent(in)  :: this
    class(pf_level_t),  intent(inout) :: lev
    real(pfdp),        intent(in)    :: t(:)
    integer, optional, intent(in)    :: flags, step

    integer :: m
        
!     which = 1
!     if(present(flags)) which = flags
    
!     mystep = 1
!     if(present(step)) mystep = step
    
    do m = 1, lev%nnodes
       call lev%ulevel%sweeper%evaluate(lev, t(m), m, flags=flags, step=step)
    end do
  end subroutine pf_generic_evaluate_all

  
  !> Generic routine to spread initial conditions
  !! Each sweeper can define its own spreadq0 or use this generic one
  subroutine pf_generic_spreadq0(this,lev, t0)
    class(pf_sweeper_t), intent(in)  :: this
    class(pf_level_t), intent(inout) :: lev  !!  Level on which to spread
    real(pfdp),       intent(in)    :: t0    !!  time at beginning of interval

    integer :: m, p

    !  Stick initial condition into first node slot
    call lev%Q(1)%copy(lev%q0)

    !  Evaluate F at first spot
    call lev%ulevel%sweeper%evaluate(lev, t0, 1)

    ! Spread F and solution to all nodes
    do m = 2, lev%nnodes
       call lev%Q(m)%copy(lev%Q(1))
       do p = 1, lev%ulevel%sweeper%npieces
         call lev%F(m,p)%copy(lev%F(1,p))
       end do
    end do
  end subroutine pf_generic_spreadq0

  subroutine pf_stop(pf_file,Nline,msg, N)
    character(len=*), intent(in) :: pf_file
    integer, intent(in):: Nline
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: N

    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *,'Stopping in File: ', pf_file    
    print *,'Line number: ', Nline
    print *,msg
    if (present(N))   print *,'value=',N
    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    stop
    
  end subroutine pf_stop


end module pf_mod_utils
