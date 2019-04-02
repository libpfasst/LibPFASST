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
  subroutine pf_residual(pf, level_index, dt, flag)
    type(pf_pfasst_t), target,  intent(inout) :: pf
    integer,            intent(in) :: level_index
    real(pfdp),         intent(in)    :: dt
    integer, optional,  intent(in)    :: flag

    real(pfdp) :: res_norms(pf%levels(level_index)%nnodes-1)    !!  Holds norms of residual
    real(pfdp) :: sol_norms(pf%levels(level_index)%nnodes)      !!  Holds norms of solution ! for adjoint: need sol at t0 as well, not only t0+dt
    integer :: m
    type(pf_level_t), pointer   :: lev
    
    call start_timer(pf, TRESIDUAL)

    lev => pf%levels(level_index)
    call lev%ulevel%sweeper%residual(pf,level_index, dt, flag)

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
  subroutine pf_generic_residual(this, pf,level_index, dt, flags)
    class(pf_sweeper_t), intent(in)  :: this
    type(pf_pfasst_t), target,  intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),        intent(in)    :: dt
    integer,  intent(in), optional  :: flags

    integer :: m
    type(pf_level_t), pointer :: lev
    lev=>pf%levels(level_index)
    
    !>  Compute the integral of F from t_n to t_m at each node
    call lev%ulevel%sweeper%integrate(pf,level_index, lev%Q, lev%F, dt, lev%I, flags)

    !> add tau if it exists
    if (lev%index < pf%state%finest_level) then    
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), flags)
!          if (m>1 .and. pf%use_Sform) call lev%I(m)%axpy(-1.0_pfdp, lev%tauQ(m-1), flags)
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
  subroutine pf_generic_evaluate_all(this, pf, level_index, t, flags, step)
    class(pf_sweeper_t), intent(in)  :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index
    real(pfdp),        intent(in)    :: t(:)
    integer, optional, intent(in)    :: flags, step

    integer :: m
    class(pf_level_t), pointer :: lev    !!  points to current level
    lev => pf%levels(level_index)   !!  Assign level pointer
!     which = 1
!     if(present(flags)) which = flags
    
!     mystep = 1
!     if(present(step)) mystep = step
    
    do m = 1, lev%nnodes
       call lev%ulevel%sweeper%evaluate(pf,level_index, t(m), m, flags=flags, step=step)
    end do
  end subroutine pf_generic_evaluate_all

  
  !> Generic routine to spread initial conditions
  !! Each sweeper can define its own spreadq0 or use this generic one
  subroutine pf_generic_spreadq0(this,pf,level_index, t0)
    class(pf_sweeper_t), intent(in)  :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to sweep
    real(pfdp),       intent(in)    :: t0    !!  time at beginning of interval

    integer :: m, p
    class(pf_level_t), pointer :: lev  !!  Level on which to spread
    lev => pf%levels(level_index)   !!  Assign level pointer
    
    !  Stick initial condition into first node slot
    call lev%Q(1)%copy(lev%q0)

    !  Evaluate F at first spot
    call lev%ulevel%sweeper%evaluate(pf,level_index, t0, 1)

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



  subroutine pf_apply_mat(dst, a, mat, src, zero_first, flags)
    !! Apply a matrix (tmat or rmat) to src and add to dst.
    !! Mathematically this is 
    !!     dst= dst + a*mat*src
    !!  Where dst and src are vectors, mat is a matrix, and a is a scalar
    !!  If the optional variable "zero" is provided and is true, then we compute
    !!     dst=  a*mat*src
    class(pf_encap_t), intent(inout) :: dst(:)       !!  destination vector
    real(pfdp),        intent(in)    :: a            !!  scalar
    real(pfdp),        intent(in)    :: mat(:, :)    !!  matrix
    class(pf_encap_t), intent(in)    :: src(:)       !!  src vector
    logical,           intent(in), optional :: zero_first   !! If true, zero out the the dst variable before computing 
    integer,           intent(in), optional :: flags  !! Used for choosing which variable to operate on 
    
    !!  Local variables
    logical :: lzero   !!  local version of input parameter zero
    integer :: which   !!  local version of flags
    integer :: n, m    !!  size of mat   
    integer :: i, j    !!  loop variables

    lzero = .true.; if (present(zero_first)) lzero = zero_first
    which = 1;      if(present(flags)) which = flags
        
    n = size(mat, dim=1)
    m = size(mat, dim=2)
        
    do i = 1, n
      if (lzero) call dst(i)%setval(0.0_pfdp, flags)
      do j = 1, m
         if (abs(a*mat(i, j)) /= 0.0_pfdp)  call dst(i)%axpy(a * mat(i, j), src(j), flags)
      end do
    end do
  end subroutine pf_apply_mat
  

  subroutine pf_apply_mat_backward(dst, a, mat, src, zero_first, flags)
    !! Apply a matrix (tmat or rmat) to src and add to dst.
    class(pf_encap_t), intent(inout) :: dst(:)       !!  destination vector
    real(pfdp),        intent(in)    :: a            !!  scalar
    real(pfdp),        intent(in)    :: mat(:, :)    !!  matrix
    class(pf_encap_t), intent(in)    :: src(:)       !!  src vector
    logical,           intent(in), optional :: zero_first   !! If true, zero out the the dst variable before computing 
    integer,           intent(in), optional :: flags  !! Used for choosing which variable to operate on 

    
    !!  Local variables
    logical :: lzero   !!  local version of input parameter zero
    integer :: which   !!  local version of flags
    integer :: n, m    !!  size of mat   
    integer :: i, j    !!  loop variables

    lzero = .true.; if (present(zero_first)) lzero = zero_first    
    which = 2;      if(present(flags)) which = flags
    
    if( which /= 2 ) &
      stop "pf_apply_mat_backward can only be used for restricting the backward integrals with which==2"

    n = size(mat, dim=1)
    m = size(mat, dim=2)
        
    do i = 1, n
      if (lzero) call dst(n+1-i)%setval(0.0_pfdp, 2)
      do j = 1, m
        if (abs(a*mat(i, j)) /= 0.0_pfdp)  call dst(n+1-i)%axpy(a * mat(i, j), src(m+1-j), 2)
      end do
    end do
  end subroutine pf_apply_mat_backward
  
end module pf_mod_utils
