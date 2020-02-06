!! Useful subroutines that don't  fit in other modules
!
! This file is part of LIBPFASST.

!> Module with useful subroutines that don't  fit in other modules
module pf_mod_utils
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_stop
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
    
    if (pf%save_timings > 1) call pf_start_timer(pf, T_RESIDUAL,level_index)

    lev => pf%levels(level_index)
    call lev%ulevel%sweeper%residual(pf,level_index, dt, flag)

    ! compute max residual norm
    !   sol_norms(1) = lev%Q(1)%norm(flag) ! for adjoint
     sol_norms = lev%Q(1)%norm(flag) ! for adjoint    
!    do m = 1, lev%nnodes-1
!       res_norms(m) = lev%R(m)%norm(flag)
       res_norms = lev%R(lev%nnodes-1)%norm(flag)
       !       sol_norms(m+1) = lev%Q(m+1)%norm(flag) ! only the value at lev%nnodes is needed for forward integration, right?
!       sol_norms(m+1) = sol_norms(1) ! only the value at lev%nnodes is needed for forward integration, right?        
!    end do
    
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

    call pf_set_resid(pf,lev%index,lev%residual)
    
    if (pf%save_timings > 1) call pf_stop_timer(pf, T_RESIDUAL,level_index)


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
       end do
    end if

    !> subtract out the solution value
    if (present(flags)) then
      do m = 1, lev%nnodes-1      
        if( (flags .eq. 0) .or. (flags .eq. 1) ) then
          call lev%R(m)%copy(lev%I(m), 1)
          !          call lev%R(m)%axpy(1.0_pfdp, lev%Q(1), 1)
          call lev%R(m)%axpy(1.0_pfdp, lev%q0, 1)          
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
        call lev%R(m)%axpy(1.0_pfdp, lev%q0)
        call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
      end do
    end if

  end subroutine pf_generic_residual

  !>  Output the current residual in the solution
  subroutine pf_echo_residual(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index
    if( pf%state%sweep < 0) print *, pf%state%sweep
    if( pf%state%iter < 0) print *, pf%state%iter
    print '("resid: time: ", f10.4," step: ",i8.8," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
         pf%state%t0+pf%state%dt,pf%state%step+1, pf%rank, pf%state%iter,level_index,pf%levels(level_index)%residual    
    
    call flush(6)
  end subroutine pf_echo_residual

  !>  Subroutine to store a residual value
  subroutine pf_set_resid(pf,level_index,resid)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index
    real(pfdp), intent(in) :: resid
    if( pf%state%sweep < 0) print *, pf%state%sweep
    if( pf%state%pfblock < 0) print *, pf%state%pfblock
    if (pf%save_residuals .and. pf%state%iter>-1)  then
       pf%results%residuals(level_index, pf%state%pfblock,pf%state%iter+1, pf%state%sweep) = resid
    end if
    
  end subroutine pf_set_resid
  !>  Subroutine to store a delta_q0 value
  subroutine pf_set_delta_q0(pf,level_index,delta)
    type(pf_pfasst_t), intent(inout)           :: pf
    integer, intent(in) :: level_index
    real(pfdp), intent(in) :: delta

    if( pf%state%sweep < 0) print *, pf%state%sweep    
    if( pf%state%pfblock < 0) print *, pf%state%pfblock
    if (pf%save_delta_q0 .and. pf%state%iter>-1)  then
       pf%results%delta_q0(level_index, pf%state%pfblock,pf%state%iter+1, pf%state%sweep) = delta
    end if
    
  end subroutine pf_set_delta_q0
  
  !>  Subroutine to store an error value
  subroutine pf_set_error(pf,level_index,error)
    type(pf_pfasst_t), intent(inout)           :: pf
    integer, intent(in) :: level_index
    real(pfdp), intent(in) :: error
    if( pf%state%sweep < 0) print*, pf%state%sweep    
    if (pf%state%iter+1 < 1)  return
    if (pf%save_errors)  then
       pf%results%errors(level_index, pf%state%pfblock,pf%state%iter+1, pf%state%sweep) = error
    end if
    
  end subroutine pf_set_error

  !>  Subroutine to set the final the iteration number for convergence
  subroutine pf_set_iter(pf,iter)
    type(pf_pfasst_t), intent(inout)           :: pf
    integer, intent(in) :: iter

    pf%results%iters(pf%state%pfblock) = iter
    
  end subroutine pf_set_iter
  

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

  subroutine pf_delta_q0(pf,level_index)
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to sweep

    class(pf_level_t), pointer :: lev  !!  Level on which to spread
    lev => pf%levels(level_index)   !!  Assign level pointer
    call lev%delta_q0%axpy(-1.0_pfdp,lev%q0)
    lev%max_delta_q0=lev%delta_q0%norm()
    call pf_set_delta_q0(pf,level_index,lev%max_delta_q0)

  end subroutine pf_delta_q0
  
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
        
    n = SIZE(mat, dim=1)
    m = SIZE(mat, dim=2)
    
    do i = 1, n
      if (lzero) call dst(i)%setval(0.0_pfdp, flags)
      do j = 1, m
         if (abs(a*mat(i, j)) /= 0.0_pfdp) then
            call dst(i)%axpy(a * mat(i, j), src(j), flags)
         end if
         
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
         call pf_stop(__FILE__,__LINE__,'pf_apply_mat_backward can only be used for restricting the backward integrals with which==2')

    n = SIZE(mat, dim=1)
    m = SIZE(mat, dim=2)
        
    do i = 1, n
      if (lzero) call dst(n+1-i)%setval(0.0_pfdp, 2)
      do j = 1, m
        if (abs(a*mat(i, j)) /= 0.0_pfdp)  call dst(n+1-i)%axpy(a * mat(i, j), src(m+1-j), 2)
      end do
    end do
  end subroutine pf_apply_mat_backward
  
  function convert_logical(q) result(q_string)
    logical,intent(in) :: q  ! true or false
    character(len=5)::  q_string
    if (q) then
       q_string=' true'
    else
       q_string='false'
    end if
  end function convert_logical

  function convert_int_array(q,n) result(q_string)
    integer,intent(in) :: n     ! length of array
    integer,intent(in) :: q(n)  ! integer array
    character(len=15)::  q_string

    character(len=15)::  f_string  !  format string
    integer i
    write(f_string,"(*(G0,:,','))") q
    q_string=adjustr('['//trim(f_string)//']')
  end function convert_int_array
  function wrap_timer_name(tname) result(q_string)
    character(len=10)::  tname
    character(len=12)::  q_string

    integer i

    q_string=adjustr('"'//trim(tname)//'"')

  end function wrap_timer_name
  

  function convert_real_array(q,n) result(q_string)
    integer,intent(in) :: n     ! length of array
    real(pfdp) , intent(in) :: q(n)  ! real array
    character(len=128)::  q_string

    character(len=128)::  f_string  !  format string
    integer i
    write(f_string,"(*(e15.6,:,','))") q
    q_string=adjustl('['//trim(f_string)//']')

  end function convert_real_array
  
  
end module pf_mod_utils
