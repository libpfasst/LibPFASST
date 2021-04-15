! File with stepper routines that are DIM independent

!>  Helper function to return stepper pointer
function as_my_stepper(stepper) result(r)
  class(pf_stepper_t), intent(inout), target :: stepper
  class(my_stepper_t), pointer :: r
  select type(stepper)
  type is (my_stepper_t)
     r => stepper
  class default
     stop
  end select
end function as_my_stepper

!>  Routine to set up stepper variables and operators
subroutine initialize(this, pf,level_index)
  use probin, only:  Ndim,dom_size,rk_order,nsteps_rk,dt
  class(my_stepper_t), intent(inout) :: this
  type(pf_pfasst_t),   intent(inout),target :: pf
  integer, intent(in) :: level_index
  
  integer ::  nx,nnodes
  class(pf_level_t), pointer :: lev
  real(pfdp) :: dt_RK
  lev => pf%levels(level_index)
  this%nx=lev%lev_shape(1)
  nx=this%nx
  nnodes=lev%nnodes
  
  !  set some info for the superclass
  this%order = rk_order(level_index)
  this%nsteps = nsteps_rk(level_index)

  ! call superclass initialize
  call this%erk_initialize(pf,level_index)

  ! allocate fft & differentiation matrices
  allocate(this%fft_tool)
  call this%fft_tool%fft_setup(lev%lev_shape,Ndim,dom_size)
  allocate(this%fft_ops)
  call this%fft_ops%init(this%fft_tool,this%nx)
  
  !  Call routine to allocate local storage
  dt_RK=dt/real(nsteps_rk(level_index),pfdp)
  call this%initialize_tmp()

  !  Compute original RK coeff
  call this%initRKCoeff(this%fft_ops%opL, dt_RK)
  this%dt_phi=dt_RK
  
end subroutine initialize

! DESTROY: calls superclass destroy and deallocates all array type parameters
subroutine destroy(this,pf,level_index)
  class(my_stepper_t), intent(inout) :: this
  type(pf_pfasst_t),   intent(inout),target :: pf
  integer,             intent(in)    :: level_index
  
  call this%erk_destroy(pf,level_index)
  deallocate(this%A_phi)
  deallocate(this%b_phi)
  deallocate(this%d_phi)
  deallocate(this%tmp)
  call this%fft_tool%fft_destroy()
  deallocate(this%fft_tool)
  call this%fft_ops%destroy()
  deallocate(this%fft_ops)
  
end subroutine destroy

!F_EVAL: evaluates the nonlinear function f(t,y(t)) at y, t.
subroutine f_eval(this, y, t, level_index, f)
  use probin, only: eq_type,splitting
  ! arguments
  class(my_stepper_t), intent(inout) :: this
  class(pf_encap_t),   intent(in)    :: y
  real(pfdp),          intent(in)    :: t
  integer,             intent(in)    :: level_index
  class(pf_encap_t),   intent(inout) :: f
  
  ! local variables
  type(pf_fft_t),     pointer :: fft
  
  this%y_encap => cast_as_zndarray(y) 
  this%f_encap => cast_as_zndarray(f) 
  call this%y_encap%get_array(this%yvec)        
  call this%f_encap%get_array(this%fvec)        
  fft => this%fft_tool
  
  call f_NL(this%yvec,this%fvec,this%fft_ops%opNL,this%tmp,fft)
  
end subroutine f_eval

  subroutine compA(this, dt, i, j, F, val)
    
    ! Arguments
    class(my_stepper_t), intent(inout) :: this
    real(pfdp), intent(in) :: dt
    integer, intent(in) :: i, j
    class(pf_encap_t), intent(in)    :: F(:,:)
    class(pf_encap_t), intent(inout) :: val
    
    ! Variables
    integer :: coeff_index
    
    if(abs(this%dt_phi-dt) .gt. 1.0d-14) then
       call this%initRKCoeff(this%fft_ops%opL, dt)
       this%dt_phi = dt
    endif
    
    coeff_index = this%AF(j,i)
    this%f_encap => cast_as_zndarray(F(j,1))
    this%y_encap => cast_as_zndarray(val)
    call this%f_encap%get_array(this%fvec)
    call this%y_encap%get_array(this%yvec)

     call this%mult_Aphi(coeff_index,val,F(j,1))            
  end subroutine compA
  
  subroutine compB(this, dt, i, F, val)
    
    ! Arguments
    class(my_stepper_t), intent(inout) :: this
    real(pfdp), intent(in) :: dt
    integer, intent(in) :: i
    class(pf_encap_t), intent(in) :: F(:,:)
    class(pf_encap_t), intent(inout) :: val

    if(abs(this%dt_phi-dt) .gt. 1.0d-14) then
       call this%initRKCoeff(this%fft_ops%opL, dt)
       this%dt_phi = dt
    endif
    call this%mult_bphi(i,val,F(i,1))    
  end subroutine compB
  
  subroutine compD(this, dt, i, y_n, val)
    
    ! Arguments
    class(my_stepper_t), intent(inout) :: this
    real(pfdp), intent(in) :: dt
    integer, intent(in) :: i
    class(pf_encap_t), intent(in) :: y_n
    class(pf_encap_t), intent(inout) :: val
    
    if(abs(this%dt_phi-dt) .gt. 1.0d-14) then
       call this%initRKCoeff(this%fft_ops%opL, dt)
       this%dt_phi = dt 
    endif

    this%f_encap => cast_as_zndarray(val)
    this%y_encap => cast_as_zndarray(y_n)
    call this%f_encap%get_array(this%fvec)
    call this%y_encap%get_array(this%yvec)

    call this%mult_dphi(i,val,y_n)
  end subroutine compD
