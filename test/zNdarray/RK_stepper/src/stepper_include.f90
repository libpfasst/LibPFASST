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
    use probin, only:  rk_order,nsteps_rk,dom_size,Ndim
    class(my_stepper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer, intent(in) :: level_index
    
    integer ::  nx,nnodes
    class(pf_level_t), pointer :: lev
    
    lev => pf%levels(level_index)
    this%nx=lev%lev_shape(1)
    nx=this%nx
    nnodes=lev%nnodes
    
    !>  Set variables for explicit and implicit parts
    this%implicit=.TRUE.
    this%explicit=.TRUE.

    ! call superclass initialize
    this%order=rk_order(level_index)
    call this%ark_initialize(pf,level_index)
    this%nsteps=nsteps_rk(level_index)
    
    ! allocate fft & differentiation matrices
    allocate(this%fft_tool)
    call this%fft_tool%fft_setup(lev%lev_shape,Ndim,dom_size)
    allocate(this%fft_ops)
    call this%fft_ops%init(this%fft_tool,this%nx)

    !  Call routine to allocate local storage
    call this%initialize_tmp()

  end subroutine initialize
  
  
  ! DESTROY: calls superclass destroy and deallocates all array type parameters
  subroutine destroy(this,pf,level_index)
    class(my_stepper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index
    
    call this%ark_destroy(pf,level_index)
    
    ! deallocate arrays
    deallocate(this%tmp)
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)
    call this%fft_ops%destroy()
    deallocate(this%fft_ops)
    
  end subroutine destroy

  !F_EVAL: evaluates the nonlinear function N(t,y(t)) at y, t.
  subroutine f_eval(this, y, t, level_index, f,piece)
    use probin, only: eq_type
    ! arguments
    class(my_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in)    :: y
    real(pfdp),          intent(in)    :: t
    integer,             intent(in)    :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    
    this%y_encap => cast_as_zndarray(y) 
    this%f_encap => cast_as_zndarray(f) 
    call this%y_encap%get_array(this%yvec)        
    call this%f_encap%get_array(this%fvec)        
    if (piece .eq. 1) then  !  Explicit
       call f_NL(this%yvec,this%fvec,this%fft_ops%opNL,this%tmp,this%fft_tool)
    else  !  Implicit
       this%fvec=this%fft_ops%opL*this%yvec
    endif

  end subroutine f_eval

  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    class(my_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    
    this%y_encap => cast_as_zndarray(y) 
    this%f_encap => cast_as_zndarray(f) 
    this%rhs_encap => cast_as_zndarray(rhs) 
    
    if (piece == 2) then
       
       call this%y_encap%get_array(this%yvec)
       call this%f_encap%get_array(this%fvec)
       call this%rhs_encap%get_array(this%rhsvec)
       
       ! Apply the inverse operator in spectral space
       this%yvec=this%rhsvec/(1.0_pfdp - dtq*this%fft_ops%opL)
       
       this%fvec = (this%yvec - this%rhsvec) / dtq
    else
       stop 'Bad value for piece in f_comp'
    end if
  end subroutine f_comp
  
