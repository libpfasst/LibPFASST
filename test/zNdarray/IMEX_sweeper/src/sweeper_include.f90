  !>  Helper function to return sweeper pointer
  function as_my_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(my_sweeper_t), pointer :: r
    select type(sweeper)
    type is (my_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_my_sweeper

  !>  Routine to set up sweeper variables and operators
  subroutine initialize(this, pf,level_index)
    use probin, only:  dom_size,Ndim
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer, intent(in) :: level_index
    
    class(pf_level_t), pointer :: lev
    
    lev => pf%levels(level_index)
    this%nx=lev%lev_shape(1)
    
    ! call superclass initialize
    call this%imex_initialize(pf,level_index)

    ! Set IMEX sweeper to have both implicit and explicit parts
    this%implicit=.TRUE.
    this%explicit=.TRUE.
    
    ! Allocate fft 
    allocate(this%fft_tool)
    call this%fft_tool%fft_setup(lev%lev_shape,Ndim,dom_size)
    ! Make differentiation matrices
    allocate(this%fft_ops)
    call this%fft_ops%init(this%fft_tool,this%nx)
    !  Call routine to allocate the temp variable
    call this%initialize_tmp()
  end subroutine initialize
  
  
  ! Call superclass destroy and deallocate all array type parameters
  subroutine destroy(this,pf,level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index
    
    call this%imex_destroy(pf,level_index)
    
    ! deallocate local stuff
    deallocate(this%tmp)
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)
    call this%fft_ops%destroy()
    deallocate(this%fft_ops)
    
  end subroutine destroy
  
  
  ! Evaluate the rhs function f(t,y(t)) at y, t.
  subroutine f_eval(this, y, t, level_index, f,piece)
    ! arguments
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in)    :: y
    real(pfdp),          intent(in)    :: t
    integer,             intent(in)    :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece  ! explicit or implicit
    
    !  Cast the generic encaps as zndarrays
    this%y_encap => cast_as_zndarray(y) 
    this%f_encap => cast_as_zndarray(f) 

    !  Get the raw arrays from the zndarray structure
    call this%y_encap%get_array(this%yvec)        
    call this%f_encap%get_array(this%fvec)        

    !  Compute the function value
    if (piece .eq. 1) then  !  Explicit
       call f_NL(this%yvec,this%fvec,this%fft_ops%opNL,this%tmp,this%fft_tool)
    else  !  Implicit
       this%fvec=this%fft_ops%opL*this%yvec
    endif
    
  end subroutine f_eval
  
  ! Solve for (I-dtq*L)y=rhs and return f=Ly also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    
    if (piece == 2) then
       !  Cast the generic encaps as zndarrays
       this%y_encap => cast_as_zndarray(y) 
       this%f_encap => cast_as_zndarray(f) 
       this%rhs_encap => cast_as_zndarray(rhs) 

       !  Get the raw arrays from the zndarray structure
       call this%y_encap%get_array(this%yvec)
       call this%f_encap%get_array(this%fvec)
       call this%rhs_encap%get_array(this%rhsvec)
       
       ! Apply the inverse operator in spectral space
       this%yvec=this%rhsvec/(1.0_pfdp - dtq*this%fft_ops%opL)
       
       ! Compute the function value
       this%fvec = (this%yvec - this%rhsvec)/dtq
    else
       stop  'Bad value for piece in f_comp'
    end if
  end subroutine f_comp
