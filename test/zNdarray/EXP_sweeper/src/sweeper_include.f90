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
    use probin, only:  Ndim,dom_size
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer, intent(in) :: level_index

    class(pf_level_t), pointer :: lev
    
    lev => pf%levels(level_index)
    this%nx=lev%lev_shape(1)
    this%nnodes=lev%nnodes
    
    ! call superclass initialize
    call this%exp_initialize(pf,level_index)
    
    
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
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index
    
    call this%exp_destroy(pf,level_index)
    
    ! deallocate arrays
    deallocate(this%tmp)
    deallocate(this%P_sweep)
    deallocate(this%W_sweep)
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)
    call this%fft_ops%destroy()
    deallocate(this%fft_ops)
       
  end subroutine destroy

  ! ==================================================================
  ! Required subroutines for any sweeper that inherits from pf_fexp_t
  ! ==================================================================

  !F_EVAL: evaluates the nonlinear function f(t,y(t)) at y, t.
  subroutine f_eval(this, y, t, level_index, f)
    use probin, only: eq_type,splitting
    ! arguments
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in)    :: y
    real(pfdp),          intent(in)    :: t
    integer,             intent(in)    :: level_index
    class(pf_encap_t),   intent(inout) :: f
    
    ! local variables
    type(pf_fft_t),     pointer :: fft
    integer :: nx,K
    
    this%y_encap => cast_as_zndarray(y) 
    this%f_encap => cast_as_zndarray(f) 
    
    call this%y_encap%get_array(this%yvec)        
    call this%f_encap%get_array(this%fvec)        
    fft => this%fft_tool

    !  Call the nonlinear function evaluation  (in utils)
    call f_NL(this%yvec,this%fvec,this%fft_ops%opNL,this%tmp,fft)
  
  end subroutine f_eval
  
  ! =================================================================================
  ! EXPSWEEPSUBSTEP: Computes the jth substep of the exponential sweeper
  !
  !     		exp(x L) y_{n,j}^{k+1} + h \varphi_1(x L) [N_{n,j}^{[k+1]} - N_{n,j}^{[k]}] + I_{n,j}
  !
  !	where x = dt * (t_{n,j+1} - t_{n,j}), I_{j} are the exponential integrals
  !
  !		I_{n,j} = x \int_{0}^{1} exp(x L(1-s)) Q_j(s) ds
  !
  !	and Q_j(s) is a Lagrange interpolating polynomial that satisfies
  !
  !		Q_j(\tau_{j,i}) = N(y(t_{n,i})) for i = 1 ... m
  !
  !	for \tau_{j,i} = (t_{n,i} - t_{n,j}) / (t{n,j+1} - t_{n,j})
  !
  !	NOTE: The operator L is not passed in as a parameter, and must be
  !       implemented appropriately within the class.
  !
  ! Arguments
  !
  !   j   	(input) Integer
  !       	substep index for determining t: t = t_{n,j+1} - t_{n,j}
  !
  !   y_jp1 (inout) pf_encap_t
  !		      stores the solution y_{n,j+1}^{k+1}
  !
  !   dt  	(input) real
  !       	stepsize
  !
  !   y_j  (input) pf_encap_t
  !		     stores the solution y_{n,j}^{k+1}
  !
  !   F    (input) pf_encap_t(:,:)
  !		     should be the lev%F matrix of dimension (nnodes-1 x 2). First component
  !        stores the nonlinear term, and the second component is used by this
  !        function to store the exponential integral terms I_{n,j}
  !
  !   Nk  	(input) pf_encap_t(:)
  !		stores N(y^{k+1}_j) for j = 1 ... m
  ! =================================================================================
  subroutine expSweepSubstep(this, y_jp1, j, dt, y_j, F, Nk)

    ! Arguments
    class(my_sweeper_t),  intent(inout) :: this
    class(pf_encap_t),  	intent(inout) :: y_jp1
    integer,            	intent(in)    :: j
    real(pfdp),         	intent(in)    :: dt
    class(pf_encap_t),  	intent(in)    :: y_j
    class(pf_encap_t),  	intent(inout)    :: F(:,:)
    class(pf_encap_t),  	intent(inout)    :: Nk(:)
    ! Variables
    integer :: i

    if(dt .ne. this%dt_sweep) then
       call this%initWSweep(dt)
       this%dt_sweep = dt
    end if

    ! Start with zero
    call y_jp1%setval(0.0_pfdp)

    ! add phi_0(h L) y_{n,j}
    call this%P_axpy(y_jp1,y_j,j,1)        

    ! add h * phi_1(h L) * (N_j^{[k+1]} - N_j^{[k]})
    call this%P_daxpy(y_jp1,F(j,1),Nk(j),j,2)
 
    ! add contribution from exponential integral I_{n,j}  and store in F(j,2)
    call F(j, 2)%setval(0.0_pfdp)
    do i = 1, size(Nk)
       call this%W_axpy(F(j,2),Nk(i),i,j)       
    end do

    ! add to new y
    call y_jp1%axpy(1.0_pfdp,F(j,2))
  end subroutine expSweepSubstep

  subroutine expResidualSubstep(this, y_jp1, j, dt, y_j, F)

    ! Arguments
    class(my_sweeper_t),  intent(inout) :: this
    class(pf_encap_t),  	intent(inout) :: y_jp1
    integer,            	intent(in)    :: j
    real(pfdp),         	intent(in)    :: dt
    class(pf_encap_t),  	intent(in)    :: y_j
    class(pf_encap_t),  	intent(in)    :: F(:,:)
    ! Variables
    integer :: i

    if(dt .ne. this%dt_sweep) then
      call this%initWSweep(dt)
      this%dt_sweep = dt
    end if

    ! Start with zero for next solution
    call y_jp1%setval(0.0_pfdp)

    ! add phi_0(h L) y_{n,j}
    call this%P_axpy(y_jp1,y_j,j,1)        

    ! add contribution from exponential integral I_{n,j} 
    do i = 1, size(F,1)
       call this%W_axpy(y_jp1,F(i,1),i,j)        
    end do

  end subroutine expResidualSubstep

