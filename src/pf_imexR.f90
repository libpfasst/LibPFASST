!
!>  Module of the  the derived sweeper class for doing IMEX sweeps for an equation of the form
!!            y' = f_1(y) + f_2(y)
!!  The f_1 piece is treated explicitly and f_2 implicitly
!!  Afer this sweeper is initialized (usually in main), the locgical flags can be changed if desired
!!
!!     explicit:  Make false if there is no explicit piece
!!
!!     implicit:  Make false if there is no implicit piece
!!
!!     use_LUqt:  Make false if backward Euler sweepers should be used instead of the 'LU trick'
!!
!!  The user needs to supply the feval and fcomp routines for a given example
!!
!!  This version uses the residual to form SDC correction terms


module pf_mod_imexR
  use pf_mod_dtype
  use pf_mod_utils

  implicit none

  !>  Define sweeper type
  type, extends(pf_sweeper_t), abstract :: pf_imexR_t
     real(pfdp), allocatable :: QtilE(:,:)   !<  Approximate explicit quadrature rule
     real(pfdp), allocatable :: QtilI(:,:)   !<  Approximate implicit quadrature rule
     real(pfdp), allocatable :: dtsdc(:)     !<  SDC step sizes
     real(pfdp), allocatable :: tsdc(:)      !<  Time at the nodes

     logical                 :: explicit = .true. !<  True if there is an explicit piece
     logical                 :: implicit = .true. !<  True if there an implicit piece
     logical                 :: use_LUq = .true.  !<  Use the LU trick if true

   contains
     procedure(pf_f_eval_p), deferred :: f_eval   !<  RHS function evaluations
     procedure(pf_f_comp_p), deferred :: f_comp   !<  Implicit solver
     !>  Set the generic functions
     procedure :: sweep      => imexR_sweep
     procedure :: initialize => imexR_initialize
     procedure :: evaluate   => imexR_evaluate
     procedure :: integrate  => imexR_integrate
     procedure :: residual   => imexR_residual
     procedure :: spreadq0   => imexR_spreadq0
     procedure :: evaluate_all => imexR_evaluate_all
     procedure :: destroy   => imexR_destroy
     procedure :: imexR_destroy
  end type pf_imexR_t

  interface
     !>  This is the interface for the routine to compute the RHS function values
     !>  Evaluate f_piece(y), where piece is one or two 
     subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
       import pf_imexR_t, pf_encap_t, pfdp
       class(pf_imexR_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y        !<  Argument for evaluation
       real(pfdp),        intent(in   ) :: t        !<  Time at evaluation
       integer,    intent(in   ) :: level_index     !<  Level index
       class(pf_encap_t), intent(inout) :: f        !<  RHS function value
       integer,    intent(in   ) :: piece           !<  Which piece to evaluate
     end subroutine pf_f_eval_p
     !>  Solve the equation y - dtq*f_2(y) =rhs
     subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       import pf_imexR_t, pf_encap_t, pfdp
       class(pf_imexR_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y      !<  Solution of implicit solve 
       real(pfdp),        intent(in   ) :: t      !<  Time of solve
       real(pfdp),        intent(in   ) :: dtq    !<  dt*quadrature weight
       class(pf_encap_t), intent(in   ) :: rhs    !<  RHS for solve
       integer,    intent(in   ) :: level_index   !<  Level index
       class(pf_encap_t), intent(inout) :: f      !<  f_2 of solution y
       integer,    intent(in   ) :: piece         !<  Which piece to evaluate
     end subroutine pf_f_comp_p
  end interface

contains

  !> Perform nsweep SDC sweeps on level Lev and set qend appropriately.
  subroutine imexR_sweep(this, pf, level_index, t0, dt, nsweeps, flags)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_imexR_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf      !<  PFASST structure
    real(pfdp),        intent(in   ) :: t0             !<  Time at beginning of time step
    real(pfdp),        intent(in   ) :: dt             !<  time step size
    integer,             intent(in)    :: level_index  !<  which level this is
    integer,             intent(in)    :: nsweeps      !<  number of sweeps to do
    integer, optional, intent(in   ) :: flags    


    class(pf_level_t), pointer :: lev    !<  points to current level

    integer     :: m, n,k   !<  Loop variables

    lev => pf%levels(level_index)   !<  Assign level pointer
    this%tsdc = t0+dt*lev%nodes

    call start_timer(pf, TLEVEL+lev%index-1)
    ! compute residual
    call pf_residual(pf, lev, dt)

    do k = 1,nsweeps   !>  Loop over sweeps
       call call_hooks(pf, level_index, PF_PRE_SWEEP)    

       ! Store all function values
       call imexR_save(lev)

       !  Recompute the first function value if this is first sweep
       if (k .eq. 1) then
          call lev%Q(1)%copy(lev%q0)
          call imexR_evaluate(this,lev,t0,1)
       end if

       ! Assign the old function value the difference in function values
       m = 1
       if (this%explicit) &
            call lev%pF(m,1)%axpy(-1.0_pfdp,lev%F(m,1))
       if (this%implicit) &
            call lev%pF(m,2)%axpy(-1.0_pfdp,lev%F(m,2))
       
       ! do the sub-stepping in sweep
       do m = 1, lev%nnodes-1
          !>  Accumulate rhs in Residual 
          do n = 1, m
             if (this%explicit) &
                  call lev%R(m)%axpy(-dt*this%QtilE(m,n), lev%pF(n,1))
             if (this%implicit) &
                  call lev%R(m)%axpy(-dt*this%QtilI(m,n), lev%pF(n,2))
          end do

          !  Add the starting old sol term
          call lev%R(m)%axpy(1.0_pfdp, lev%Q(m+1))

          !  Solve for the implicit piece
          if (this%implicit) then
             !  subtract old F
             call lev%R(m)%axpy(-dt*this%QtilI(m,m+1), lev%pF(m+1,2))
             call this%f_comp(lev%Q(m+1), this%tsdc(m+1), dt*this%QtilI(m,m+1), lev%R(m), lev%index,lev%F(m+1,2),2)
          else
             call lev%Q(m+1)%copy(lev%R(m))
          end if
          !  Compute explicit function on new value
          if (this%explicit) &
               call this%f_eval(lev%Q(m+1), this%tsdc(m+1), lev%index, lev%F(m+1,1),1)
          !  Make next pF the difference
          if (this%explicit) &
               call lev%pF(m+1,1)%axpy(-1.0_pfdp,lev%F(m+1,1))
          if (this%implicit) &
               call lev%pF(m+1,2)%axpy(-1.0_pfdp,lev%F(m+1,2))
          
       end do
       call pf_residual(pf, lev, dt)


       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do  !>  End loop on sweeps

    call lev%qend%copy(lev%Q(lev%nnodes))
    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine imexR_sweep

  !> Subroutine to initialize matrices and space for sweeper
  subroutine imexR_initialize(this, lev)
    class(pf_imexR_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev    !<  Current level

    integer    :: m,n, nnodes

    this%npieces = 2

    nnodes = lev%nnodes
    allocate(this%QtilE(nnodes-1,nnodes))  
    allocate(this%QtilI(nnodes-1,nnodes))  
    allocate(this%dtsdc(nnodes-1))  
    allocate(this%tsdc(nnodes))  

    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    !>  Array of substep sizes
    this%dtsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)

    ! Implicit matrix
    if (this%use_LUq) then 
       ! Get the LU
       call myLUq(lev%qmat,lev%LUmat,lev%nnodes,0)
       this%QtilI = lev%LUmat
    else 
       this%QtilI =  lev%qmatBE
    end if

    ! Explicit matrix
    this%QtilE =  lev%qmatFE
  end subroutine imexR_initialize

  !>  Subroutine to deallocate sweeper
  subroutine imexR_destroy(this, lev)
    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev   !<  Current level

    deallocate(this%QtilE)
    deallocate(this%QtilI)
    deallocate(this%dtsdc)
    deallocate(this%tsdc)

  end subroutine imexR_destroy


  !> Subroutine to compute  Picard integral of function values
  subroutine imexR_integrate(this, lev, qSDC, fSDC, dt, fintSDC, flags)
    class(pf_imexR_t), intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev          !<  Current level
    class(pf_encap_t), intent(in   ) :: qSDC(:)      !<  Solution values
    class(pf_encap_t), intent(in   ) :: fSDC(:, :)   !<  RHS Function values
    real(pfdp),        intent(in   ) :: dt           !<  Time step
    class(pf_encap_t), intent(inout) :: fintSDC(:)   !<  Integral from t_n to t_m
    integer, optional, intent(in   ) :: flags    

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          if (this%explicit) &
               call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,1))
          if (this%implicit) &
               call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,2))
       end do
    end do
  end subroutine imexR_integrate

  !> Subroutine to compute  Residual
  subroutine imexR_residual(this, lev, dt, flags)
    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev  !<  Current level
    real(pfdp),        intent(in   ) :: dt   !<  Time step
    integer, optional, intent(in   ) :: flags    
    call pf_generic_residual(this, lev, dt)
  end subroutine imexR_residual

  subroutine imexR_spreadq0(this, lev, t0, flags, step)
    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t0
    integer, optional, intent(in   ) :: flags    
    integer, optional, intent(in   ) :: step    
    call pf_generic_spreadq0(this, lev, t0)
  end subroutine imexR_spreadq0

  !> Subroutine to evaluate function value at node m
  subroutine imexR_evaluate(this, lev, t, m, flags, step)

    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev  !<  Current level
    real(pfdp),        intent(in   ) :: t    !<  Time at which to evaluate
    integer,           intent(in   ) :: m    !<  Node at which to evaluate
    integer, optional, intent(in   ) :: flags    
    integer, optional, intent(in   ) :: step    

    if (this%explicit) &
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1),1)
    if (this%implicit) &
         call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2),2)
  end subroutine imexR_evaluate

  !> Subroutine to evaluate the function values at all nodes
  subroutine imexR_evaluate_all(this, lev, t, flags, step)
    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev   !<  Current level
    real(pfdp),        intent(in   ) :: t(:)  !<  Array of times at each node
    integer, optional, intent(in   ) :: flags    
    integer, optional, intent(in   ) :: step    
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine imexR_evaluate_all

  !>  Save function values so that difference can be computed
  subroutine imexR_save(lev)
    class(pf_level_t), intent(inout) :: lev  !<  Level to save on
    
    integer :: m, p
    
    do m = 1, lev%nnodes
       do p = 1,size(lev%F(1,:))
          call lev%pF(m,p)%copy(lev%F(m,p))
       end do
    end do
    
  end subroutine imexR_save
  
end module pf_mod_imexR
