!! Verlet type sweeoer for 2nd order problems
!
!  This will work for
! $$   q'=p, p'=f(q) $$ 
!
module pf_mod_verlet
  use pf_mod_dtype
  use pf_mod_utils

  implicit none

  !> Define sweeper type
  type, extends(pf_sweeper_t), abstract :: pf_verlet_t
     integer :: whichQQ
     integer :: doLU
     real(pfdp) :: Htol, H0

     !  Matrices
     real(pfdp), ALLOCATABLE :: QQmat(:,:)
     real(pfdp), ALLOCATABLE :: Qmat(:,:)     
     real(pfdp), ALLOCATABLE :: QQtil(:,:)
     real(pfdp), ALLOCATABLE :: Qtil(:,:)     
     real(pfdp), ALLOCATABLE :: QQLU(:,:)
     real(pfdp), ALLOCATABLE :: DQQmat(:,:)
     real(pfdp), ALLOCATABLE :: DQQtil(:,:)
     real(pfdp), ALLOCATABLE :: DQtil(:,:)     
     real(pfdp), ALLOCATABLE :: DQQLU(:,:)
     real(pfdp), ALLOCATABLE :: bvec(:)
     real(pfdp), ALLOCATABLE :: bbarvec(:)
     real(pfdp), allocatable :: dtsdc(:)     !!  SDC step sizes
     real(pfdp), allocatable :: tsdc(:)     !!  SDC step sizes     
     logical :: iqend  !  Decide whether to set qend by another Picard

     class(pf_encap_t), allocatable :: rhs   !! holds rhs for implicit solve     
   contains
     procedure(pf_f_eval_p), deferred :: f_eval        !!  RHS function evaluations
     procedure(pf_f_comp_p), deferred :: f_comp        !!  Implicit solver
     procedure(pf_f_eval_p), deferred :: hamiltonian   !!  Hamiltonian
     !>  Set the generic functions
     procedure :: sweep      => verlet_sweep
     procedure :: initialize => verlet_initialize
     procedure :: evaluate   => verlet_evaluate
     procedure :: integrate  => verlet_integrate
     procedure :: residual   => verlet_residual
     procedure :: spreadq0   => verlet_spreadq0
     procedure :: evaluate_all => verlet_evaluate_all
     procedure :: destroy   => verlet_destroy
     procedure :: verlet_destroy
  end type pf_verlet_t

  interface
     !>  This is the interface for the routine to compute the RHS function values
     !>  Evaluate f_piece(y), where piece is one or two 
     subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
       !>  Evaluate f_piece(y), where piece is one or two 
       import pf_verlet_t, pf_encap_t, pfdp
       class(pf_verlet_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y        !!  Argument for evaluation
       real(pfdp),        intent(in   ) :: t        !!  Time at evaluation
       integer,    intent(in   ) :: level_index     !!  Level index
       class(pf_encap_t), intent(inout) :: f        !!  RHS function value
       integer,    intent(in   ) :: piece           !!  Which piece to evaluate
     end subroutine pf_f_eval_p
     !>  Solve the equation y - dtq*f_2(y) =rhs
     subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       import pf_verlet_t, pf_encap_t, pfdp
       class(pf_verlet_t),  intent(inout) :: this


       class(pf_encap_t), intent(inout) :: y      !!  Solution of implicit solve 
       real(pfdp),        intent(in   ) :: t      !!  Time of solve
       real(pfdp),        intent(in   ) :: dtq    !!  dt*quadrature weight
       class(pf_encap_t), intent(in   ) :: rhs    !!  RHS for solve
       integer,    intent(in   ) :: level_index   !!  Level index
       class(pf_encap_t), intent(inout) :: f      !!  f_2 of solution y
       integer,    intent(in   ) :: piece         !!  Which piece to evaluate
     end subroutine pf_f_comp_p
  end interface
contains

  !-----------------------------------------------------------------------------
  ! Perform one SDC sweep on level lev_index and set qend appropriately
  !
  subroutine verlet_sweep(this, pf, level_index, t0, dt,nsweeps, flags)  
    use pf_mod_timer
    use pf_mod_hooks

    !>  Inputs
    class(pf_verlet_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  which level to sweep on
    real(pfdp),        intent(in   ) :: t0           !!  Time at beginning of time step
    real(pfdp),        intent(in   ) :: dt           !!  time step size
    integer,           intent(in)    :: nsweeps      !!  number of sweeps to do
    integer, optional, intent(in   ) :: flags    

    !>  Local variables
    class(pf_level_t), pointer :: lev    !!  points to current level
    integer     :: k,m, n,nnodes
    real(pfdp)  :: t,dtmhalf,dtsq
    real(pfdp)  :: H


    nnodes = Lev%nnodes
    lev => pf%levels(level_index)   !!  Assign level pointer


    call start_timer(pf, TLEVEL+lev%index-1)

    !
    ! check hamiltonian
    !
!    call this%hamiltonian(t0+dt, Lev%qend, encapctx%m,H)
!    print *,'Ham=',H,this%Htol,this%H0


!    if ((pf%state%iter > 1) .and. (abs(H-this%H0) < this%Htol)) then
!       call Lev%encap%copy(Lev%qend, Lev%Q(nnodes))
!          print *, 'Skipping SDC sweep'
!       return
!    end if

    !
    ! compute integrals and add fas correction
    !
    dtsq = dt*dt    
    do k = 1,nsweeps

       do m = 1, lev%nnodes-1
          call lev%I(m)%setval(0.0_pfdp)
          if (pf%state%iter .eq. 1)  then  !  Do verlet on the first iteration
             do n = 1, nnodes
                call lev%I(m)%axpy(dt*this%DQtil(m+1,n), lev%F(n,1),1)
                call lev%I(m)%axpy(dtsq*this%DQQtil(m+1,n), lev%F(n,1),2)
             end do
          else
             do n = 1, nnodes
                call lev%I(m)%axpy(dt*this%DQtil(m+1,n), lev%F(n,1),1)
                call lev%I(m)%axpy(dtsq*this%DQQLU(m+1,n), lev%F(n,1),2)
             end do
          end if
          if (allocated(lev%tauQ)) then
             call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))  
          end if
       end do
       !  Recompute the first function value if this is first sweep
       if (k .eq. 1) then
          call lev%Q(1)%copy(lev%q0)
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1),1)
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2),2)
       end if

       t = t0
       ! do the sub-stepping in sweep    
       do m = 1, nnodes-1
          t = t + dt*this%dtsdc(m)
          dtmhalf = 0.5d0*dt*this%dtsdc(m)
          
          call this%rhs%setval(0.0_pfdp)
          !  Lower triangular verlet to  new piece
          if (pf%state%iter .eq. 1)  then
             do n = 1, m
                call this%rhs%axpy(dtsq*this%QQtil(m+1,n), lev%F(n,1), 2)  
             end do
          else
             do n = 1, m
                call this%rhs%axpy(dtsq*this%QQLU(m+1,n), lev%F(n,1), 2)  
             end do
          endif
          
          !>  Add the integral term
          call this%rhs%axpy(1.0_pfdp, lev%I(m),2)

          !>  Add the starting value
          call this%rhs%axpy(1.0_pfdp, lev%Q(1),2)

          !>  Add the dt*v_0
          call this%rhs%axpy(t-t0, lev%Q(1),12)

          !  Update position term (trapezoid rule)
          call lev%Q(m+1)%copy(this%rhs,2)
          
          !  update function values
          call this%f_eval(Lev%Q(m+1), t, lev%index, Lev%F(m+1,1),2)  
          
          call this%rhs%setval(0.0_pfdp)          
          !  Lower triangular verlet to  new piece
          do n = 1, m+1
             call this%rhs%axpy(dt*this%Qtil(m+1,n), Lev%F(n,1), 1)  
          end do
       
          call this%rhs%axpy(1.0_pfdp,Lev%I(m),1);
          call this%rhs%axpy(1.0_pfdp, Lev%Q(1),1) !  Start m+1 with value from m
          call lev%Q(m+1)%copy(this%rhs,1)
       end do  !!  End substep loop

       !  Set the value of qend
       !  If Gauss nodes, we must do integration
       !  unless the sweep was an initial Verlet
       !  For Lobatto nodes, we have a choice of whether to just use the
       !  value at the last node, or recompute it.
!!$       if (this%iqend .and. pf%state%iter .gt. 1) then
!!$          call Lev%encap%copy(Lev%qend, Lev%Q(1))
!!$          call Lev%encap%axpy(Lev%qend, dt, Lev%Q(1), 12)  !  Add the dt*v_0 term
!!$          m = nnodes
!!$          do n = 1, nnodes
!!$             call Lev%encap%axpy(Lev%qend, dt*this%Qmat(m,n), Lev%F(n,1),1)
!!$             call Lev%encap%axpy(Lev%qend, dtsq*this%QQmat(m,n), Lev%F(n,1),2)
!!$          end do
!!$          if (associated(Lev%tauQ)) then
!!$             call Lev%encap%axpy(Lev%qend, 1.0_pfdp, Lev%tauQ(nnodes-1))  
!!$             !          print *,'XXXXXXXXXXX  need code in verlet.f90'
!!$          end if
!!$       else
!!$          call Lev%encap%copy(Lev%qend, Lev%Q(nnodes))
!!$       end if
       
       call pf_residual(pf, lev, dt)
       call lev%qend%copy(lev%Q(lev%nnodes))
       
       call call_hooks(pf, level_index, PF_POST_SWEEP)
    
    end do ! end loop on sweeps


    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine verlet_sweep


  !-----------------------------------------------------------------------------
  ! Initialize integration matrices
  !
  subroutine verlet_initialize(this,lev)
    class(pf_verlet_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer :: i,j,nnodes,ierr

    this%npieces = 2
    nnodes = Lev%nnodes
    
    allocate(this%Qmat(nnodes,nnodes),stat=ierr)  !  0  to node integral
    if (ierr /=0) stop "allocate fail in verlet_initialize for Qmat"
    allocate(this%QQmat(nnodes,nnodes),stat=ierr)  !  0 to node double integral (like Qmat*Qmat)
    if (ierr /=0) stop "allocate fail in verlet_initialize for QQmat"


    allocate(this%Qtil(nnodes,nnodes),stat=ierr)  !  0  to node integral  approximation of Qmat (trap)
    if (ierr /=0) stop "allocate fail in verlet_initialize for Qtil"
    allocate(this%QQtil(nnodes-1,nnodes),stat=ierr)  !  0 to node QQmat  aproximation
    if (ierr /=0) stop "allocate fail in verlet_initialize for QQtil"
    allocate(this%QQLU(nnodes-1,nnodes),stat=ierr)  !  0 to node QQLU  aproximation
    if (ierr /=0) stop "allocate fail in verlet_initialize for QQLU"

    allocate(this%DQtil(nnodes-1,nnodes),stat=ierr)  !  0 to node Qmat-Qtil
    if (ierr /=0) stop "allocate fail in verlet_initialize for DQtil"
    allocate(this%DQQtil(nnodes-1,nnodes),stat=ierr)  !  node to node QQmat-QQtil
    if (ierr /=0) stop "allocate fail in verlet_initialize for DQQtil"

    allocate(this%bvec(nnodes),stat=ierr)  !  Integration rule for v
    if (ierr /=0) stop "allocate fail in verlet_initialize for DStil"
    allocate(this%bbarvec(nnodes),stat=ierr)  !  Integration rule for x
    if (ierr /=0) stop "allocate fail in verlet_initialize for DSStil"

    allocate(this%dtsdc(nnodes-1),stat=ierr)
    if (ierr /=0) stop "allocate fail in verlet_initialize for dtsdc"    
    allocate(this%tsdc(nnodes),stat=ierr)
    if (ierr /=0) stop "allocate fail in verlet_initialize for tsdc"    
    
    !>  Array of substep sizes
    this%dtsdc = lev%sdcmats%qnodes(2:nnodes) - lev%sdcmats%qnodes(1:nnodes-1)
    this%tsdc = lev%sdcmats%qnodes - lev%sdcmats%qnodes(1)

    !  Build Q from qmat
    this%Qmat=lev%sdcmats%qmat     !   I just use qmat now?

    !  The quadrature rule is the last row of Q
    this%bvec=this%Qmat(nnodes,:);

    !  form the QQ matrix depending on what you want
    select case (this%whichQQ)
    case (0)  !  Collocation (make it the product)
       print *,'Making QQ by collocation Q*Q'
       this%QQmat = matmul(this%Qmat,this%Qmat)
    case (1)  !  Make the pair like in Lobatto A/B pair
       print *,'Making QQ by collocation Lobatto pair'
       do i = 1,nnodes
          do j = 1,nnodes
             this%QQmat(i,j) =  this%bvec(j)*(1.0_pfdp-this%Qmat(j,i)/this%bvec(i))
          end do
       end do
       this%QQmat = matmul(this%Qmat,this%QQmat)
    case (2)  !  Make the pair like in Lobatto B/A pair
       print *,'Making QQ by collocation Lobatto pair'
       do i = 1,nnodes
          do j = 1,nnodes
             this%QQmat(i,j) =  this%bvec(j)*(1.0_pfdp-this%Qmat(j,i)/this%bvec(i))
          end do
       end do
       this%QQmat = matmul(this%QQmat,this%Qmat)
    case default
       write(*,*) 'bad case in verlet.f90, whichQQ = ',this%whichQQ
    end select


   ! 0 to node
   this%Qtil = 0.5_pfdp*(lev%sdcmats%qmatFE+lev%sdcmats%qmatBE)
   this%QQtil = matmul(lev%sdcmats%qmatFE,this%Qtil) + 0.5_pfdp*lev%sdcmats%qmatFE*lev%sdcmats%qmatFE

   !  Get LU matrices if desired
!!$    if (this%use_LUq .eq. 1) then 
!!$       print *,'Doing LU with doLU=',this%doLU
!!$       call myLUq(SDCmats%qmat,SDCmats%qmatLU,nnodes,0)
!!$       call pf_myLUexp(this%QQmat,L,U,nnodes,this%doLU)
!!$      this%QQLU=U
!!$      print *, 'U from LU',this%QQLU
!!$   else
        this%QQLU=this%QQtil  !  Normal verlet all the time     
!!$   end if
!!$
    !  Make differences
    this%DQQtil = this%QQmat-this%QQtil
    this%DQQLU = this%QQmat-this%QQLU
    this%DQtil = this%Qmat-this%Qtil

  end subroutine verlet_initialize
  
  !-----------------------------------------------------------------------------
  ! Integrate (t_n to node)
  subroutine verlet_integrate(this, lev, qSDC, fSDC, dt, fintSDC,flags)
    class(pf_verlet_t), intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev          !!  Current level
    class(pf_encap_t), intent(in   ) :: qSDC(:)      !!  Solution values
    class(pf_encap_t), intent(in   ) :: fSDC(:, :)   !!  RHS Function values
    real(pfdp),        intent(in   ) :: dt           !!  Time step
    class(pf_encap_t), intent(inout) :: fintSDC(:)   !!  Integral from t_n to t_m
    integer, optional, intent(in   ) :: flags    

    integer :: n, m

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       call fintSDC(n)%axpy(dt*this%tsdc(n+1), qSDC(1), 12)  !  Add the dt*v_0 term
       do m = 1, lev%nnodes
          call fintSDC(n)%axpy(dt*this%Qmat(n+1,m), fSDC(m,1),1)
          call fintSDC(n)%axpy( dt*dt*this%QQmat(n+1,m), fSDC(m,1),2)
       end do
    end do


  end subroutine verlet_integrate
  !-----------------------------------------------------------------------------
  ! Compute residual (t_n to node)
  subroutine verlet_residual(this, lev, dt, flags)
    class(pf_verlet_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev  !!  Current level
    real(pfdp),        intent(in   ) :: dt   !!  Time step
    integer, intent(in), optional   :: flags

    integer :: n, m


    call this%integrate(lev, lev%Q, lev%F, dt, lev%I, flags)
    

    ! add tau 
    if (allocated(lev%tauQ)) then
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), flags)
       end do
    end if
    do m = 1, lev%nnodes-1      
       call lev%R(m)%copy(lev%I(m))
       call lev%R(m)%axpy(1.0_pfdp, lev%Q(1))
       call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
    end do

 
  end subroutine verlet_residual

  !-----------------------------------------------------------------------------
  ! Integrate to fill qend
  !
!!$  subroutine verlet_qend_integrate(Lev, dt)
!!$    type(pf_level_t), intent(in) :: Lev
!!$    real(pfdp),       intent(in) :: dt
!!$
!!$
!!$    real(pfdp) :: dtsdc(1:Lev%nnodes-1)
!!$    integer :: nnodes, m
!!$    type(pf_verlet_t), pointer :: verlet
!!$    call c_f_pointer(Lev%sweeper%sweeperctx, verlet)
!!$
!!$    nnodes = Lev%nnodes
!!$  
!!$    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
!!$    call Lev%encap%copy(Lev%qend, Lev%Q(1))        
!!$    call Lev%encap%axpy(Lev%qend, dt, Lev%Q(1), 12)  !  Add the dt*v_0 term
!!$    do m = 1, Lev%nnodes
!!$       call Lev%encap%axpy(Lev%qend, dt*Lev%qmat(nnodes,m), Lev%F(m,1),1)
!!$       call Lev%encap%axpy(Lev%qend, dt*dt*thisSSmat(nnodes,m), Lev%F(m,1),2)
!!$    end do
!!$  end subroutine verlet_qend_integrate


  !-----------------------------------------------------------------------------
  ! Create Verlet sweeper
  !


  subroutine verlet_destroy(this,lev)
    class(pf_verlet_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev   !!  Current level

    deallocate(this%Qmat)
    deallocate(this%QQmat)


    deallocate(this%Qtil)
    deallocate(this%QQtil)
    deallocate(this%QQLU)    


    deallocate(this%DQtil)
    deallocate(this%DQQtil)

    deallocate(this%bvec)
    deallocate(this%bbarvec)

    deallocate(this%dtsdc)
    deallocate(this%tsdc)    


    call lev%ulevel%factory%destroy_single(this%rhs, lev%index,  lev%shape)    
  end subroutine verlet_destroy

  subroutine verlet_spreadq0(this, lev, t0, flags, step)
    class(pf_verlet_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t0
    integer, optional,   intent(in)    :: flags, step
    call pf_generic_spreadq0(this, lev, t0)
  end subroutine verlet_spreadq0

  !> Subroutine to evaluate function value at node m
  subroutine verlet_evaluate(this, lev, t, m, flags, step)

    class(pf_verlet_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev  !!  Current level
    real(pfdp),        intent(in   ) :: t    !!  Time at which to evaluate
    integer,           intent(in   ) :: m    !!  Node at which to evaluate
    integer, intent(in), optional   :: flags, step

       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1),1)
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2),2)
  end subroutine verlet_evaluate

  !> Subroutine to evaluate the function values at all nodes
  subroutine verlet_evaluate_all(this, lev, t, flags, step)
    class(pf_verlet_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev   !!  Current level
    real(pfdp),        intent(in   ) :: t(:)  !!  Array of times at each node
    integer, intent(in), optional   :: flags, step
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine verlet_evaluate_all
  

end module pf_mod_verlet

