!! Multi-implicit sweeper module
!
! This file is part of LIBPFASST.
!
!>  Module of the  the derived sweeper class for doing MISDC sweeps for an equation of the form
!!       $$     y' = f_1(y) + f_2(y) + f_3(y) $$
!!  The \(f_1\) piece is treated explicitly and \(f_2\) and \(f_3\) implicitly
!!  Afer this sweeper is initialized (usually in main), the locgical flags can be changed if desired
module pf_mod_misdcQ
  use pf_mod_dtype
  use pf_mod_utils
  
  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_misdcQ_t
     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
     real(pfdp), allocatable :: dtsdc(:)

     class(pf_encap_t), allocatable :: I3(:)
     class(pf_encap_t), allocatable :: rhs
   contains 
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: sweep        => misdcQ_sweep
     procedure :: initialize   => misdcQ_initialize
     procedure :: integrate    => misdcQ_integrate
     procedure :: residual   => misdcQ_residual
     procedure :: spreadq0   => misdcQ_spreadq0
     procedure :: evaluate_all => misdcQ_evaluate_all
     procedure :: evaluate   => misdcQ_evaluate
     procedure :: destroy      => misdcQ_destroy
     procedure :: misdcQ_destroy
  end type pf_misdcQ_t

  interface
     !>  This is the interface for the routine to compute the RHS function values
     !>  Evaluate f_piece(y), where piece is one or two 
     subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
       !>  Evaluate f_piece(y), where piece is one or two 
       import pf_misdcQ_t, pf_encap_t, pfdp
       class(pf_misdcQ_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y        !!  Argument for evaluation
       real(pfdp),        intent(in   ) :: t        !!  Time at evaluation
       integer,    intent(in   ) :: level_index     !!  Level index
       class(pf_encap_t), intent(inout) :: f        !!  RHS function value
       integer,    intent(in   ) :: piece           !!  Which piece to evaluate
     end subroutine pf_f_eval_p
     !>  Solve the equation \(y - dtq*f_n(y) =rhs \) where n is given by the argument piece
     subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       import pf_misdcQ_t, pf_encap_t, pfdp
       class(pf_misdcQ_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y      !!  Solution of implicit solve 
       real(pfdp),        intent(in   ) :: t      !!  Time of solve
       real(pfdp),        intent(in   ) :: dtq    !!  dt*quadrature weight
       class(pf_encap_t), intent(in   ) :: rhs    !!  RHS for solve
       integer,    intent(in   ) :: level_index   !!  Level index
       class(pf_encap_t), intent(inout) :: f      !!  f_n of solution y
       integer,    intent(in   ) :: piece         !!  Which piece to evaluate
     end subroutine pf_f_comp_p
  end interface

contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine misdcQ_sweep(this, pf, level_index, t0, dt, nsweeps,flags)
    use pf_mod_timer
    use pf_mod_hooks    
    class(pf_misdcQ_t),      intent(inout) :: this
    type(pf_pfasst_t),target,intent(inout) :: pf  !!  PFASST structure
    integer,          intent(in) :: level_index   !!  which level to sweep on
    real(pfdp),       intent(in) :: t0            !!  Time at beginning of time step
    real(pfdp),       intent(in) :: dt            !!  time step size
    integer,          intent(in) :: nsweeps       !!  number of sweeps to do
    integer, optional, intent(in   ) :: flags    
    !>  Local variables    
    class(pf_level_t), pointer:: lev
    integer                        :: m, n,k
    real(pfdp)                     :: t

    call start_timer(pf, TLEVEL+lev%index-1)

    lev => pf%levels(level_index)   !!  Assign level pointer

    do k = 1,nsweeps   !!  Loop over sweeps
       call call_hooks(pf, level_index, PF_PRE_SWEEP)    

       ! compute integrals and add fas correction
       
       do m = 1, lev%nnodes-1
          
          call lev%I(m)%setval(0.0_pfdp)
          call this%I3(m)%setval(0.0d0)
          do n = 1, lev%nnodes
             call lev%I(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
             call lev%I(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2))
             call lev%I(m)%axpy(dt*lev%sdcmats%qmat(m,n),    lev%F(n,3))
             call this%I3(m)%axpy(dt*this%QtilI(m,n),     lev%F(n,3))
             !  Note we have to leave off the -dt*Qtil here and put it in after f2comp
          end do
          if (allocated(lev%tauQ)) then
             call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))
          end if
       end do
       
       ! do the time-stepping
       if (k .eq. 1) then
          call lev%Q(1)%copy(lev%q0)
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1),1)
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2),2)
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,3),3)
       endif
       
       t = t0
       do m = 1, lev%nnodes-1
          t = t + dt*this%dtsdc(m)
          
          call this%rhs%setval(0.0_pfdp)
          do n = 1, m
             call this%rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1))  
             call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2))  
          end do
          !  Add the tau term
          call this%rhs%axpy(1.0_pfdp, lev%I(m))
          !  Add the starting value
          call this%rhs%axpy(1.0_pfdp, lev%Q(1))
          
          call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index, lev%F(m+1,2),2)
          
          !  Now we need to do the final subtraction for the f3 piece
          call this%rhs%copy(Lev%Q(m+1))       
          do n = 1, m
             call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,3))  
          end do
          
          call this%rhs%axpy(-1.0_pfdp, this%I3(m))
          
          call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index, lev%F(m+1,3),3)
          call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1),1)
          call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,2),2)
       end do
       
       call pf_residual(pf, lev, dt)
       call lev%qend%copy(lev%Q(lev%nnodes))
       
       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do  !  End loop on sweeps
    call end_timer(pf, TLEVEL+lev%index-1)

  end subroutine misdcQ_sweep
     
  ! Initialize matrices
  subroutine misdcQ_initialize(this, lev)
    class(pf_misdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m, n, nnodes

    this%npieces = 3

    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes)) ! S-FE
    allocate(this%QdiffI(nnodes-1,nnodes)) ! S-BE 
    allocate(this%QtilE(nnodes-1,nnodes)) ! S-FE
    allocate(this%QtilI(nnodes-1,nnodes)) ! S-BE
    allocate(this%dtsdc(nnodes-1))
    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp
    
        !>  Array of substep sizes
    this%dtsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    ! Implicit matrix
    if (this%use_LUq) then 
       ! Get the LU
       this%QtilI = lev%sdcmats%qmatLU

    else 
       this%QtilI = lev%sdcmats%qmatBE
    end if
    ! Explicit matrix    
    this%QtilE=lev%sdcmats%qmatFE

    this%QdiffE = lev%sdcmats%qmat-this%QtilE
    this%QdiffI = lev%sdcmats%qmat-this%QtilI

    !>  Make space for rhs
    call lev%ulevel%factory%create_single(this%rhs, lev%index,   lev%shape)

    !>  Make space for extra integration piece
    call lev%ulevel%factory%create_array(this%I3,lev%nnodes-1,lev%index,lev%shape)

  end subroutine misdcQ_initialize

  subroutine misdcQ_destroy(this, lev)
    class(pf_misdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
    deallocate(this%dtsdc)


    call lev%ulevel%factory%destroy_array(this%I3,lev%nnodes-1,lev%index,lev%shape)
    call lev%ulevel%factory%destroy_single(this%rhs, lev%index,  lev%shape)
  end subroutine misdcQ_destroy

  ! Compute SDC integral
  subroutine misdcQ_integrate(this, lev, qSDC, fSDC, dt, fintSDC,flags)
    class(pf_misdcQ_t),  intent(inout) :: this
    class(pf_level_t),  intent(in)    :: lev
    class(pf_encap_t), intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in)    :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)
    integer, optional, intent(in   ) :: flags    

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%sdcmats%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do    
  end subroutine misdcQ_integrate

  !> Subroutine to evaluate function value at node m
  subroutine misdcQ_evaluate(this, lev, t, m, flags, step)

    class(pf_misdcQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev  !!  Current level
    real(pfdp),        intent(in   ) :: t    !!  Time at which to evaluate
    integer,           intent(in   ) :: m    !!  Node at which to evaluate
    integer, intent(in), optional   :: flags, step

       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1),1)
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2),2)
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,3),3)       
  end subroutine misdcQ_evaluate

  !> Subroutine to evaluate the function values at all nodes
  subroutine misdcQ_evaluate_all(this, lev, t, flags, step)
    class(pf_misdcQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev   !!  Current level
    real(pfdp),        intent(in   ) :: t(:)  !!  Array of times at each node
    integer, intent(in), optional   :: flags, step
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine misdcQ_evaluate_all

  !> Subroutine to compute  Residual
  subroutine misdcQ_residual(this, lev, dt, flags)
    class(pf_misdcQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev  !!  Current level
    real(pfdp),        intent(in   ) :: dt   !!  Time step
    integer, intent(in), optional   :: flags
    call pf_generic_residual(this, lev, dt)
  end subroutine misdcQ_residual

  subroutine misdcQ_spreadq0(this, lev, t0, flags, step)
    class(pf_misdcQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t0
    integer, optional,   intent(in)    :: flags, step
    call pf_generic_spreadq0(this, lev, t0)
  end subroutine misdcQ_spreadq0

end module pf_mod_misdcQ
