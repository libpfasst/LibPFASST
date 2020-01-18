!!  Old style multi-implicit sweeper
!
! This file is part of LIBPFASST.
!
!>  Module for multi-implicit SDC sweeper type (old style), extends abstract sweeper class for an equation
!!   of the form
!!         $$   y' = f_1(y) + f_2(y) + f_2(y)  $$

module pf_mod_misdc
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  !>  MISDC sweeper type, extends abstract sweeper
  type, extends(pf_sweeper_t), abstract :: pf_misdc_t
     real(pfdp), allocatable :: SdiffE(:,:)
     real(pfdp), allocatable :: SdiffI(:,:)
   contains
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: sweep        => misdc_sweep
     procedure :: initialize   => misdc_initialize
     procedure :: evaluate     => misdc_evaluate
     procedure :: integrate    => misdc_integrate
     procedure :: residual     => misdc_residual
     procedure :: evaluate_all => misdc_evaluate_all
     procedure :: destroy      => misdc_destroy
     procedure :: misdc_destroy
  end type pf_misdc_t
  interface
     subroutine pf_f_eval_p(this,y, t, level, f, piece)
       import pf_misdc_t, pf_encap_t, pfdp
       class(pf_misdc_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y
       real(pfdp),        intent(in   ) :: t
       integer,    intent(in   ) :: level
       class(pf_encap_t), intent(inout) :: f
       integer,    intent(in   ) :: piece
     end subroutine pf_f_eval_p
      subroutine pf_f_comp_p(this,y, t, dt, rhs, level, f, piece)
       import pf_misdc_t, pf_encap_t, pfdp
       class(pf_misdc_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y
       real(pfdp),        intent(in   ) :: t
       real(pfdp),        intent(in   ) :: dt
       class(pf_encap_t), intent(in   ) :: rhs
       integer,    intent(in   ) :: level
       class(pf_encap_t), intent(inout) :: f
       integer,    intent(in   ) :: piece
     end subroutine pf_f_comp_p
  end interface

contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine misdc_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer
    class(pf_misdc_t), intent(inout)    :: this
    type(pf_pfasst_t),    intent(inout) :: pf
    real(pfdp),           intent(in)    :: dt, t0
    class(pf_level_t),     intent(inout) :: lev

    integer                        :: m, n
    real(pfdp)                     :: t
    real(pfdp)                     :: dtsdc(1:lev%nnodes-1)
    class(pf_encap_t), allocatable :: S3(:)
    class(pf_encap_t), allocatable :: rhs

    if (pf%save_timings > 1) call pf_start_timer(pf, T_SWEEP,lev%index)
    
    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1
       call lev%S(m)%setval(0.0_pfdp)
       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%SdiffE(m,n), lev%F(n,1))
          call lev%S(m)%axpy(dt*this%SdiffI(m,n), lev%F(n,2))
          call lev%S(m)%axpy(dt*lev%s0mat(m,n),   lev%F(n,3))
       end do
       if (allocated(lev%tau)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%copy(lev%q0)

    call misdc_evaluate(this, lev, t, 1)
    call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1),1)
    call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2),2)
    call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,3),3)

    call lev%ulevel%factory%create_single(rhs, lev%index,   lev%lev_shape)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       call rhs%copy(lev%Q(m))
       call rhs%axpy(dtsdc(m), lev%F(m,1))
       call rhs%axpy(1.0_pfdp, lev%S(m))

       call this%f_comp(lev%Q(m+1), t, dtsdc(m), rhs, lev%index, lev%F(m+1,2),2)

       !  Now we need to do the final subtraction for the f3 piece
       call rhs%copy(Lev%Q(m+1))       
       call rhs%axpy(-1.0_pfdp*dtsdc(m), lev%F(m+1,3))

       call this%f_comp(lev%Q(m+1), t, dtsdc(m), rhs, lev%index, lev%F(m+1,3),3)
       call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1),1)
       call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,2),2)
    end do
                         
    call lev%qend%copy(lev%Q(lev%nnodes))

    ! done
    call lev%ulevel%factory%destroy_single(rhs, lev%index,  lev%lev_shape)

    if (pf%save_timings > 1) call pf_stop_timer(pf, T_SWEEP,lev%index)

  end subroutine misdc_sweep
     

  ! Evaluate function values
  subroutine misdc_evaluate(this, lev, t, m)
    use pf_mod_dtype
    class(pf_misdc_t), intent(inout) :: this
    real(pfdp),        intent(in)    :: t
    integer,           intent(in)    :: m
    class(pf_level_t),  intent(inout) :: lev

    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1),1)
    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2),2)
    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,3),3)
  end subroutine misdc_evaluate


  ! Initialize matrices
  subroutine misdc_initialize(this, lev)
    class(pf_misdc_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m, nnodes,ierr

    this%npieces = 3

    nnodes = lev%nnodes
    allocate(this%SdiffE(nnodes-1,nnodes),stat=ierr)  !  S-FE
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)       
    allocate(this%SdiffI(nnodes-1,nnodes),stat=ierr)  !  S-BE
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)       

    this%SdiffE = lev%s0mat
    this%SdiffI = lev%s0mat

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       this%SdiffE(m,m)   = this%SdiffE(m,m)   - dsdc(m)
       this%SdiffI(m,m+1) = this%SdiffI(m,m+1) - dsdc(m)
    end do
  end subroutine misdc_initialize

  subroutine misdc_destroy(this, lev)
    class(pf_misdc_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%SdiffE)
    deallocate(this%SdiffI)
  end subroutine misdc_destroy
  
  ! Compute SDC integral
  subroutine misdc_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_misdc_t),  intent(inout) :: this
    class(pf_level_t),  intent(in)    :: lev
    class(pf_encap_t), intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in)    :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do    
  end subroutine misdc_integrate

  subroutine misdc_residual(this, lev, dt)
    class(pf_misdc_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),       intent(in)    :: dt

    call pf_generic_residual(this, lev, dt)
  end subroutine misdc_residual
  
  subroutine misdc_evaluate_all(this, lev, t)
    class(pf_misdc_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),       intent(in)    :: t(:)

    call pf_generic_evaluate_all(this, lev, t)
  end subroutine misdc_evaluate_all
  
end module pf_mod_misdc
