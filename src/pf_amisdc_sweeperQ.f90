!!  Asynchronous MISDC sweeper
!
! This file is part of LIBPFASST.
!
!> Module for Asynchronous multi-implicit sweeper
module pf_mod_amisdcQ
  use pf_mod_amisdc
  implicit none

  !>  Asynchronous multi-implicit sweeper type 
  type, extends(pf_amisdc_t), abstract :: pf_amisdcQ_t
     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
     logical                 :: use_LUq_ = .true.
   contains 
     procedure :: sweep        => amisdcQ_sweep
     procedure :: initialize   => amisdcQ_initialize
     procedure :: integrate    => amisdcQ_integrate
     procedure :: destroy      => amisdcQ_destroy
     procedure :: sweep_coupled_implicit_terms
     procedure :: sweep_decoupled_implicit_terms
     procedure :: amisdcQ_destroy
     procedure :: amisdcQ_initialize
  end type pf_amisdcQ_t

contains

  ! Perform an SDC sweep on level lev and set qend appropriately.
  ! In the asynchronous updates, the two implicit parts are coupled
  subroutine sweep_coupled_implicit_terms(this, pf, lev, t0, dt)
    use pf_mod_timer
    class(pf_amisdcQ_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout) :: pf
    real(pfdp),          intent(in)    :: dt, t0
    class(pf_level_t),   intent(inout) :: lev

    integer                        :: m, n
    real(pfdp)                     :: t
    real(pfdp)                     :: dtsdc(1:lev%nnodes-1)
    class(pf_encap_t), allocatable :: rhsA, rhsB, QA, QB
    class(pf_encap_t), allocatable :: S2(:), S3(:)


    call pf_start_timer(pf, T_SWEEP,lev%index)
   
    call lev%ulevel%factory%create_array(S2,lev%nnodes-1,lev%index,lev%lev_shape)
    call lev%ulevel%factory%create_array(S3,lev%nnodes-1,lev%index,lev%lev_shape)
    
    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1

       call lev%S(m)%setval(0.0_pfdp)
       call S2(m)%setval(0.0d0)
       call S3(m)%setval(0.0d0)

       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
          call lev%S(m)%axpy(1.0_pfdp*dt*lev%qmat(m,n),    lev%F(n,2))
          call lev%S(m)%axpy(1.0_pfdp*dt*lev%qmat(m,n),    lev%F(n,3))
          call S2(m)%axpy(2.0_pfdp*dt*this%QtilI(m,n),     lev%F(n,2))
          call S3(m)%axpy(2.0_pfdp*dt*this%QtilI(m,n),     lev%F(n,3))
       end do
       if (allocated(lev%tauQ)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%copy(lev%q0)

    call this%f1eval(lev%Q(1), t0, lev%index, lev%F(1,1))
    call this%f2eval(lev%Q(1), t0, lev%index, lev%F(1,2))
    call this%f3eval(lev%Q(1), t0, lev%index, lev%F(1,3))

    call lev%ulevel%factory%create_single(rhsA, lev%index,   lev%lev_shape)
    call lev%ulevel%factory%create_single(rhsB, lev%index,   lev%lev_shape)
    call lev%ulevel%factory%create_single(QA,   lev%index,   lev%lev_shape)
    call lev%ulevel%factory%create_single(QB,   lev%index,   lev%lev_shape)

    call QA%setval(0.0_pfdp)
    call QB%setval(0.0_pfdp)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)
             
       call rhsA%copy(lev%Q(1))
       ! First compute the explicit part of the right-hand side
       do n = 1, m
          call rhsA%axpy(dt*this%QtilE(m,n), lev%F(n,1))  
       end do
       call rhsA%axpy(1.0_pfdp, lev%S(m))
   
       ! Save the right-hand side with only the explicit contribution
       call rhsB%copy(rhsA)

       ! Add the first implicit part to the right-hand side and solve for the first asynchronous update
       do n = 1, m
          call rhsA%axpy(2.0_pfdp*dt*this%QtilI(m,n), lev%F(n,2))  
       end do
       call rhsA%axpy(-1.0_pfdp, S2(m))  
       call this%f2comp(QA, t, 2.0_pfdp*dt*this%QtilI(m,m+1), rhsA, lev%index, lev%F(m+1,2))

       ! Add the second implicit part to the right-hand side and solve for the second asynchronous update
       do n = 1, m
          call rhsB%axpy(2.0_pfdp*dt*this%QtilI(m,n), lev%F(n,3))  
       end do
       call rhsB%axpy(-1.0_pfdp, S3(m))  
       call this%f3comp(QB, t, 2.0_pfdp*dt*this%QtilI(m,m+1), rhsB, lev%index, lev%F(m+1,3))

       ! Now we average the two asynchronous updates
       call lev%Q(m+1)%setval(0.0_pfdp)
       call lev%Q(m+1)%axpy(0.5_pfdp, QA)
       call lev%Q(m+1)%axpy(0.5_pfdp, QB)

       ! Evaluate the three right-hand sides with the updated variables
       call this%f1eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1))
       call this%f2eval(lev%Q(m+1), t, lev%index, lev%F(m+1,2))
       call this%f3eval(lev%Q(m+1), t, lev%index, lev%F(m+1,3))
    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    call lev%ulevel%factory%destroy_array(S2)
    call lev%ulevel%factory%destroy_array(S3)
    call lev%ulevel%factory%destroy_single(rhsA)
    call lev%ulevel%factory%destroy_single(rhsB)
    call lev%ulevel%factory%destroy_single(QA)
    call lev%ulevel%factory%destroy_single(QB)

    call pf_stop_timer(pf, T_SWEEP,lev%index)


  end subroutine sweep_coupled_implicit_terms

  
  ! Perform an SDC sweep on level lev and set qend appropriately.
  ! In the asynchronous updates, the two implicit parts are decoupled

  ! (in progress)
  subroutine sweep_decoupled_implicit_terms(this, pf, lev, t0, dt)
    use pf_mod_timer
    class(pf_amisdcQ_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout) :: pf
    real(pfdp),          intent(in)    :: dt, t0
    class(pf_level_t),   intent(inout) :: lev

    ! integer                        :: m, n
    ! real(pfdp)                     :: t
    ! real(pfdp)                     :: dtsdc(1:lev%nnodes-1)
    ! class(pf_encap_t), allocatable :: rhsA, rhsB, QA, QB
    ! class(pf_encap_t), allocatable :: S2(:), S3(:)

    ! call start_timer(pf, TLEVEL+lev%index-1)
   
    ! call lev%ulevel%factory%create_array(S2,lev%nnodes-1,lev%index,lev%lev_shape)
    ! call lev%ulevel%factory%create_array(S3,lev%nnodes-1,lev%index,lev%lev_shape)
    
    ! ! compute integrals and add fas correction
    ! do m = 1, lev%nnodes-1

    !    call lev%S(m)%setval(0.0_pfdp)
    !    call S2(m)%setval(0.0d0)
    !    call S3(m)%setval(0.0d0)

    !    do n = 1, lev%nnodes
    !       call lev%S(m)%axpy(dt*this%QdiffE(m,n),       lev%F(n,1))
    !       call S2(m)%axpy( 2.0_pfdp*dt*lev%qmat(m,n),   lev%F(n,2))
    !       call S3(m)%axpy( 2.0_pfdp*dt*lev%qmat(m,n),   lev%F(n,3))
    !       call S2(m)%axpy(-2.0_pfdp*dt*this%QtilI(m,n), lev%F(n,2))
    !       call S3(m)%axpy(-2.0_pfdp*dt*this%QtilI(m,n), lev%F(n,3))
    !    end do
    !    if (allocated(lev%tauQ)) then
    !       call lev%S(m)%axpy(1.0_pfdp, lev%tauQ(m))
    !    end if
    ! end do

    ! ! do the time-stepping
    ! call lev%Q(1)%unpack(lev%q0)

    ! call this%f1eval(lev%Q(1), t0, lev%index, lev%F(1,1))
    ! call this%f2eval(lev%Q(1), t0, lev%index, lev%F(1,2))
    ! call this%f3eval(lev%Q(1), t0, lev%index, lev%F(1,3))

    ! call lev%ulevel%factory%create_single(rhsA, lev%index,  lev%lev_shape)
    ! call lev%ulevel%factory%create_single(rhsB, lev%index,  lev%lev_shape)
    ! call lev%ulevel%factory%create_single(QA,   lev%index,  lev%lev_shape)
    ! call lev%ulevel%factory%create_single(QB,   lev%index,  lev%lev_shape)

    ! call QA%setval(0.0_pfdp)
    ! call QB%setval(0.0_pfdp)

    ! t = t0
    ! dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    ! do m = 1, lev%nnodes-1
    !    t = t + dtsdc(m)
             
    !    call rhsA%copy(lev%Q(1))
    !    ! First compute the explicit part of the right-hand side
    !    do n = 1, m
    !       call rhsA%axpy(dt*this%QtilE(m,n), lev%F(n,1))  
    !    end do
    !    call rhsA%axpy(1.0_pfdp, lev%S(m))
   
    !    ! Save the right-hand side with only the explicit contribution
    !    call rhsB%copy(rhsA)

    !    ! Add the first implicit part to the right-hand side and solve for the first asynchronous update
    !    do n = 1, m
    !       call rhsA%axpy(2.0_pfdp*dt*this%QtilI(m,n), lev%F(n,2))  
    !    end do
    !    call rhsA%axpy(1.0_pfdp, S2(m))  
    !    call this%f2comp(QA, t, 2.0_pfdp*dt*this%QtilI(m,m+1), rhsA, lev%index, lev%F(m+1,2))

    !    ! Add the second implicit part to the right-hand side and solve for the second asynchronous update
    !    do n = 1, m
    !       call rhsB%axpy(2.0_pfdp*dt*this%QtilI(m,n), lev%F(n,3))  
    !    end do
    !    call rhsB%axpy(1.0_pfdp, S3(m))  
    !    call this%f3comp(QB, t, 2.0_pfdp*dt*this%QtilI(m,m+1), rhsB, lev%index, lev%F(m+1,3))

    !    ! Now we average the two asynchronous updates
    !    call lev%Q(m+1)%setval(0.0_pfdp)
    !    call lev%Q(m+1)%axpy(0.5_pfdp, QA)
    !    call lev%Q(m+1)%axpy(0.5_pfdp, QB)

    !    ! Evaluate the three right-hand sides with the updated variables
    !    call this%f1eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1))
    !    call this%f2eval(lev%Q(m+1), t, lev%index, lev%F(m+1,2))
    !    call this%f3eval(lev%Q(m+1), t, lev%index, lev%F(m+1,3))
    ! end do

    ! call lev%qend%copy(lev%Q(lev%nnodes))

    ! call lev%ulevel%factory%destroy_array(S2,lev%nnodes-1,lev%index,lev%lev_shape)
    ! call lev%ulevel%factory%destroy_array(S3,lev%nnodes-1,lev%index,lev%lev_shape)
    ! call lev%ulevel%factory%destroy_single(rhsA, lev%index,   lev%lev_shape)
    ! call lev%ulevel%factory%destroy_single(rhsB, lev%index,   lev%lev_shape)
    ! call lev%ulevel%factory%destroy_single(QA,   lev%index,   lev%lev_shape)
    ! call lev%ulevel%factory%destroy_single(QB,   lev%index,   lev%lev_shape)

        
  end subroutine sweep_decoupled_implicit_terms


  ! Perform an SDC sweep on level lev and set qend appropriately.
  subroutine amisdcQ_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer
    class(pf_amisdcQ_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout) :: pf
    real(pfdp),          intent(in)    :: dt, t0
    class(pf_level_t),   intent(inout) :: lev

    call sweep_coupled_implicit_terms(this, pf, lev, t0, dt)

  end subroutine amisdcQ_sweep
    
  ! Initialize matrices
  subroutine amisdcQ_initialize(this, lev)
    class(pf_amisdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m, n, nnodes,ierr

    this%npieces = 3

    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes),stat=ierr)  !  S-FE
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
    allocate(this%QdiffI(nnodes-1,nnodes),stat=ierr)  !  S-BE
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
    allocate(this%QtilE(nnodes-1,nnodes),stat=ierr)  !  S-FE
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
    allocate(this%QtilI(nnodes-1,nnodes),stat=ierr)  !  S-BE
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)

    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    ! Implicit matrix
    if (this%use_LUq_) then 
       ! Get the LU
       call myLUq(lev%qmat,lev%LUmat,lev%nnodes,1)
       this%QtilI = lev%LUmat
    else 
       do m = 1, nnodes-1
          do n = 1,m
             this%QtilI(m,n+1) =  dsdc(n)
          end do
       end do
    end if
    ! Explicit matrix
    do m = 1, nnodes-1
       do n = 1,m
          this%QtilE(m,n)   =  dsdc(n)
       end do
    end do

    this%QdiffE = lev%qmat-this%QtilE
    this%QdiffI = lev%qmat-this%QtilI
  end subroutine amisdcQ_initialize

  ! Destroy the matrices
  subroutine amisdcQ_destroy(this, lev)
    class(pf_amisdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
  end subroutine amisdcQ_destroy


  ! Compute SDC integral
  subroutine amisdcQ_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_amisdcQ_t),  intent(inout) :: this
    class(pf_level_t),  intent(in)      :: lev
    class(pf_encap_t), intent(in)       :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in)       :: dt
    class(pf_encap_t), intent(inout)    :: fintSDC(:)

    integer :: n, m, p
    
    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do  
  end subroutine amisdcQ_integrate
 
end module pf_mod_amisdcQ
