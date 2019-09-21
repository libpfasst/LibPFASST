!! Multi-implicit forward-backward sweeper module
!
! This file is part of LIBPFASST.
!
!>  Module of the  the derived sweeper class for doing MISDC sweeps for an equation of the form
!!       $$     y' = f_1(y) + f_2(y) + f_3(y) $$
!!  The \(f_1\) piece is treated explicitly and \(f_2\) and \(f_3\) implicitly
!!  Afer this sweeper is initialized (usually in main), the locgical flags can be changed if desired
module pf_mod_misdcQ_oc
  use pf_mod_dtype
  use pf_mod_utils

  implicit none

  !>  Multi-implicit SDC sweeper type for optimal control, extends abstract sweeper
  type, extends(pf_sweeper_t), abstract :: pf_misdcQ_oc_t
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
     procedure :: sweep        => misdcQ_oc_sweep
     procedure :: initialize   => misdcQ_oc_initialize
     procedure :: integrate    => misdcQ_oc_integrate
     procedure :: residual     => misdcQ_oc_residual
     procedure :: spreadq0     => misdcQ_oc_spreadq0
     procedure :: evaluate_all => misdcQ_oc_evaluate_all
     procedure :: evaluate     => misdcQ_oc_evaluate
     procedure :: destroy      => misdcQ_oc_destroy
     procedure :: misdcQ_oc_destroy
     procedure :: misdcQ_oc_initialize
  end type pf_misdcQ_oc_t

  interface
     !>  This is the interface for the routine to compute the RHS function values
     !>  Evaluate f_piece(y), where piece is one or two
     subroutine pf_f_eval_p(this,y, t, level_index, f, piece, flags, idx, step)
       !>  Evaluate f_piece(y), where piece is one or two
       import pf_misdcQ_oc_t, pf_encap_t, pfdp
       class(pf_misdcQ_oc_t), intent(inout) :: this
       class(pf_encap_t),     intent(in   ) :: y            !!  Argument for evaluation
       real(pfdp),            intent(in   ) :: t            !!  Time at evaluation
       integer,               intent(in   ) :: level_index  !!  Level index
       class(pf_encap_t),     intent(inout) :: f            !!  RHS function value
       integer,               intent(in   ) :: piece        !!  Which piece to evaluate
       integer,               intent(in   ) :: flags        !!  forward or backward
       integer,   intent(in), optional    :: idx            !! index of quadrature node
       integer,   intent(in), optional    :: step           !! time step for sequential version
     end subroutine pf_f_eval_p
     !>  Solve the equation \(y - dtq*f_n(y) =rhs \) where n is given by the argument piece
     subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece, flags)
       import pf_misdcQ_oc_t, pf_encap_t, pfdp
       class(pf_misdcQ_oc_t), intent(inout) :: this
       class(pf_encap_t),     intent(inout) :: y           !!  Solution of implicit solve
       real(pfdp),            intent(in   ) :: t           !!  Time of solve
       real(pfdp),            intent(in   ) :: dtq         !!  dt*quadrature weight
       class(pf_encap_t),     intent(in   ) :: rhs         !!  RHS for solve
       integer,               intent(in   ) :: level_index !!  Level index
       class(pf_encap_t),     intent(inout) :: f           !!  f_n of solution y
       integer,               intent(in   ) :: piece       !!  Which piece to evaluate
       integer,               intent(in   ) :: flags       !!  forward or backward
     end subroutine pf_f_comp_p
  end interface

contains

  ! Perform one forward and/or backward SDC sweep on level and set qend/q0 appropriately.
  subroutine misdcQ_oc_sweep(this, pf, level_index, t0, dt, nsweeps, flags)
    use pf_mod_timer
    use pf_mod_hooks
    class(pf_misdcQ_oc_t),      intent(inout) :: this
    type(pf_pfasst_t), target,  intent(inout) :: pf
    integer,                    intent(in   ) :: level_index, nsweeps
    real(pfdp),                 intent(in   ) :: dt, t0
    integer,          optional, intent(in   ) :: flags
    !>  Local variables
    class(pf_level_t), pointer :: lev
    integer     :: k, m, n, which, Nnodes
    real(pfdp)  :: t, tend
    logical     :: sweep_y, sweep_p
    integer     :: step

    lev => pf%levels(level_index)   !  Assign level pointer

    call start_timer(pf, TLEVEL+lev%index-1)
    step = pf%state%step+1

    which = 0
    if (present(flags)) which = flags

    if (which .eq. 1) then
        sweep_y = .true.
        sweep_p = .false.
    else if (which .eq. 2) then
        sweep_y = .false.
        sweep_p = .true.
    else
       sweep_y = .true.
       sweep_p = .true.
    end if

    Nnodes = lev%nnodes
    tend = t0+dt

    do k = 1,nsweeps
       pf%state%sweep=k
       call call_hooks(pf, level_index, PF_PRE_SWEEP)

       ! compute integrals and add fas correction
        if( sweep_y ) then
        !  Forward in y
          do m = 1, Nnodes-1
            call lev%I(m)%setval(0.0_pfdp,1)
            call this%I3(m)%setval(0.0_pfdp,1)
            
            do n = 1, Nnodes
               call lev%I(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1), 1)
               call lev%I(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2), 1)
               call lev%I(m)%axpy(dt*lev%sdcmats%qmat(m,n), lev%F(n,3), 1)
               call this%I3(m)%axpy(dt*this%QtilI(m,n), lev%F(n,3), 1)
            end do
            if (level_index < pf%state%finest_level) then
               call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m),1)
            end if
          end do
        end if

        if( sweep_p ) then
        !  Backward in p
          do m =  Nnodes-1,1,-1
            call lev%I(m)%setval(0.0_pfdp,2)
            call this%I3(m)%setval(0.0_pfdp,2)

            do n = Nnodes,1,-1
               call lev%I(m)%axpy(dt*this%QdiffE(Nnodes-m,Nnodes+1-n), lev%F(n,1), 2)
               call lev%I(m)%axpy(dt*this%QdiffI(Nnodes-m,Nnodes+1-n), lev%F(n,2), 2)
               call lev%I(m)%axpy(dt*lev%sdcmats%qmat(Nnodes-m,Nnodes+1-n), lev%F(n,3), 2)
               call this%I3(m)%axpy(dt*this%QtilI(Nnodes-m,Nnodes+1-n), lev%F(n,3), 2)
            end do
            if (level_index < pf%state%finest_level) then
               call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), 2)
            end if
          end do
        end if

       ! do the time-stepping
       if (k .eq. 1) then
         if( sweep_y ) then
            call lev%Q(1)%copy(lev%q0, 1)
            call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1), 1, 1, 1, step)
            call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2), 2, 1, 1, step)
            call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,3), 3, 1, 1, step)
         end if
         if( sweep_p ) then
            call lev%Q(Nnodes)%copy(lev%qend, 2)
            call this%f_eval(lev%Q(Nnodes), tend, lev%index, lev%F(Nnodes,1), 1, 2, Nnodes, step)
            call this%f_eval(lev%Q(Nnodes), tend, lev%index, lev%F(Nnodes,2), 2, 2, Nnodes, step)
            call this%f_eval(lev%Q(Nnodes), tend, lev%index, lev%F(Nnodes,3), 3, 2, Nnodes, step)
         end if
       end if ! k .eq. 1

       if (sweep_y) then
         t = t0
         do m = 1, lev%nnodes-1
            t = t + dt*this%dtsdc(m)

            call this%rhs%setval(0.0_pfdp,1)
            do n = 1, m
               call this%rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1), 1)
               call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2), 1)
            end do
            !  Add the tau term
            call this%rhs%axpy(1.0_pfdp, lev%I(m), 1)
            !  Add the starting value
            call this%rhs%axpy(1.0_pfdp, lev%Q(1), 1)

            call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index, lev%F(m+1,2), 2, 1)

            !  Now we need to do the final subtraction for the f3 piece
            call this%rhs%copy(lev%Q(m+1), 1)
            do n = 1, m
               call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,3), 1)
            end do

            call this%rhs%axpy(-1.0_pfdp, this%I3(m), 1)
            call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index, lev%F(m+1,3), 3, 1)
            call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1), 1, 1, m+1, step)
            call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,2), 2, 1, m+1, step)
         end do
         !call pf_residual(pf, level_index, dt, 1)
         call lev%qend%copy(lev%Q(lev%nnodes), 1)
       end if ! sweep_y

       if (sweep_p) then
         t = tend
         do m = Nnodes-1, 1, -1
            t = t - dt*this%dtsdc(m)
            call this%rhs%setval(0.0_pfdp,2)
            do n = Nnodes, m+1,-1
               call this%rhs%axpy(dt*this%QtilE(Nnodes-m,Nnodes-n+1), lev%F(n,1), 2)
               call this%rhs%axpy(dt*this%QtilI(Nnodes-m,Nnodes-n+1), lev%F(n,2), 2)
            end do
            !  Add the tau term
            call this%rhs%axpy(1.0_pfdp, lev%I(m), 2)
            !  Add the starting value
            call this%rhs%axpy(1.0_pfdp, lev%Q(Nnodes), 2)

            call this%f_comp(lev%Q(m), t, dt*this%QtilI(Nnodes-m,Nnodes-m+1), this%rhs, lev%index, lev%F(m,2), 2, 2)

            !  Now we need to do the final subtraction for the f3 piece
            call this%rhs%copy(lev%Q(m), 2)
            do n = Nnodes, m+1,-1
               call this%rhs%axpy(dt*this%QtilI(Nnodes-m,Nnodes-n+1), lev%F(n,3), 2)
            end do

            call this%rhs%axpy(-1.0_pfdp, this%I3(m), 2)

            call this%f_comp(lev%Q(m), t, dt*this%QtilI(Nnodes-m,Nnodes-m+1), this%rhs, lev%index, lev%F(m,3), 3, 2)
            call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1), 1, 2, m, step)
            call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2), 2, 2, m, step)
         end do
         !call pf_residual(pf, level_index, dt, 2)
         call lev%q0%copy(lev%Q(1), 2)
       end if  ! sweep_p

       if( sweep_p .and. sweep_y ) then
         call pf_residual(pf, level_index, dt, 0)
       else if( sweep_y ) then
         call pf_residual(pf, level_index, dt, 1)
       else if (sweep_p ) then
         call pf_residual(pf, level_index, dt, 2)
       else
         stop "neither sweep on p nor on y : that should not happen"
       end if
       ! done
       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do ! k=1,nsweeps
  end subroutine misdcQ_oc_sweep

  ! Initialize matrices
  subroutine misdcQ_oc_initialize(this, pf,level_index)
    class(pf_misdcQ_oc_t), intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    integer    :: m,n, nnodes
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    this%npieces = 3

    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(this%QtilE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QtilI(nnodes-1,nnodes))  !  S-BE
    allocate(this%dtsdc(nnodes-1))
    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

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
    call lev%ulevel%factory%create_single(this%rhs, lev%index, lev%lev_shape)

    !>  Make space for extra integration piece
    call lev%ulevel%factory%create_array(this%I3,lev%nnodes-1,lev%index,lev%lev_shape)

  end subroutine misdcQ_oc_initialize

  subroutine misdcQ_oc_destroy(this, pf,level_index)
    class(pf_misdcQ_oc_t), intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
    deallocate(this%dtsdc)


    call lev%ulevel%factory%destroy_array(this%I3)
    call lev%ulevel%factory%destroy_single(this%rhs)
  end subroutine misdcQ_oc_destroy

  ! Compute SDC integral
  subroutine misdcQ_oc_integrate(this, pf,level_index, qSDC, fSDC, dt, fintSDC, flags)
    class(pf_misdcQ_oc_t),  intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    class(pf_encap_t),      intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),             intent(in)    :: dt
    class(pf_encap_t),       intent(inout) :: fintSDC(:)
    integer,      optional, intent(in)    :: flags

    integer :: n, m, p, which
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    which = 0
    if(present(flags)) which = flags

    do n = 1, lev%nnodes-1
       if( (which .eq. 0) .or. (which .eq. 1) ) then
          call fintSDC(n)%setval(0.0_pfdp, 1)
          do m = 1, lev%nnodes
            do p = 1, this%npieces
              call fintSDC(n)%axpy(dt*lev%sdcmats%qmat(n,m), fSDC(m,p), 1)
            end do
          end do
       end if

       !  Backward in p
       if( (which .eq. 0) .or. (which .eq. 2) ) then
          call fintSDC(lev%nnodes-n)%setval(0.0_pfdp, 2)
          do m = 1, lev%nnodes
            do p = 1, this%npieces
              call fintSDC(lev%nnodes-n)%axpy(dt*lev%sdcmats%qmat(n,m), fSDC(lev%nnodes+1-m,p), 2)
            end do
          end do
       end if
    end do
  end subroutine misdcQ_oc_integrate

  ! Evaluate function values
  subroutine misdcQ_oc_evaluate(this, pf,level_index, t, m, flags, step)
    class(pf_misdcQ_oc_t),  intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),             intent(in   ) :: t
    integer,                intent(in   ) :: m
    integer,      optional, intent(in   ) :: flags, step

    integer  :: which, mystep
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    
    which = 0
    if (present(flags)) which = flags

    mystep = 1
    if(present(step)) mystep = step

    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1), 1, which, m, step)
    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2), 2, which, m, step)
    call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,3), 3, which, m, step)
  end subroutine misdcQ_oc_evaluate


  subroutine misdcQ_oc_evaluate_all(this, pf,level_index, t, flags, step)
    !! Evaluate all function values
    class(pf_misdcQ_oc_t),  intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),        intent(in   ) :: t(:)
    integer, intent(in), optional   :: flags, step

    integer :: m
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer

    if (.not.present(flags)) stop "MISDCQ_OC EVAL_ALL WITHOUT FLAGS"
    if (.not.present(step)) stop "MISDCQ_OC EVAL_ALL WITHOUT step"


    do m = 1, lev%nnodes
      call this%evaluate(pf,level_index, t(m), m, flags, step)
    end do
  end subroutine misdcQ_oc_evaluate_all

  subroutine misdcQ_oc_residual(this, pf, level_index, dt, flags)
    class(pf_misdcQ_oc_t),  intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),             intent(in)    :: dt
    integer,      optional, intent(in)    :: flags

    integer :: m, n, which
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer
    
    which = 0
    if(present(flags)) which = flags

    call this%integrate(pf,level_index, pf%levels(level_index)%Q, pf%levels(level_index)%F, dt, &
         pf%levels(level_index)%I, which)

    
    ! add tau (which is 'node to node')
    if (level_index < pf%state%finest_level) then
       do m = 1, pf%levels(level_index)%nnodes-1
          call pf%levels(level_index)%I(m)%axpy(1.0_pfdp, pf%levels(level_index)%tauQ(m), which)
       end do
    end if

    ! subtract out Q
    do m = 1, pf%levels(level_index)%nnodes-1
       if( (which .eq. 0) .or. (which .eq. 1) ) then
          call pf%levels(level_index)%R(m)%copy(pf%levels(level_index)%I(m), 1)
          call pf%levels(level_index)%R(m)%axpy(1.0_pfdp, pf%levels(level_index)%Q(1), 1)
          call pf%levels(level_index)%R(m)%axpy(-1.0_pfdp, pf%levels(level_index)%Q(m+1), 1)
       end if
       if( (which .eq. 0) .or. (which .eq. 2) ) then
          call pf%levels(level_index)%R(m)%copy(pf%levels(level_index)%I(m), 2)
          call pf%levels(level_index)%R(m)%axpy(1.0_pfdp, pf%levels(level_index)%Q(pf%levels(level_index)%nnodes), 2)
          call pf%levels(level_index)%R(m)%axpy(-1.0_pfdp, pf%levels(level_index)%Q(m), 2)
       end if
    end do

  end subroutine misdcQ_oc_residual

  subroutine misdcQ_oc_spreadq0(this, pf,level_index, t0, flags, step)
    class(pf_misdcQ_oc_t), intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),            intent(in   ) :: t0
    integer,     optional, intent(in)    :: flags, step

    integer :: m, p, which, mystep
    type(pf_level_t), pointer :: lev
    lev => pf%levels(level_index)   !  Assign level pointer


    which = 3
    if(present(flags)) which = flags
    if (.not.present(flags)) stop "IMEXQ_OC SPREADQ0 WITHOUT FLAGS"

    mystep = 1
    if(present(step))  then
      mystep = step !needed for sequential version
    else
      print *, "step not present in spreadq0", which
      stop
    end if

    select case(which)
      case(1)
        !  Stick initial condition into first node slot
        call lev%Q(1)%copy(lev%q0, 1)
        !  Evaluate F at first spot
        call lev%ulevel%sweeper%evaluate(pf,level_index, t0, 1, 1, mystep)
        ! Spread F and solution to all nodes
        do m = 2, lev%nnodes
          call lev%Q(m)%copy(lev%Q(1), 1)
          do p = 1, lev%ulevel%sweeper%npieces
            call lev%F(m,p)%copy(lev%F(1,p), 1)
          end do
        end do
      case(2)
        !  Stick terminal condition into last node slot
        call lev%Q(lev%nnodes)%copy(lev%qend, 2)
        !  Evaluate F at first spot
        call lev%ulevel%sweeper%evaluate(pf,level_index, t0, lev%nnodes, 2, mystep)
        ! Spread F and solution to all nodes
        do m = lev%nnodes-1, 1, -1
          call lev%Q(m)%copy(lev%Q(lev%nnodes), 2)
          do p = 1, lev%ulevel%sweeper%npieces
            call lev%F(m,p)%copy(lev%F(lev%nnodes,p), 2)
          end do
        end do
      case default
         call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',which)
        stop
    end select

  end subroutine misdcQ_oc_spreadq0

end module pf_mod_misdcQ_oc

