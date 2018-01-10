module pf_mod_rkstepper
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  type, extends(pf_stepper_t), abstract :: pf_ARK_t
     real(pfdp), allocatable :: AmatI(:,:)
     real(pfdp), allocatable :: AmatE(:,:)
     real(pfdp), allocatable :: cvec(:)
     real(pfdp), allocatable :: bvecI(:)
     real(pfdp), allocatable :: bvecE(:)
     real(pfdp), allocatable :: QtilI(:,:)
     logical                 :: explicit = .true.
     logical                 :: implicit = .true.
     integer                 :: Nstages
   contains
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: do_N_steps    => ARK_Nsteps
     procedure :: initialize => ARK_initialize
     procedure :: destroy   => ARK_destroy
  end type pf_ARK_t

  interface
     subroutine pf_f_eval_p(this,y, t, level, f, piece)
       import pf_ARK_t, pf_encap_t, pfdp
       class(pf_ARK_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y
       real(pfdp),        intent(in   ) :: t
       integer,    intent(in   ) :: level
       class(pf_encap_t), intent(inout) :: f
       integer,    intent(in   ) :: piece
     end subroutine pf_f_eval_p
      subroutine pf_f_comp_p(this,y, t, dt, rhs, level, f, piece)
       import pf_ARK_t, pf_encap_t, pfdp
       class(pf_ARK_t),  intent(inout) :: this
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
    !> Perform N step of ARK on level lev_index and set qend appropriately.
    subroutine ARK_Nsteps(this, pf, level_index, t0, big_dt,nsteps)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_ARK_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    real(pfdp),        intent(in   ) :: t0             !  Time at start of time interval
    real(pfdp),        intent(in   ) :: big_dt         !  Size of time interval to integrato on
    integer,             intent(in)    :: level_index  !  Level of the index to step on
    integer,             intent(in)    :: nsteps       !  Number of steps to use

    class(pf_level_t), pointer :: lev

    integer     :: i, j,m,n  !  Loop counters
    real(pfdp)  :: t,dt    !  size of each ARK step

    class(pf_encap_t), allocatable :: rhs   !  The accumulated right hand side for implicit solves

    lev => pf%levels(level_index)   !  Assign pointer to appropriate level

    dt = big_dt/real(nsteps, pfdp)


    do n=1,nsteps      !  Loop over time steps
       do m=1,this%Nstages  !  Loop over stage values
          t = t + dt*this%cvec(m)

          call rhs%setval(0.0_pfdp)
          do j = 1, m
             if (this%explicit) &
                  call rhs%axpy(dt*this%AmatE(m,j), lev%F(n,j))
             if (this%implicit) &
                  call rhs%axpy(dt*this%AmatI(m,j), lev%F(n,j))

             call rhs%axpy(1.0_pfdp, lev%Q(1))
             if (this%implicit) then
                call this%f_comp(lev%Q(m+1), t, dt*this%AmatI(m,m+1), rhs, lev%index,lev%F(m+1,2),2)
             else
                call lev%Q(m+1)%copy(rhs)
             end if
             if (this%explicit) &
                  call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1),1)

          end do
          

          
       end do          !  Loop over stage values
       

       !  Compute final value using quadrature rule
       call lev%Q(lev%nnodes)%setval(0.0_pfdp)
       do j = 1, this%Nstages
          if (this%explicit) &
               call lev%Q(lev%nnodes)%axpy(dt*this%bvecE(j), lev%F(j,1))
          if (this%implicit) &
               call lev%Q(lev%nnodes)%axpy(dt*this%bvecI(j), lev%F(j,2))
       end do
    end do             !  Loop over time steps: n
    
    !  Assign final value to end of time step
    call lev%qend%copy(lev%Q(lev%nnodes))
  end subroutine ARK_Nsteps
  

  subroutine ARK_initialize(this, lev)
    class(pf_ARK_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer  :: Nstages
    real(pfdp) :: gamma, delta, TWO

    !  Ascher-Ruuth-Spiteri
    this%explicit = .TRUE.
    this%implicit = .TRUE.
    Nstages = 3

    this%Nstages = Nstages
    allocate(this%AmatE(Nstages,Nstages))  !  Explicit Butcher matrix
    allocate(this%AmatI(Nstages,Nstages))  !  Implicit Butcher matrix
    allocate(this%cvec(Nstages) ) !  stage times
    allocate(this%bvecE(Nstages))  !  quadrature weights on explicit
    allocate(this%bvecI(Nstages))  !  quadrature weights on implicit

    this%AmatE = 0.0_pfdp
    this%AmatI = 0.0_pfdp
    this%bvecE = 0.0_pfdp
    this%bvecI = 0.0_pfdp
    this%cvec=0.0_pfdp


    gamma = (TWO - sqrt(TWO))/TWO
    delta = -2*sqrt(TWO)/THREE
    this%AmatE(2,1) = gamma
    this%AmatE(3,1) = delta
    this%AmatE(3,2) = ONE-delta

    this%AmatI(2,2) = gamma
    this%AmatI(3,2) = ONE-gamma
    this%AmatI(3,3) = gamma

    this%cvec =(/ZERO, gamma, ONE/)
    this%bvecE =(/ZERO, ONE-gamma, gamma/)
    this%bvecI =this%bvecE
  end subroutine ARK_initialize

  subroutine ARK_destroy(this, lev)
    class(pf_ARK_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%AmatE)
    deallocate(this%AmatI)
    deallocate(this%bvecE)
    deallocate(this%bvecI)
    deallocate(this%cvec)
  end subroutine ARK_destroy

end module pf_mod_rkstepper
