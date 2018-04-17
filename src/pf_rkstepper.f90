module pf_mod_rkstepper
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  type, extends(pf_stepper_t), abstract :: pf_ark_t
     real(pfdp), allocatable :: AmatI(:,:)
     real(pfdp), allocatable :: AmatE(:,:)
     real(pfdp), allocatable :: cvec(:)
     real(pfdp), allocatable :: bvecI(:)
     real(pfdp), allocatable :: bvecE(:)
     real(pfdp), allocatable :: QtilI(:,:)
     logical                 :: explicit = .true.
     logical                 :: implicit = .true.
     integer                 :: nstages
   contains
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: do_n_steps  => ark_do_n_steps
     procedure :: initialize  => ark_initialize
     procedure :: destroy     => ark_destroy
  end type pf_ark_t

  interface

     subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
       import pf_ark_t, pf_encap_t, pfdp
       class(pf_ark_t),   intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y
       real(pfdp),        intent(in   ) :: t
       integer,           intent(in   ) :: level_index
       class(pf_encap_t), intent(inout) :: f
       integer,           intent(in   ) :: piece
     end subroutine pf_f_eval_p

      subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       import pf_ark_t, pf_encap_t, pfdp
       class(pf_ark_t),   intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y
       real(pfdp),        intent(in   ) :: t
       real(pfdp),        intent(in   ) :: dtq
       class(pf_encap_t), intent(in   ) :: rhs
       integer,           intent(in   ) :: level_index
       class(pf_encap_t), intent(inout) :: f
       integer,           intent(in   ) :: piece
     end subroutine pf_f_comp_p

  end interface

contains

  !> Perform N step of ark on level level_index and set qend appropriately.
  subroutine ark_do_n_steps(this, pf, level_index, t0, big_dt, nsteps_rk)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_ark_t),   intent(inout)         :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in   )         :: t0           !  Time at start of time interval
    real(pfdp),        intent(in   )         :: big_dt       !  Size of time interval to integrato on
    integer,           intent(in)            :: level_index  !  Level of the index to step on
    integer,           intent(in)            :: nsteps_rk    !  Number of steps to use

    class(pf_level_t), pointer               :: lev          !  Pointer to level level_index
    class(pf_encap_t), allocatable           :: rhs          !  Accumulated right hand side for implicit solves
    integer                                  :: j, m, n   !  Loop counters
    real(pfdp)                               :: t, dt        !  Size of each ark step

    ! Assign pointer to appropriate level
    lev => pf%levels(level_index)   
    
    ! Set the internal time step size based on the number of rk steps
    dt = big_dt/real(nsteps_rk, pfdp)

    ! Allocate space for the right-hand side
    call lev%ulevel%factory%create_single(rhs, lev%index, lev%shape)

    ! Loop over time steps
    do n = 1, nsteps_rk      

       ! Recompute the first explicit function value 
       if (n == 1) then
          call lev%Q(1)%copy(lev%q0)
       else
          call lev%Q(1)%copy(lev%Q(lev%nnodes))
       end if

       if (this%explicit) &
            call this%f_eval(lev%Q(1), t0+dt*(n-1)+dt*this%cvec(1), lev%index, lev%F(1,1),1)
       if (this%implicit) &
            call lev%F(1,2)%setval(0.0_pfdp)
     
       ! Loop over stage values
       do m = 1, this%nstages-1  
          
          ! Set current time
          t = t0 + dt*(n-1) + dt*this%cvec(m+1)

          ! Initialize the right-hand size for each stage
          call rhs%copy(lev%Q(1))

          do j = 1, m

             ! Add explicit rhs
             if (this%explicit) &
                  call rhs%axpy(dt*this%AmatE(m+1,j), lev%F(j,1))

             ! Add implicit rhs
             if (this%implicit) &
                  call rhs%axpy(dt*this%AmatI(m+1,j), lev%F(j,2))

          end do

          ! Solve the implicit system
          if (this%implicit .and. this%AmatI(m+1,m+1) /= 0) then
             call this%f_comp(lev%Q(m+1), t, dt*this%AmatI(m+1,m+1), rhs, lev%index,lev%F(m+1,2), 2)
          else
             call lev%Q(m+1)%copy(rhs)
          end if
                    
          ! Reevaluate explicit rhs with the new solution
          if (this%explicit) &
               call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1), 1)
          
       end do  ! End loop over stage values
       
       ! Compute final value using quadrature rule
       call lev%Q(lev%nnodes)%copy(lev%Q(1))

       ! Loop over stage values one more time
       do j = 1, this%nstages

          ! Add explicit terms
          if (this%explicit) &
               call lev%Q(lev%nnodes)%axpy(dt*this%bvecE(j), lev%F(j,1))

          ! Add implicit terms
          if (this%implicit) &
               call lev%Q(lev%nnodes)%axpy(dt*this%bvecI(j), lev%F(j,2))

       end do ! End loop over stage values

    end do ! End Loop over time steps
    
    ! Assign final value to end of time step
    call lev%qend%copy(lev%Q(lev%nnodes))

  end subroutine ark_do_n_steps
  

  subroutine ark_initialize(this, lev)
    class(pf_ark_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer    :: nstages
    real(pfdp) :: gamma, delta

    !  Ascher-Ruuth-Spiteri
    this%explicit = .true.
    this%implicit = .true.

    nstages = 3

    this%nstages = nstages
    allocate(this%AmatE(nstages,nstages))  !  Explicit Butcher matrix
    allocate(this%AmatI(nstages,nstages))  !  Implicit Butcher matrix
    allocate(this%cvec(nstages))           !  stage times
    allocate(this%bvecE(nstages))          !  quadrature weights on explicit
    allocate(this%bvecI(nstages))          !  quadrature weights on implicit

    this%AmatE = 0.0_pfdp
    this%AmatI = 0.0_pfdp
    this%bvecE = 0.0_pfdp
    this%bvecI = 0.0_pfdp
    this%cvec  = 0.0_pfdp
    
    gamma           = (TWO - sqrt(TWO))/TWO
    delta           = -TWO*sqrt(TWO)/THREE

    this%AmatE(2,1) = gamma
    this%AmatE(3,1) = delta
    this%AmatE(3,2) = ONE-delta

    print *, this%AmatE

    this%AmatI(2,2) = gamma
    this%AmatI(3,2) = ONE-gamma
    this%AmatI(3,3) = gamma

    print *, this%AmatI

    this%cvec       = (/ ZERO, gamma, ONE /)
    this%bvecE      = (/ ZERO, ONE-gamma, gamma /)
    this%bvecI      = this%bvecE

  end subroutine ark_initialize

  subroutine ark_destroy(this, lev)
    class(pf_ark_t),   intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%AmatE)
    deallocate(this%AmatI)
    deallocate(this%bvecE)
    deallocate(this%bvecI)
    deallocate(this%cvec)
  end subroutine ark_destroy

end module pf_mod_rkstepper
