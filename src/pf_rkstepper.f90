!!  Runge-Kutta time steppers
!
! This file is part of LIBPFASST.
!
!>  Module to do Runge-Kutta stepping
module pf_mod_rkstepper
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  !>  IMEX or additive or semi-implicit Runge-Kutta stepper  type
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

     subroutine pf_f_eval_p(this,y, t, level_index, f, piece, flags, idx, step)
       import pf_ark_t, pf_encap_t, pfdp
       class(pf_ark_t),   intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y
       real(pfdp),        intent(in   ) :: t
       integer,           intent(in   ) :: level_index
       class(pf_encap_t), intent(inout) :: f
       integer,           intent(in   ) :: piece
       integer, intent(in), optional    :: flags
       integer, intent(in), optional    :: idx       ! index of quadrature node
       integer, intent(in), optional    :: step      ! time step for sequential version
     end subroutine pf_f_eval_p

      subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece, flags)
       import pf_ark_t, pf_encap_t, pfdp
       class(pf_ark_t),   intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y
       real(pfdp),        intent(in   ) :: t
       real(pfdp),        intent(in   ) :: dtq
       class(pf_encap_t), intent(in   ) :: rhs
       integer,           intent(in   ) :: level_index
       class(pf_encap_t), intent(inout) :: f
       integer,           intent(in   ) :: piece
       integer,           intent(in), optional    :: flags

     end subroutine pf_f_comp_p

  end interface

contains

  !> Perform N steps of ark on level level_index and set qend appropriately.
  subroutine ark_do_n_steps(this, pf, level_index, t0, big_dt, nsteps_rk, flags)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_ark_t),   intent(inout)         :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in   )         :: t0           !!  Time at start of time interval
    real(pfdp),        intent(in   )         :: big_dt       !!  Size of time interval to integrate on
    integer,           intent(in)            :: level_index  !!  Level of the index to step on
    integer,           intent(in)            :: nsteps_rk    !!  Number of steps to use
    integer,           intent(in), optional    :: flags  

    class(pf_level_t), pointer               :: lev          !!  Pointer to level level_index
    class(pf_encap_t), allocatable           :: rhs          !!  Accumulated right hand side for implicit solves
    integer                                  :: j, m, n      !!  Loop counters
    real(pfdp)                               :: t            !!  Time
    real(pfdp)                               :: dt           !!  Size of each ark step

    
    lev => pf%levels(level_index)   !! Assign pointer to appropriate level
    

    dt = big_dt/real(nsteps_rk, pfdp)   ! Set the internal time step size based on the number of rk steps

    ! Allocate space for the right-hand side
    call lev%ulevel%factory%create_single(rhs, lev%index,  lev%shape)

    
    do n = 1, nsteps_rk      ! Loop over time steps

       ! Recompute the first explicit function value 
       if (n == 1) then
          call lev%Q(1)%copy(lev%q0)
       else
          call lev%Q(1)%copy(lev%Q(lev%nnodes))
       end if


       ! this assumes that cvec(1) == 0
       if (this%explicit) &
            call this%f_eval(lev%Q(1), t0+dt*(n-1)+dt*this%cvec(1), lev%index, lev%F(1,1),1)
       if (this%implicit) &
            call this%f_eval(lev%Q(1), t0+dt*(n-1)+dt*this%cvec(1), lev%index, lev%F(1,2),2)
     
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
  
  
    !> Perform N steps backward in time of ark on level level_index and set q0 appropriately.
  subroutine ark_do_n_steps_backward(this, pf, level_index, t0, big_dt, nsteps_rk)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_ark_t),   intent(inout)         :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in   )         :: t0           !!  Time at beginning of time interval
    real(pfdp),        intent(in   )         :: big_dt       !!  Size of time interval to integrate on; time interval is [tend, tend-big_dt]
    integer,           intent(in)            :: level_index  !!  Level of the index to step on
    integer,           intent(in)            :: nsteps_rk    !!  Number of steps to use

    class(pf_level_t), pointer               :: lev          !!  Pointer to level level_index
    class(pf_encap_t), allocatable           :: rhs          !!  Accumulated right hand side for implicit solves
    integer                                  :: j, m, n      !!  Loop counters
    real(pfdp)                               :: t            !!  Time
    real(pfdp)                               :: dt           !!  Size of each ark step
    real(pfdp)                               :: tend         !!  Final time of time interval
    
    lev => pf%levels(level_index)   !! Assign pointer to appropriate level
    
    tend = t0 + big_dt
    dt = big_dt/real(nsteps_rk, pfdp)   ! Set the internal time step size based on the number of rk steps

    ! Allocate space for the right-hand side
    call lev%ulevel%factory%create_single(rhs, lev%index,  lev%shape)

    
    do n = 1, nsteps_rk      ! Loop over time steps

       ! Recompute the first explicit function value 
       if (n == 1) then
          call lev%Q(lev%nnodes)%copy(lev%qend,2)
       else
          call lev%Q(lev%nnodes)%copy(lev%Q(1),2)
       end if


       ! this assumes that cvec(1) == 0
       if (this%explicit) &
            call this%f_eval(lev%Q(lev%nnodes), tend-dt*(n-1)-dt*this%cvec(1), lev%index, lev%F(lev%nnodes,1),1,flags=2,idx=lev%nnodes)
       if (this%implicit) &
            call this%f_eval(lev%Q(lev%nnodes), tend-dt*(n-1)-dt*this%cvec(1), lev%index, lev%F(lev%nnodes,2),2,flags=2,idx=lev%nnodes)
     
       ! Loop over stage values
       do m = 1, this%nstages-1  
          
          ! Set current time
          t = tend - dt*(n-1) - dt*this%cvec(m+1)

          ! Initialize the right-hand size for each stage
          call rhs%copy(lev%Q(lev%nnodes),flags=2)

          do j = 1, m

             ! Add explicit rhs
             if (this%explicit) &
                  call rhs%axpy(dt*this%AmatE(m+1,j), lev%F(j,1),flags=2)

             ! Add implicit rhs
             if (this%implicit) &
                  call rhs%axpy(dt*this%AmatI(m+1,j), lev%F(j,2),flags=2)

          end do

          ! Solve the implicit system
          if (this%implicit .and. this%AmatI(m+1,m+1) /= 0) then
             call this%f_comp(lev%Q(m+1), t, dt*this%AmatI(m+1,m+1), rhs, lev%index,lev%F(m+1,2), 2, flags=2)
          else
             call lev%Q(m+1)%copy(rhs,flags=2)
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

  end subroutine ark_do_n_steps_backward

  subroutine ark_initialize(this, lev)
    class(pf_ark_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer    :: nstages
    real(pfdp) :: gamma, delta

    this%explicit = .true.
    this%implicit = .true.

    if (this%order == 2) then 

       !  Ascher-Ruuth-Spiteri

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
       
       this%AmatI(2,2) = gamma
       this%AmatI(3,2) = ONE-gamma
       this%AmatI(3,3) = gamma
       
       this%cvec       = (/ ZERO, gamma, ONE /)
       this%bvecE      = (/ ZERO, ONE-gamma, gamma /)
       this%bvecI      = this%bvecE

    else if (this%order == 3) then
      
       ! Third-order Kennedy-Carpenter

       nstages = 4
      
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

       this%AmatE(2,1) =   1767732205903.0_pfdp  / 2027836641118.0_pfdp
       this%AmatE(3,1) =   5535828885825.0_pfdp  / 10492691773637.0_pfdp
       this%AmatE(3,2) =   788022342437.0_pfdp   / 10882634858940.0_pfdp
       this%AmatE(4,1) =   6485989280629.0_pfdp  / 16251701735622.0_pfdp
       this%AmatE(4,2) = - 4246266847089.0_pfdp  / 9704473918619.0_pfdp
       this%AmatE(4,3) =   10755448449292.0_pfdp / 10357097424841.0_pfdp

       this%AmatI(2,1) =   1767732205903.0_pfdp  / 4055673282236.0_pfdp
       this%AmatI(2,2) =   1767732205903.0_pfdp  / 4055673282236.0_pfdp
       this%AmatI(3,1) =   2746238789719.0_pfdp  / 10658868560708.0_pfdp
       this%AmatI(3,2) = - 640167445237.0_pfdp   / 6845629431997.0_pfdp
       this%AmatI(3,3) =   1767732205903.0_pfdp  / 4055673282236.0_pfdp
       this%AmatI(4,1) =   1471266399579.0_pfdp  / 7840856788654.0_pfdp
       this%AmatI(4,2) = - 4482444167858.0_pfdp  / 7529755066697.0_pfdp
       this%AmatI(4,3) =   11266239266428.0_pfdp / 11593286722821.0_pfdp
       this%AmatI(4,4) =   1767732205903.0_pfdp  / 4055673282236.0_pfdp

       this%cvec       = (/ 0.0_pfdp, 1767732205903.0_pfdp / 2027836641118.0_pfdp, 3.0_pfdp / 5.0_pfdp, 1.0_pfdp /)
       this%bvecE      = (/ 1471266399579.0_pfdp  / 7840856788654.0_pfdp,  - 4482444167858.0_pfdp / 7529755066697.0_pfdp,&
                            11266239266428.0_pfdp / 11593286722821.0_pfdp,   1767732205903.0_pfdp / 4055673282236.0_pfdp /)
       this%bvecI      = this%bvecE   

    else if (this%order == 4) then

       ! Fourth-order Kennedy-Carpenter
       
       nstages = 6 

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

       this%AmatE(2,1) =   0.5_pfdp
       this%AmatE(3,1) =   13861.0_pfdp          / 62500.0_pfdp
       this%AmatE(3,2) =   6889.0_pfdp           / 62500.0_pfdp
       this%AmatE(4,1) = - 116923316275.0_pfdp   / 2393684061468.0_pfdp
       this%AmatE(4,2) = - 2731218467317.0_pfdp  / 15368042101831.0_pfdp
       this%AmatE(4,3) =   9408046702089.0_pfdp  / 11113171139209.0_pfdp
       this%AmatE(5,1) = - 451086348788.0_pfdp   / 2902428689909.0_pfdp
       this%AmatE(5,2) = - 2682348792572.0_pfdp  / 7519795681897.0_pfdp
       this%AmatE(5,3) =   12662868775082.0_pfdp / 11960479115383.0_pfdp
       this%AmatE(5,4) =   3355817975965.0_pfdp  / 11060851509271.0_pfdp
       this%AmatE(6,1) =   647845179188.0_pfdp   / 3216320057751.0_pfdp
       this%AmatE(6,2) =   73281519250.0_pfdp    / 8382639484533.0_pfdp
       this%AmatE(6,3) =   552539513391.0_pfdp   / 3454668386233.0_pfdp
       this%AmatE(6,4) =   3354512671639.0_pfdp  / 8306763924573.0_pfdp 
       this%AmatE(6,5) =   4040.0_pfdp           / 17871.0_pfdp
       
       this%AmatI(2,1) =   0.25_pfdp
       this%AmatI(2,2) =   0.25_pfdp
       this%AmatI(3,1) =   8611.0_pfdp        / 62500.0_pfdp
       this%AmatI(3,2) = - 1743.0_pfdp        / 31250.0_pfdp
       this%AmatI(3,3) =   0.25_pfdp
       this%AmatI(4,1) =   5012029.0_pfdp     / 34652500.0_pfdp
       this%AmatI(4,2) = - 654441.0_pfdp      / 2922500.0_pfdp
       this%AmatI(4,3) =   174375.0_pfdp      / 388108.0_pfdp
       this%AmatI(4,4) =   0.25_pfdp
       this%AmatI(5,1) =   15267082809.0_pfdp / 155376265600.0_pfdp
       this%AmatI(5,2) = - 71443401.0_pfdp    / 120774400.0_pfdp
       this%AmatI(5,3) =   730878875.0_pfdp   / 902184768.0_pfdp
       this%AmatI(5,4) =   2285395.0_pfdp     / 8070912.0_pfdp
       this%AmatI(5,5) =   0.25_pfdp     
       this%AmatI(6,1) =   82889.0_pfdp       / 524892.0_pfdp
       this%AmatI(6,2) =   0.0_pfdp
       this%AmatI(6,3) =   15625.0_pfdp       / 83664.0_pfdp
       this%AmatI(6,4) =   69875.0_pfdp       / 102672.0_pfdp
       this%AmatI(6,5) = - 2260.0_pfdp        / 8211.0_pfdp
       this%AmatI(6,6) =   0.25_pfdp
      
       this%cvec       = (/ 0.0_pfdp,                     0.5_pfdp,                  83.0_pfdp / 250.0_pfdp,       &
                            31.0_pfdp / 50.0_pfdp,        17.0_pfdp / 20.0_pfdp,     1.0_pfdp /)
       this%bvecE      = (/ 82889.0_pfdp / 524892.0_pfdp, 0.0_pfdp,                  15625.0_pfdp /  83664.0_pfdp, &
                            69875.0_pfdp / 102672.0_pfdp, - 2260.0_pfdp / 8211.0_pfdp, 0.25_pfdp /)
       this%bvecI      = this%bvecE   

    else
       stop "ark_initialize: This RK order is not supported"
    end if

    if (lev%nnodes < this%nstages + 1)  &
         stop "ark_initialize: With RK, lev%nnodes should be equal to rkstepper%nstages + 1"

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
