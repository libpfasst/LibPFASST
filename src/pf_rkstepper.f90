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
  type, extends(pf_stepper_t), abstract :: pf_ark_stepper_t
     real(pfdp), allocatable :: AmatI(:,:)
     real(pfdp), allocatable :: AmatE(:,:)
     real(pfdp), allocatable :: cvec(:)
     real(pfdp), allocatable :: bvecI(:)
     real(pfdp), allocatable :: bvecE(:)
     real(pfdp), allocatable :: QtilI(:,:)
     logical                 :: explicit = .true.
     logical                 :: implicit = .true.
     integer                 :: nstages

     !  Local storage (not optimal)
     class(pf_encap_t), allocatable           :: rhs          !!  Accumulated right hand side for implicit solves
     class(pf_encap_t), allocatable           :: ytemp        !!  Temp for y
     class(pf_encap_t), allocatable           :: y0           !!  Local y0
     class(pf_encap_t), pointer               :: F(:,:)       !!  Pointer to F

   contains
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: do_n_steps  => ark_do_n_steps
     procedure :: initialize  => ark_initialize
     procedure :: destroy     => ark_destroy
     procedure :: ark_initialize 
     procedure :: ark_destroy    
  end type pf_ark_stepper_t

  interface

     subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
       import pf_ark_stepper_t, pf_encap_t, pfdp
       class(pf_ark_stepper_t),   intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y
       real(pfdp),        intent(in   ) :: t
       integer,           intent(in   ) :: level_index
       class(pf_encap_t), intent(inout) :: f
       integer,           intent(in   ) :: piece
     end subroutine pf_f_eval_p

      subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       import pf_ark_stepper_t, pf_encap_t, pfdp
       class(pf_ark_stepper_t),   intent(inout) :: this
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

  !> Perform N steps of ark on level level_index and set yend appropriately.  Note that input y0 is not disturbed
  subroutine ark_do_n_steps(this, pf, level_index, t0, y0,yend,big_dt, nsteps_rk)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_ark_stepper_t),   intent(inout)         :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in   )         :: t0           !!  Time at start of time interval
    class(pf_encap_t), intent(in   )         :: y0           !!  Starting value
    class(pf_encap_t), intent(inout)         :: yend         !!  Final value
    real(pfdp),        intent(in   )         :: big_dt       !!  Size of time interval to integrato on
    integer,           intent(in)            :: level_index  !!  Level of the index to step on
    integer,           intent(in)            :: nsteps_rk    !!  Number of steps to use

    class(pf_level_t), pointer               :: lev          !!  Pointer to level level_index
    integer                                  :: j, m, n      !!  Loop counters
    real(pfdp)                               :: tn           !!  Time at beginning of RKstep
    real(pfdp)                               :: tc           !!  Time at  RK stage    
    real(pfdp)                               :: dt           !!  Size of each ark step

    
    lev => pf%levels(level_index)   !! Assign pointer to appropriate level

    dt = big_dt/real(nsteps_rk, pfdp)   ! Set the internal time step size based on the number of rk steps

    call this%y0%copy(y0)
    
    do n = 1, nsteps_rk      ! Loop over time steps
       tn=t0+dt*real(n-1,pfdp)
       ! Reset initial condition
       if (n > 1) call this%y0%copy(this%ytemp)


       ! this assumes that cvec(1) == 0
       if (this%explicit) &
            call this%f_eval(this%y0, tn+dt*this%cvec(1), level_index, this%F(1,1),1)
       if (this%implicit) &
            call this%f_eval(this%y0, tn+dt*this%cvec(1), level_index, this%F(1,2),2)
     
       ! Loop over stage values
       do m = 1, this%nstages-1  
          
          ! Set current time
          tc = tn + dt*this%cvec(m+1)

          ! Initialize the right-hand size for each stage
          call this%rhs%copy(this%y0)

          do j = 1, m

             ! Add explicit rhs
             if (this%explicit) &
                  call this%rhs%axpy(dt*this%AmatE(m+1,j), this%F(j,1))

             ! Add implicit rhs
             if (this%implicit) &
                  call this%rhs%axpy(dt*this%AmatI(m+1,j), this%F(j,2))

          end do

          ! Solve the implicit system
          if (this%implicit .and. this%AmatI(m+1,m+1) /= 0) then
             call this%f_comp(this%ytemp, tc, dt*this%AmatI(m+1,m+1), this%rhs, level_index,this%F(m+1,2), 2)
          else
             call this%ytemp%copy(this%rhs)
          end if
                    
          ! Reevaluate explicit rhs with the new solution
          if (this%explicit) &
               call this%f_eval(this%ytemp, tc, level_index, this%F(m+1,1), 1)
          
       end do  ! End loop over stage values
       
       ! Compute final value using quadrature rule
       call this%ytemp%copy(this%y0)

       ! Loop over stage values one more time
       do j = 1, this%nstages

          ! Add explicit terms
          if (this%explicit) &
               call this%ytemp%axpy(dt*this%bvecE(j), this%F(j,1))

          ! Add implicit terms
          if (this%implicit) &
               call this%ytemp%axpy(dt*this%bvecI(j), this%F(j,2))

       end do ! End loop over stage values

    end do ! End Loop over time steps
    call yend%copy(this%ytemp)
  end subroutine ark_do_n_steps
  

  subroutine ark_initialize(this, pf, level_index)
    class(pf_ark_stepper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize

    integer    :: nstages,npieces   !  Local copies for convenience
    integer    :: i,ierr   
    real(pfdp) :: gamma, delta
    type(pf_level_t), pointer  :: lev    !  Current level
    lev => pf%levels(level_index)   !  Assign level pointer

    !  The implicit and explicit flags should be set before calling initialize
    npieces = 1
    if (this%implicit .eqv. .true. .and. this%explicit .eqv. .true.) npieces=2

    select case (this%order)
    case (1)  !  Forward-backward Euler
       nstages = 2
    case (2)  !  Ascher-Ruuth-Spiteri
       nstages = 3
    case (3) ! Third-order Kennedy-Carpenter
       nstages = 4
    case (4) ! Fourth-order Kennedy-Carpenter
       nstages = 6 
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'RK order is not supported ,order=',this%order)              
    end select
    this%nstages = nstages
    this%npieces = npieces

    allocate(this%AmatE(nstages,nstages),stat=ierr)  !  Explicit Butcher matrix
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)               
    allocate(this%AmatI(nstages,nstages),stat=ierr)  !  Implicit Butcher matrix
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)               
    allocate(this%cvec(nstages),stat=ierr)           !  stage times
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)               
    allocate(this%bvecE(nstages),stat=ierr)          !  quadrature weights on explicit
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)               
    allocate(this%bvecI(nstages),stat=ierr)          !  quadrature weights on implicit
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)           
    
    this%AmatE = 0.0_pfdp
    this%AmatI = 0.0_pfdp
    this%bvecE = 0.0_pfdp
    this%bvecI = 0.0_pfdp
    this%cvec  = 0.0_pfdp

    select case (this%order)
    case (1)  !  Forward-backward Euler
       this%AmatE(2,1) = ONE
       this%AmatI(2,2) = ONE
       
       this%cvec       = (/ ZERO, ONE /)
       this%bvecE      = (/  ONE, ZERO /)
       this%bvecI      = (/ ZERO, ONE /)
    case (2)  !  Ascher-Ruuth-Spiteri
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

    case (3) ! Third-order Kennedy-Carpenter
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

    case (4) ! Fourth-order Kennedy-Carpenter
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

    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'RK order is not supported ,order=',this%order)              
    end select

    ! Allocate space for local variables
    call lev%ulevel%factory%create_single(this%rhs, level_index,  lev%lev_shape)
    call lev%ulevel%factory%create_single(this%y0, level_index,  lev%lev_shape)
    call lev%ulevel%factory%create_single(this%ytemp, level_index,  lev%lev_shape)            
    call lev%ulevel%factory%create_array(lev%Frkflt, nstages*npieces, level_index,  lev%lev_shape)
    do i = 1, nstages*npieces
       call lev%Frkflt(i)%setval(0.0_pfdp, 0)
    end do
    
    this%F(1:nstages,1:npieces) => lev%Frkflt


  end subroutine ark_initialize

  subroutine ark_destroy(this, pf,level_index)
    class(pf_ark_stepper_t),   intent(inout) :: this
    type(pf_pfasst_t),  target,  intent(inout) :: pf
    integer,              intent(in)    :: level_index 

    type(pf_level_t), pointer  :: lev        !  Current level
    lev => pf%levels(level_index)   !  Assign level pointer

    deallocate(this%AmatE)
    deallocate(this%AmatI)
    deallocate(this%bvecE)
    deallocate(this%bvecI)
    deallocate(this%cvec)
    call lev%ulevel%factory%destroy_single(this%y0)
    call lev%ulevel%factory%destroy_single(this%rhs)
    call lev%ulevel%factory%destroy_single(this%ytemp)
    call lev%ulevel%factory%destroy_array(lev%Frkflt)
  end subroutine ark_destroy

end module pf_mod_rkstepper
