!!  Runge-Kutta time steppers
!
! This file is part of LIBPFASST.
!
!>  Module to do IMEX Euler stepping
module pf_mod_eulerstepper
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  !>  IMEX or additive or semi-implicit Runge-Kutta stepper  type
  type, extends(pf_stepper_t), abstract :: pf_eulerstepper_t
     logical                 :: explicit = .true.
     logical                 :: implicit = .true.
!      real(pfdp)              :: cvec(2)
   contains
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: do_n_steps  => euler_do_n_steps
     procedure :: initialize  => euler_initialize
     procedure :: destroy     => euler_destroy
  end type pf_eulerstepper_t

  interface

     subroutine pf_f_eval_p(this,y, t, level_index, f,  piece, flags, idx, step)
       import pf_eulerstepper_t, pf_encap_t, pfdp
       class(pf_eulerstepper_t),   intent(inout) :: this
       class(pf_encap_t),   intent(in   ) :: y
       real(pfdp),          intent(in   ) :: t
       integer,             intent(in   ) :: level_index
       class(pf_encap_t),   intent(inout) :: f
       integer,             intent(in   ) :: piece
       integer, intent(in)              :: flags
       integer, intent(in), optional    :: idx       ! index of quadrature node
       integer, intent(in), optional    :: step      ! time step for sequential version
     end subroutine pf_f_eval_p

      subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece, flags)
       import pf_eulerstepper_t, pf_encap_t, pfdp
       class(pf_eulerstepper_t),   intent(inout) :: this
       class(pf_encap_t),   intent(inout) :: y
       real(pfdp),          intent(in   ) :: t
       real(pfdp),          intent(in   ) :: dtq
       class(pf_encap_t),   intent(in   ) :: rhs
       integer,             intent(in   ) :: level_index
       class(pf_encap_t),   intent(inout) :: f
       integer,             intent(in   ) :: piece
      integer,              intent(in)    :: flags

     end subroutine pf_f_comp_p

  end interface

contains

  !> Perform N fwd or bwd steps of IMEX Euler on level level_index and set qend or q0 appropriately.
  subroutine euler_do_n_steps(this, pf, level_index, t0, q0, qend, big_dt, nsteps_rk, &
                              state, adjoint, flags)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_eulerstepper_t), intent(inout)           :: this
    type(pf_pfasst_t), intent(inout), target   :: pf
    real(pfdp),        intent(in   )           :: t0           !!  Time at start of time interval
    class(pf_encap_t), intent(inout)         :: q0           !!  Starting value
    class(pf_encap_t), intent(inout)         :: qend !!  Final value
    real(pfdp),        intent(in   )           :: big_dt       !!  Length of time interval to integrate on
    integer,           intent(in)              :: level_index  !!  Level of the index to step on
    integer,           intent(in)              :: nsteps_rk       !!  Number of steps to use
!     real(pfdp),        intent(inout), optional :: state(:,:,:), adjoint(:,:,:)  !! store the solution at all time steps (size has to be nsteps times nvars)
    real(pfdp),        intent(inout), optional :: state(:,:), adjoint(:,:)  !! store the solution at all time steps (size has to be nsteps times nvars)
    integer,           intent(in), optional    :: flags        !!  which component to compute
    
    class(pf_level_t), pointer               :: lev          !!  Pointer to level level_index
    class(pf_encap_t), allocatable           :: rhs          !!  Accumulated right hand side for implicit solves
    integer                                  :: j, m, n      !!  Loop counters
    real(pfdp)                               :: t, tend      !!  Time
    real(pfdp)                               :: dt           !!  Size of each ark step

    integer                                  :: which        !!  which component to compute
    integer                                  :: mystep       !!  for the time transform in backward equation
    integer :: nstart, nend, inc
    
    which = 1;
    if(present(flags)) which = flags ! which component to compute
    if (which == 0)    which = 1     ! mixed integration not supported
    
    lev => pf%levels(level_index)   !! Assign pointer to appropriate level
    

    tend = t0 + big_dt
    dt = big_dt/real(nsteps_rk, pfdp)   ! Set the internal time step size based on the number of rk steps
    if (which == 2) dt = -dt ! backward stepping
    
    ! Allocate space for the right-hand side
    call lev%ulevel%factory%create_single(rhs, lev%index,  lev%shape)

    if (which == 2) then
      nstart = nsteps_rk+1
      nend = 2
      inc = -1
    else
      nstart = 1
      nend = nsteps_rk
      inc=1
    end if
   
!print *, which, nstart, nend, inc
 
    if (which .eq. 1) then ! forward solve
      call lev%Q(1)%copy(lev%q0, 1)
!       call lev%Q(1)%pack(state(1,:), 1)      
    else
      call lev%Q(1)%copy(lev%qend, which)
!       call lev%Q(1)%unpack(state(nsteps_rk+1,:), 1)
!       call lev%Q(1)%pack(adjoint(nsteps_rk+1,:), 2)
      !               call solution(nsteps_rk)%copy(lev%Q(1),which) ! store terminal solution
    end if
    
    do n = nstart, nend, inc !1, nsteps_rk      ! Loop over time steps

       if(which==2) then
!          call lev%Q(1)%unpack(state(n,lev%nnodes,:), 1)
!          call lev%Q(1)%pack(adjoint(n,lev%nnodes,:), 2)
         call lev%Q(1)%unpack(state(n,:), 1)
         call lev%Q(1)%pack(adjoint(n,:), 2)
       else
!          call lev%Q(1)%pack(state(n,1,:), 1)  
         call lev%Q(1)%pack(state(n,:), 1)  
       end if
       ! Recompute the first explicit function value 
        
       ! t is the _end_ of the current local time step
       if (which .eq. 2) then ! backward solve
          t = tend + dt*n    ! dt is negative here  
          mystep = n
          ! here we need to unpack the state solution into Q(1)
!           call lev%Q(n)%copy(solution(mystep), 1)
       else
          t = t0 + dt*n
          mystep = n
       end if
       
       if (this%explicit) &
            call this%f_eval(lev%Q(1), t-dt, lev%index, lev%F(1,1), 1, &
                            flags=which, idx=1, step=mystep )
!        if (this%implicit) &
!             call this%f_eval(lev%Q(1), t, lev%index, lev%F(1,2),2, &
!                              flags=which, idx=1, step=mystep )
     
       ! set up rhs
       call rhs%copy(lev%Q(1), which)
       if (this%explicit) &
          call rhs%axpy(inc*dt, lev%F(1,1), which)
          !print *, n+inc, n
      ! solve
      !print *, t, inc*dt, n+inc
      if (this%implicit) then
        call this%f_comp(lev%Q(2), t, inc*dt, rhs, lev%index, lev%F(2,2), 2, which)
      else
        call lev%Q(2)%copy(rhs,which)
      end if
      
      call lev%Q(1)%copy(lev%Q(2),which)
      call lev%F(1,2)%copy(lev%F(2,2),which)
 
      ! store solution
       if (which .eq. 2) then ! backward solve
         call lev%Q(1)%pack(adjoint(n+inc,:),2) ! store solution value
         call lev%Q(1)%unpack(state(n+inc,:), 1)
!          call lev%Q(1)%pack(adjoint(n,1,:),2) ! store solution value
!          call lev%Q(1)%unpack(state(n,1,:), 1)
       else
!          call lev%Q(1)%pack(state(n,lev%nnodes,:), 1)
         call lev%Q(1)%pack(state(n+inc,:), 1)
       end if

    end do ! End Loop over time steps
    
    ! Assign final value to end of time step
    if (which .eq. 2) then ! backward solve
      call lev%q0%copy(lev%Q(1), which)
    else
      call lev%qend%copy(lev%Q(1), which)
    end if

  end subroutine euler_do_n_steps
  
  
  

  subroutine euler_initialize(this, pf, level_index)
    class(pf_eulerstepper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize

    this%explicit = .true.
    this%implicit = .true.
    allocate(this%cvec(2))
    this%cvec = (/0.0_pfdp, 1.0_pfdp/)
  end subroutine euler_initialize

  subroutine euler_destroy(this, pf, level_index)
    class(pf_eulerstepper_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
  end subroutine euler_destroy

end module pf_mod_eulerstepper
