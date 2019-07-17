module pf_mod_parareal
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_pfasst
  use pf_mod_comm
  
  implicit none
  !>  The main PFASST data type which includes pretty much everythingl
  
contains
  !>  Do the parareal algorithm
  subroutine pf_parareal_run(pf, q0, dt, tend, nsteps, qend)
    type(pf_pfasst_t), intent(inout), target   :: pf   !<  The complete PFASST structure
    class(pf_encap_t), intent(in   )           :: q0   !<  The initial condition
    real(pfdp),        intent(in   )           :: dt   !<  The time step for each processor
    real(pfdp),        intent(in   )           :: tend !<  The final time of run
    integer,           intent(in   ), optional :: nsteps  !<  The number of time steps
    class(pf_encap_t), intent(inout), optional :: qend    !<  The computed solution at tend
    !  Local variables

    integer :: nproc  !!  Total number of processors
    integer :: nsteps_loc  !!  local number of time steps
    real(pfdp) :: tend_loc !!  The final time of run


    ! make a local copy of nproc
    nproc = pf%comm%nproc

    !>  Set the number of time steps to do
    !!  The user can either pass in the number of time steps or
    !!  pass in the time step size and length of run
    if (present(nsteps)) then
      nsteps_loc = nsteps
      tend_loc=dble(nsteps_loc*dt)
    else
      nsteps_loc = ceiling(tend/dt)
      !  Do  sanity check on steps
      if (abs(real(nsteps_loc,pfdp)-tend/dt) > dt/100.0) then
        print *,'dt=',dt
        print *,'nsteps=',nsteps_loc
        print *,'tend=',tend
        stop "Invalid nsteps"
      end if
    end if
    pf%state%nsteps = nsteps_loc

    !>  Allocate stuff for holding results
    call pf_initialize_results(pf)

    !  do sanity checks on Nproc
    if (mod(nsteps,nproc) > 0) stop "ERROR: nsteps must be multiple of nproc (pf_parallel.f90)."

    if (present(qend)) then
       call pf_parareal_block_run(pf, q0, dt, nsteps_loc,qend=qend)
    else
       call pf_parareal_block_run(pf, q0, dt,  nsteps_loc)
    end if


    call pf_dump_results(pf)

    !>   deallocate results data
    call pf_destroy_results(pf)

  end subroutine pf_parareal_run

  !>  parareal controller for block mode
  subroutine pf_parareal_block_run(pf, q0, dt, nsteps, qend,flags)
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt
    integer,           intent(in   )           :: nsteps
    class(pf_encap_t), intent(inout), optional :: qend
    integer,           intent(in   ), optional :: flags(:)

    class(pf_level_t), pointer :: lev_p  !!  pointer to the one level we are operating on
    integer                   :: j, k
    integer                   :: nblocks !!  The number of blocks of steps to do
    integer                   :: nproc   !!  The number of processors being used
    integer                   :: level_index_c !!  Coarsest level in V (Lambda)-cycle
    integer                   :: level_max_depth !!  Finest level in V-cycle
    integer::  nsteps_c,nsteps_f  

    call start_timer(pf, TTOTAL)

    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt

    ! set finest level to visit in the following run
    pf%state%finest_level = pf%nlevels

    !  pointer to finest  level to start
    lev_p => pf%levels(pf%state%finest_level)

    !  Stick the initial condition into q0 (will happen on all processors)
    call lev_p%q0%copy(q0, flags=0)


    nproc = pf%comm%nproc
    nblocks = nsteps/nproc

    !  Decide what the coarsest level in the V-cycle is
    level_index_c=1
    if (.not. pf%Vcycle)     level_index_c=pf%state%finest_level

    do k = 1, nblocks   !  Loop over blocks of time steps
       ! print *,'Starting  step=',pf%state%step,'  block k=',k
       ! Each block will consist of
       !  1.  predictor
       !  2.  Vcycle until max iterations, or tolerances met
       !  3.  Move solution to next block

       !  Reset some flags
       !>  When starting a new block, broadcast new initial conditions to all procs
       !>  For initial block, this is done when initial conditions are set

       !> Reset some flags
       pf%state%iter    = -1
       pf%state%itcnt   = 0
       pf%state%mysteps = 0
       pf%state%status  = PF_STATUS_PREDICTOR
       pf%state%pstatus = PF_STATUS_PREDICTOR
       pf%comm%statreq  = -66
       pf%state%pfblock = k
       pf%state%sweep = 1   !  Needed for compatibility of residual storage       


       if (k > 1) then
          if (nproc > 1)  then
             call lev_p%qend%pack(lev_p%send)    !!  Pack away your last solution
             call pf_broadcast(pf, lev_p%send, lev_p%mpibuflen, pf%comm%nproc-1)
             call lev_p%q0%unpack(lev_p%send)    !!  Everyone resets their q0
          else
             call lev_p%q0%copy(lev_p%qend, flags=0)    !!  Just stick qend in q0
          end if

          !>  Update the step and t0 variables for new block
          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt
       end if

       !> Call the predictor to get an initial guess on all levels and all processors
       call pf_parareal_predictor(pf, pf%state%t0, dt, flags)

       !>  Start the parareal iterations
       call start_timer(pf, TITERATION)
       do j = 1, pf%niters

          call call_hooks(pf, -1, PF_PRE_ITERATION)

          pf%state%iter = j

          !  Do a v_cycle
          call pf_parareal_v_cycle(pf, k, pf%state%t0, dt, 1,2)

          !  Check for convergence
          call pf_check_convergence_block(pf, pf%state%finest_level, send_tag=1111*k+j)

          call call_hooks(pf, -1, PF_POST_ITERATION)

          !  If we are converged, exit block
          if (pf%state%status == PF_STATUS_CONVERGED)  exit
       end do  !  Loop over the iteration in this bloc
       call call_hooks(pf, -1, PF_POST_CONVERGENCE)
       call end_timer(pf, TITERATION)
       call call_hooks(pf, -1, PF_POST_STEP)
    end do !  Loop over the blocks

    call end_timer(pf, TTOTAL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       call qend%copy(lev_p%qend, flags=0)
    end if
  end subroutine pf_parareal_block_run
  !>  The parareal predictor does a serial integration on the coarse level followed
  !>  by a fine integration if there is a fine level
  
  subroutine pf_parareal_predictor(pf, t0, dt, flags)
    type(pf_pfasst_t), intent(inout), target :: pf     !! PFASST main data structure
    real(pfdp),        intent(in   )         :: t0     !! Initial time of this processor
    real(pfdp),        intent(in   )         :: dt     !! time step
    integer,           intent(in   ), optional :: flags(:)  !!  User defined flags

    class(pf_level_t), pointer :: c_lev_p
    class(pf_level_t), pointer :: f_lev_p     !!
    integer                   :: k,n               !!  Loop indices
    integer                   :: nsteps_c,nsteps_f    !!  Number of RK  steps
    integer                   :: level_index     !!  Local variable for looping over levels
    real(pfdp)                :: t0k             !!  Initial time at time step k
    pf%state%iter = 0          

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call start_timer(pf, TPREDICTOR)

    !  This is for one two levels only or one if only RK is done
    c_lev_p => pf%levels(1)
    f_lev_p => pf%levels(pf%state%finest_level)

    if (pf%debug) print*, 'DEBUG --', pf%rank, 'beginning parareal predictor'

    !! Step 1. Getting the initial condition on the coarsest level
    if (pf%state%finest_level > 1) then
       if (pf%q0_style < 2) then  !  Copy coarse
          call c_lev_p%q0%copy(f_lev_p%q0)
       end if
    end if
    level_index = 1

    !!
    !! Step 2. Do coarse level integration, no communication necessary
    nsteps_c= c_lev_p%ulevel%stepper%nsteps  !  Each processor integrates alone
    do n=1,pf%rank+1
       if (n .gt. 1) call c_lev_p%q0%copy(c_lev_p%qend)       
       t0k      = dt*real(n-1,pfdp)
       call c_lev_p%ulevel%stepper%do_n_steps(pf, 1, t0k, c_lev_p%q0,c_lev_p%qend,dt, nsteps_c)
    end do
    ! Save the coarse level value
    call c_lev_p%Q(2)%copy(c_lev_p%qend, flags=0)     


    !!  Step 3:  Return to fine level and Step there
!    if(pf%state%finest_level > 1) then  !  Will do nothing with one level
!       call f_lev_p%q0%copy(c_lev_p%q0, flags=0)       !  Get fine initial condition
!       nsteps_f= f_lev_p%ulevel%stepper%nsteps  !  Each processor integrates alone
!       call f_lev_p%ulevel%stepper%do_n_steps(pf, level_index,pf%state%t0, f_lev_p%q0,f_lev_p%qend, dt, nsteps_f)
!    endif

    call end_timer(pf, TPREDICTOR)

    pf%state%iter   = 1
    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING
    if (pf%debug) print*,  'DEBUG --', pf%rank, 'ending predictor'
    call call_hooks(pf, -1, PF_POST_PREDICTOR)

  end subroutine pf_parareal_predictor

  !> Execute a V-cycle between levels nfine and ncoarse
  subroutine pf_parareal_v_cycle(pf, iteration, t0, dt,level_index_c,level_index_f, flags)


    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer,           intent(in)    :: iteration
    integer,           intent(in)    :: level_index_c  !! Coarsest level of V-cycle
    integer,           intent(in)    :: level_index_f  !! Finest level of V-cycle
    integer, optional, intent(in)    :: flags

    type(pf_level_t), pointer :: f_lev_p, c_lev_p
    integer :: level_index, j,nsteps_f,nsteps_c

    if (pf%nlevels <2) return
    
    !  This is for two levels only
    c_lev_p => pf%levels(1)
    f_lev_p => pf%levels(2)
    nsteps_c= c_lev_p%ulevel%stepper%nsteps
    nsteps_f= f_lev_p%ulevel%stepper%nsteps  


    !  Do fine steps with old initial condition
    if (pf%rank /= 0) then
       call f_lev_p%q0%copy(c_lev_p%q0, flags=0)       !  Get fine initial condition
    end if
    call f_lev_p%ulevel%stepper%do_n_steps(pf, 2,pf%state%t0, f_lev_p%q0,f_lev_p%qend, dt, nsteps_f)
    
    ! Get a new initial condition on coarse
    call pf_recv(pf, c_lev_p, 10000+iteration, .true.)

    !  Step on coarse
    call c_lev_p%ulevel%stepper%do_n_steps(pf, 1,pf%state%t0, c_lev_p%q0,c_lev_p%qend, dt, nsteps_c)

    !  Compute the correction (store in Q(1))
    call c_lev_p%Q(1)%copy(f_lev_p%qend, flags=0)  !  Current 
    call c_lev_p%Q(1)%axpy(-1.0_pfdp,c_lev_p%Q(2)) !       

    ! Save the result of the coarse sweep
    call c_lev_p%Q(2)%copy(c_lev_p%qend, flags=0)     

    ! correct coarse level solution at end (the parareal correction)
    call c_lev_p%qend%axpy(1.0_pfdp,c_lev_p%Q(1))    

    !  Send coarse forward  (nonblocking)
    call pf_send(pf, c_lev_p, 10000+iteration, .false.)
    
    !  Compute the jump in the initial condition
    call f_lev_p%q0_delta%copy(c_lev_p%q0, flags=0)
    call f_lev_p%q0_delta%axpy(-1.0d0,f_lev_p%q0, flags=0)
    f_lev_p%residual=f_lev_p%q0_delta%norm(flags=0)
    call pf_set_resid(pf,2,f_lev_p%residual)

!!$    print *,'after fine steps '
!!$    call f_lev_p%qend%eprint()
!!$    call c_lev_p%qend%eprint()

  end subroutine pf_parareal_v_cycle
  

  !> Subroutine to check if the current processor has converged and
  !> to update the next processor on the status
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor can't be converged yet either
  subroutine pf_check_convergence_block(pf, level_index, send_tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    integer,           intent(in)    :: send_tag  !! identifier for status send and receive

    logical           :: residual_converged, converged


    ! Shortcut for fixed iteration mode
    if (pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)

    !> Check to see if tolerances are met
    call pf_check_residual(pf, level_index, residual_converged)

    !>  Until I hear the previous processor is done, recieve it's status
    if (pf%state%pstatus /= PF_STATUS_CONVERGED) call pf_recv_status(pf, send_tag)

    !>  Check to see if I am converged
    converged = .false.
    if (residual_converged) then
       if (pf%rank == 0) then
          converged = .true.
       else  !  I am not the first processor, so I need to check the previous one
          if (pf%state%pstatus == PF_STATUS_CONVERGED) converged = .true.
       end if
    end if ! (residual_converged)

    !>  For parareal, Proc N is converged after iteration N
    if (pf%rank .lt. pf%state%iter) then
       converged = .true.
    end if
    
    !> Assign status and send it forward
    if (converged) then
       if (pf%state%status == PF_STATUS_ITERATING) then
          !  If I am converged for the first time
          !  then flip my flag and send the last status update
          pf%state%status = PF_STATUS_CONVERGED
          call pf_send_status(pf, send_tag)
       end if
    else
       !  I am not converged, send the news
       pf%state%status = PF_STATUS_ITERATING
       call pf_send_status(pf, send_tag)
    end if

  end subroutine pf_check_convergence_block

  !> Subroutine to test residuals to determine if the current processor has converged.
  subroutine pf_check_residual(pf, level_index, residual_converged)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    logical,           intent(out)   :: residual_converged  !! Return true if residual is below tolerances

    residual_converged = .false.

    ! Check to see if absolute tolerance is met
    if   (pf%levels(level_index)%residual     < pf%abs_res_tol)  then
       if (pf%debug) print*, 'DEBUG --',pf%rank, 'residual tol met',pf%levels(level_index)%residual
       residual_converged = .true.
    end if

  end subroutine pf_check_residual

end module pf_mod_parareal
