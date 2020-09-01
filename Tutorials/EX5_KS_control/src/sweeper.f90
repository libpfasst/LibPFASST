!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 1-D K-S equation.
!>     u_t + v*u_x = nu*u_xx
module my_sweeper
  use pf_mod_dtype
  use pf_mod_ndarray_oc
  use pf_mod_imexQ_oc
  use pf_mod_fftpackage
  use pf_mod_restrict
    
  use probin
  
  implicit none

  !>  extend the imex_oc sweeper type with stuff we need to compute rhs
  type, extends(pf_imexQ_oc_t) :: my_sweeper_t
     integer ::     nx   !  Grid size

     !>  FFT and Spectral derivatives
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: lap(:) ! Explicit spectral operator
     complex(pfdp), allocatable :: ddx(:) ! Implicit spectral operator
     
     real(pfdp), allocatable :: ydesired(:,:,:) !(time step, quadrature node, nx)
     integer                 :: nsteps_per_rank, nproc, myrank

   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves
     procedure :: initialize  !  Bypasses base sweeper initialize
     procedure :: destroy     !  Bypasses base sweeper destroy

  end type my_sweeper_t

contains

  !>  Helper function to return sweeper pointer
  function as_my_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(my_sweeper_t), pointer :: r
    select type(sweeper)
    type is (my_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_my_sweeper


  !>  Routine to initialize sweeper (bypasses imex sweeper initialize)
  subroutine initialize(this, pf, level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index

    integer     :: nx

    !  Call the imex sweeper initialize
    call this%imexQ_oc_initialize(pf,level_index)    

    this%implicit=.TRUE.
    this%explicit=.TRUE.
  
    nx=pf%levels(level_index)%lev_shape(1)  !  local convenient grid size

    !>  Set up the FFT 
    allocate(this%fft_tool)
    call this%fft_tool%fft_setup([nx],1)

    !>  Define spectral derivatitive operators
    allocate(this%lap(nx))
    allocate(this%ddx(nx))
    call this%fft_tool%make_lap(this%lap)
    call this%fft_tool%make_deriv(this%ddx)
    
        
  end subroutine initialize
  

  !>  Destroy sweeper (bypasses base sweeper destroy)
  subroutine destroy(this,pf,level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),  target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    !>  Call base sweeper destroy
    call this%imexQ_oc_destroy(pf,level_index)

    !> Nuke the FFT operators 
    deallocate(this%lap)
    deallocate(this%ddx)

    !>  Free up FFT stuff
    call this%fft_tool%fft_destroy()
    deallocate(this%fft_tool)
    
    !> Free up optimization data
    if(allocated(this%ydesired)) &
       deallocate(this%ydesired)

  end subroutine destroy
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece, flags, idx, step)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece  !  Which piece to solve for
    integer,             intent(in   ) :: flags
    integer,             intent(in   ), optional :: idx    ! index of quadrature node
    integer,             intent(in   ), optional :: step   ! time step for sequential version
    
    real(pfdp),       pointer     :: yvec(:), fvec(:), pvec(:) ! state, rhs, adjoint
    real(pfdp),       allocatable :: tmp(:)
    type(pf_fft_t),   pointer     :: fft
    
    integer :: thisstep, mystep, nx
    

    ! which time step are we on (to access the correct data of the optimization problem)
    thisstep = 1
    if(present(step)) thisstep = step

    ! the step passed here is pf%state%step+1 which is the actual time step to be computed
    ! the step that the respective processor perform is given as (note that thisstep = pf%state%step+1!):
    mystep = (thisstep - 1 - this%myrank)/this%nproc + 1

    fft => this%fft_tool
    
    nx = size(this%ydesired,3)
    allocate(tmp(nx))
        
    ! Apply spectral operators using the FFT convolution function
    select case (piece)
    case (1)  ! Explicit piece
       select case (flags) ! are we solving for state or adjoint?
       case(1) ! State
          !  Grab the arrays from the encap
          yvec => get_array1d_oc(y,1) ! ( ,1) means get state component
          fvec => get_array1d_oc(f,1)
        
          ! explicit piece: -y y_x
          tmp = -0.5_pfdp*yvec*yvec
          call fft%conv(tmp, this%ddx, fvec)
       
       case(2) ! Adjoint
          !  Grab the arrays from the encap
          yvec => get_array1d_oc(y,1) ! state sol required for rhs (y-y_desired)
          fvec => get_array1d_oc(f,2)
          pvec => get_array1d_oc(y,2)
        
        ! explicit piece: +y y_x
        tmp = 0.5_pfdp*yvec*yvec
        call fft%conv(tmp, this%ddx, fvec)
        
        ! source term from distributed tracking objective
        fvec = fvec + (yvec-this%ydesired(mystep, idx, :))
       
       case DEFAULT
        print *, "ERROR in f_eval: only 1, 2 allowed as flags", flags
        call exit(1)
       end select
       
    case (2)  ! Implicit piece
       select case (flags)
       case (1)
          yvec => get_array1d_oc(y,1)
          fvec => get_array1d_oc(f,1)
          call fft%conv(yvec,-this%lap-this%lap*this%lap,fvec)
       case (2)
          pvec => get_array1d_oc(y,2)
          fvec => get_array1d_oc(f,2)
          call fft%conv(pvec,-this%lap-this%lap*this%lap,fvec)
       case default
          print *, "ERROR in f_eval: only 1, 2 allowed as flags", flags
          call exit(1)
       end select
       
    case DEFAULT
       print *,'Bad case for piece in f_eval ', piece
       call exit(1)
    end select

    deallocate(tmp)
    
  end subroutine f_eval

  ! Solve for y and return f2 also
  !   y-dtq*f(y,t) = rhs
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f, piece, flags)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y       !  The solution we seek
    real(pfdp),          intent(in   ) :: t       !  Equation time of implicit solve
    real(pfdp),          intent(in   ) :: dtq     !  The 
    class(pf_encap_t),   intent(in   ) :: rhs     !  The right hand side of the solve
    integer,             intent(in   ) :: level_index !  Which level this is
    class(pf_encap_t),   intent(inout) :: f       !  The function value
    integer,             intent(in   ) :: piece   !  Designates which piece to solve for (here implicit)
    integer,             intent(in   ) :: flags   !  Which equation we are solving (1 = state, 2 = adjoint)

    real(pfdp),      pointer :: yvec(:), rhsvec(:), fvec(:)
    type(pf_fft_t),     pointer :: fft

    !  Grab the arrays from the encaps
    yvec  => get_array1d_oc(y,flags)
    rhsvec => get_array1d_oc(rhs,flags)
    fvec => get_array1d_oc(f,flags)

    ! Grab the fft workspace
    fft => this%fft_tool

    if (piece == 2) then
       ! Apply the inverse operator with the FFT convolution
       call fft%conv(rhsvec,1.0_pfdp/(1.0_pfdp - dtq*(-this%lap-this%lap*this%lap)),yvec)

       !  The function is easy to derive
       fvec = (yvec - rhsvec) / dtq
    else
       print *,'Bad piece in f_comp ',piece
       call exit(0)
    end if
  end subroutine f_comp



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial_sol(q0)
    real(pfdp), intent(inout) :: q0(:)
    
    integer :: nx, i
    real(pfdp) :: dx, x
    
    q0= 0.0_pfdp
        
    nx = size(q0,1)
    dx = Lx/dble(nx)
    
    do i=1,nx
       x = (i-1)*dx
       q0(i) = cos(x/16.0_pfdp)*(1+sin(x/16.0_pfdp))
    end do
    
  end subroutine initial_sol

  
  
  !> To set up optimization problem data (like desired state)
  subroutine initialize_ocp(sweeper, pf, level_index, nsteps_per_rank)
    class(pf_sweeper_t), intent(inout) :: sweeper
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index, nsteps_per_rank
    
    class(my_sweeper_t), pointer :: this
    this => as_my_sweeper(sweeper)
   
    this%nproc = pf%comm%nproc
    this%nsteps_per_rank = nsteps_per_rank
    this%myrank = pf%rank  ! does this work, or do we have to set it later?

    allocate(this%ydesired(nsteps_per_rank,pf%levels(level_index)%nnodes,pf%levels(level_index)%lev_shape(1)))
    this%ydesired = 0.0_pfdp
    
  end subroutine initialize_ocp
  
  
  subroutine objective_function(s, sol, shape, m, objective, step)
  ! actually just the tracking part of the objective
    class(pf_sweeper_t), intent(inout) :: s
    class(pf_encap_t), intent(in   )   :: sol
    integer,             intent(in)    :: shape(1), m, step
    real(pfdp),          intent(out)   :: objective

    real(pfdp),  pointer   :: y(:), f(:) !, obj(:)
    integer                :: nx,ny,nz,i,j,k !, nnodes
    

    class(my_sweeper_t), pointer :: sweeper
    sweeper => as_my_sweeper(s)

    nx = shape(1)
    allocate(f(nx))

       y => get_array1d_oc(sol, 1)
!        f = (y -sweeper%ydesiredT(:,:,:))
       f = (y -sweeper%ydesired(step,m,:))
       objective = 0.0_pfdp
       objective = sum(f**2)
       objective = objective * Lx / dble(nx)

    deallocate(f)
  end subroutine objective_function
  
  !> Set ydesired in sweeper to some space-time target state
  subroutine set_ydesired(sweeper, targetState)
     class(pf_sweeper_t), intent(inout) :: sweeper
     real(pfdp),          intent(in   ) :: targetState(:,:,:)
     
     class(my_sweeper_t), pointer :: this
     this => as_my_sweeper(sweeper)
     
     this%ydesired = targetState
  end subroutine set_ydesired
  
  
  
  !> Restrict desired state from fine to coarse level
  subroutine restrict_ydesired(sweeperC, sweeperF)
    class(pf_sweeper_t), intent(inout) :: sweeperC, sweeperF

    real(pfdp), pointer  :: ydesiredF(:,:,:), ydesiredC(:,:,:)
    integer :: nvarF, nvarC, xrat, nnodesF, nnodesC, trat, m
    
    class(my_sweeper_t), pointer :: sC, sF
    sC => as_my_sweeper(sweeperC)
    sF => as_my_sweeper(sweeperF)
    
    ydesiredF => sF%ydesired
    ydesiredC => sC%ydesired

    nnodesF = size(ydesiredF,2)
    nnodesC = size(ydesiredC,2)
    nvarF = size(ydesiredF,3)
    nvarC = size(ydesiredC,3)

    xrat  = nvarF / nvarC
    trat  = ceiling(real(nnodesF) / real(nnodesC))

    ydesiredC(:,:,:) = ydesiredF(:,::trat,::xrat)
  end subroutine restrict_ydesired
  
  subroutine restrict_for_adjoint(pf, t0, dt, which, step)
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: t0, dt
    integer, intent(in) :: which, step
    real(pfdp), pointer :: tF(:), tG(:) !zF(:), zG(:)
    integer :: l, m !, nnodesF, nnodesG, nvarsF, nvarsG

    do l = pf%state%finest_level, 2, -1
      call restrict_ts(pf%levels(l), pf%levels(l-1), pf%levels(l)%Q, pf%levels(l-1)%Q, t0+dt*pf%levels(l)%nodes, which)
    end do
  end subroutine restrict_for_adjoint
  
  !> compute an L2(Omega) scalar product between two functions (no time dependence here!)
  function compute_scalar_prod(f, g) result(r)
    real(pfdp), intent(in)  :: f(:), g(:)
    real(pfdp)              :: r
    integer                 :: m, i,j,k,nx,ny,nz, nsteps, n

    r = dot_product(f,g)
    r = r*Lx/dble(nx)
  end function compute_scalar_prod

end module my_sweeper
