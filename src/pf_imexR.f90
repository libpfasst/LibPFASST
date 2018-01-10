!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!
!>  This is the derived sweeper class for doing IMEX sweeps for an equation of the form
!>            y' = f_1(y) + f_2(y)
!>  The f_1 piece is treated explicitly and f_2 implicitly
!>  Afer this sweeper is initialized (usually in main), the locgical flags can be changed if desired
!>  The user need to supply the feval and fcomp routines for a given example   
!>     explicit:  Make false if there is no explicit piece
!>     implicit:  Make false if there is no implicit piece
!>     use_LUqt:  Make false if backward Euler sweepers should be used instead of the 'LU trick'
module pf_mod_imexR
  use pf_mod_dtype
  use pf_mod_utils

  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_imexR_t
     real(pfdp), allocatable :: QtilE(:,:)   !>  Approximate explcit quadrature rule
     real(pfdp), allocatable :: QtilI(:,:)   !>  Approximate implcit quadrature rule
     real(pfdp), allocatable :: dtsdc(:)     !>  SDC step sizes
     real(pfdp), allocatable :: QdiffE(:,:)  !>  qmat-QtilE
     real(pfdp), allocatable :: QdiffI(:,:)  !>  qmat-QtilI

     logical                 :: explicit = .true. !>  Is there an explicit piece
     logical                 :: implicit = .true. !>  Is there an implicit piece
     logical                 :: use_LUq = .true.  !>  Use the LU trick

     class(pf_encap_t), allocatable :: rhs   !> holds rhs for implicit solve

   contains
     procedure(pf_f_eval_p), deferred :: f_eval   !>  RHS function evaluations
     procedure(pf_f_comp_p), deferred :: f_comp   !>  Implicit solver
     !>  Set the generic functions
     procedure :: sweep      => imexR_sweep
     procedure :: initialize => imexR_initialize
     procedure :: evaluate   => imexR_evaluate
     procedure :: integrate  => imexR_integrate
     procedure :: residual   => imexR_residual
     procedure :: evaluate_all => imexR_evaluate_all
     procedure :: destroy   => imexR_destroy
     procedure :: imexR_destroy
  end type pf_imexR_t

  interface
     !>  This is the interface for the routine to compute the RHS function values
     subroutine pf_f_eval_p(this,y, t, level_index, f, piece)
       !>  Evaluae f_piece(y), where piece is one or two 
       import pf_imexR_t, pf_encap_t, pfdp
       class(pf_imexR_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y        !>  Argument for evaluation
       real(pfdp),        intent(in   ) :: t        !>  Time at evaluation
       integer,    intent(in   ) :: level_index     !>  Level index
       class(pf_encap_t), intent(inout) :: f        !>  RHS function value
       integer,    intent(in   ) :: piece           !>  Which piece to evaluate
     end subroutine pf_f_eval_p
     subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       !>  Solve the equation y - dtq*f_2(y) =rhs
       import pf_imexR_t, pf_encap_t, pfdp
       class(pf_imexR_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y      !>  Solution of implicit solve 
       real(pfdp),        intent(in   ) :: t      !>  Time of solve
       real(pfdp),        intent(in   ) :: dtq    !>  dt*quadrature weight
       class(pf_encap_t), intent(in   ) :: rhs    !>  RHS for solve
       integer,    intent(in   ) :: level_index   !>  Level index
       class(pf_encap_t), intent(inout) :: f      !>  f_2 of solution y
       integer,    intent(in   ) :: piece         !>  Which piece to evaluate
     end subroutine pf_f_comp_p
  end interface

contains

  !> Perform nsweep SDC sweep on level Lev and set qend appropriately.
  subroutine imexR_sweep(this, pf, level_index, t0, dt,nsweeps)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_imexR_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf  !>  PFASST structure
    real(pfdp),        intent(in   ) :: t0         !>  Time at beginning of time step
        real(pfdp),        intent(in   ) :: dt     !>  time step size
    integer,             intent(in)    :: level_index  !>  which level this is
    integer,             intent(in)    :: nsweeps      !>  number of sweeps to do

    class(pf_level_t), pointer :: lev    !>  points to current level

    integer     :: m, n,k   !>  Loop variables
    real(pfdp)  :: t        !>  Time at nodes
    
    lev => pf%levels(level_index)   !>  Assign level pointer

    call start_timer(pf, TLEVEL+lev%index-1)

    do k = 1,nsweeps   !>  Loop over sweeps
       call call_hooks(pf, level_index, PF_PRE_SWEEP)    

       ! compute integrals and add fas correction
       do m = 1, lev%nnodes-1
          call lev%I(m)%setval(0.0_pfdp)
          if (this%explicit) then
             do n = 1, lev%nnodes
                call lev%I(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
             end do
          end if
          if (this%implicit) then
             do n = 1, lev%nnodes
                call lev%I(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2))
             end do
          end if
          if (allocated(lev%tauQ)) then
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))
          end if
       end do
       !  Recompute the first function value if this is first sweep
       if (k .eq. 1) then
          call lev%Q(1)%copy(lev%q0)
          if (this%explicit) &
               call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1),1)
          if (this%implicit) &
               call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2),2)
       end if
       
       t = t0
       ! do the sub-stepping in sweep
       do m = 1, lev%nnodes-1
          t = t + dt*this%dtsdc(m)

          !>  Accumulate rhs
          call this%rhs%setval(0.0_pfdp)
          do n = 1, m
             if (this%explicit) &
                  call this%rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1))
             if (this%implicit) &
                  call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2))
          end do
          !  Add the tau term
          call this%rhs%axpy(1.0_pfdp, lev%I(m))

          !  Add the starting value
          call this%rhs%axpy(1.0_pfdp, lev%Q(1))

          !  Solve for the implicit piece
          if (this%implicit) then
             call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index,lev%F(m+1,2),2)
          else
             call lev%Q(m+1)%copy(this%rhs)
          end if
          !  Compute explicit function on new value
          if (this%explicit) &
               call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1),1)

          
       end do
       call pf_residual(pf, lev, dt)
       call lev%qend%copy(lev%Q(lev%nnodes))

       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do  !>  End loop on sweeps

    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine imexR_sweep

  !> Initialize matrices and space for sweeper
  subroutine imexR_initialize(this, lev)
    class(pf_imexR_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    integer    :: m,n, nnodes

    this%npieces = 2

    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes))  
    allocate(this%QdiffI(nnodes-1,nnodes))  
    allocate(this%QtilE(nnodes-1,nnodes))  
    allocate(this%QtilI(nnodes-1,nnodes))  
    allocate(this%dtsdc(nnodes-1))  

    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    !>  Array of substep sizes
    this%dtsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)

    ! Implicit matrix
    if (this%use_LUq) then 
       ! Get the LU
       call myLUq(lev%qmat,lev%LUmat,lev%nnodes,0)
       this%QtilI = lev%LUmat
!       print *,'LU',lev%LUmat
!       print *,'BE',lev%qmatBE
    else 
       this%QtilI =  lev%qmatBE
    end if

    ! Explicit matrix
    this%QtilE =  lev%qmatFE

    this%QdiffE = lev%qmat-this%QtilE
    this%QdiffI = lev%qmat-this%QtilI

    !>  Make space for rhs
    call lev%ulevel%factory%create_single(this%rhs, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

  end subroutine imexR_initialize

  subroutine imexR_destroy(this, lev)
    !>  deallocate
    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
    deallocate(this%dtsdc)

    call lev%ulevel%factory%destroy_single(this%rhs, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)


  end subroutine imexR_destroy


  !> Compute  Picard integral of function values
  subroutine imexR_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_imexR_t), intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in   ) :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          if (this%explicit) &
               call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,1))
          if (this%implicit) &
               call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,2))
       end do
    end do
  end subroutine imexR_integrate

  subroutine imexR_residual(this, lev, dt)
    !> Compute  Residual
    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: dt
    call pf_generic_residual(this, lev, dt)
  end subroutine imexR_residual
  
  subroutine imexR_evaluate(this, lev, t, m)
    !> Evaluate function value at node m

    class(pf_imexR_t),  intent(inout) :: this
    real(pfdp),        intent(in   ) :: t
    integer,           intent(in   ) :: m
    class(pf_level_t), intent(inout) :: lev
    if (this%explicit) &
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1),1)
    if (this%implicit) &
         call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2),2)
  end subroutine imexR_evaluate

  subroutine imexR_evaluate_all(this, lev, t)
    !> Evaluate all function values
    class(pf_imexR_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t(:)
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine imexR_evaluate_all

end module pf_mod_imexR
