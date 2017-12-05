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

module pf_mod_imexQ
  use pf_mod_dtype
  use pf_mod_utils

  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_imexQ_t
     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
     logical                 :: use_LUq = .false.
     logical                 :: explicit = .true.
     logical                 :: implicit = .true.
   contains
     procedure(pf_f_eval_p), deferred :: f_eval
     procedure(pf_f_comp_p), deferred :: f_comp
     procedure :: sweep      => imexQ_sweep
     procedure :: initialize => imexQ_initialize
     procedure :: evaluate   => imexQ_evaluate
     procedure :: integrate  => imexQ_integrate
     procedure :: residual   => imexQ_residual
     procedure :: evaluate_all => imexQ_evaluate_all
     procedure :: destroy   => imexQ_destroy
     procedure :: imexQ_destroy
  end type pf_imexQ_t

    interface
     subroutine pf_f_eval_p(this,y, t, level, f, piece)
       import pf_imexQ_t, pf_encap_t, pfdp
       class(pf_imexQ_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y
       real(pfdp),        intent(in   ) :: t
       integer,    intent(in   ) :: level
       class(pf_encap_t), intent(inout) :: f
       integer,    intent(in   ) :: piece
     end subroutine pf_f_eval_p
      subroutine pf_f_comp_p(this,y, t, dt, rhs, level, f, piece)
       import pf_imexQ_t, pf_encap_t, pfdp
       class(pf_imexQ_t),  intent(inout) :: this
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

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine imexQ_sweep(this, pf, level_index, t0, dt,nsweeps)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_imexQ_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    integer,             intent(in)    :: level_index
    integer,             intent(in)    :: nsweeps

    class(pf_level_t), pointer :: lev

    integer     :: m, n,k
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:pf%levels(level_index)%nnodes-1)
    class(pf_encap_t), allocatable :: rhs

    lev => pf%levels(level_index)
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    call lev%ulevel%factory%create_single(rhs, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call start_timer(pf, TLEVEL+lev%index-1)
    print *,'sweeping',this%explicit,this%implicit
    do k = 1,nsweeps
       call call_hooks(pf, level_index, PF_PRE_SWEEP)    
       ! compute integrals and add fas correction
       do m = 1, lev%nnodes-1
          call lev%S(m)%setval(0.0_pfdp)
          if (this%explicit) then
             do n = 1, lev%nnodes
                call lev%S(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
             end do
          end if
          if (this%implicit) then
             do n = 1, lev%nnodes
                call lev%S(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2))
             end do
          end if
          if (allocated(lev%tauQ)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tauQ(m))
          end if
       end do
       !  Recompute the first function value
       if (k .eq. 1) then
          call lev%Q(1)%copy(lev%q0)
          if (this%explicit) &
               call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1),1)
          if (this%implicit) &
               call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2),2)
       end if
       
       t = t0
       ! do the time-stepping
       do m = 1, lev%nnodes-1
          t = t + dtsdc(m)
          
          call rhs%setval(0.0_pfdp)
          do n = 1, m
             if (this%explicit) &
                  call rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1))
             if (this%implicit) &
                  call rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2))
          end do
          !  Add the tau term
          call rhs%axpy(1.0_pfdp, lev%S(m))
          !  Add the starting value
          call rhs%axpy(1.0_pfdp, lev%Q(1))
          if (this%implicit) then
             call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), rhs, lev%index,lev%F(m+1,2),2)
          else
             call lev%Q(m+1)%copy(rhs)
          end if
          if (this%explicit) &
               call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1),1)

          
       end do
       call pf_residual(pf, lev, dt)
       call lev%qend%copy(lev%Q(lev%nnodes))

       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do

    call lev%ulevel%factory%destroy_single(rhs, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine imexQ_sweep

  ! Initialize matrices
  subroutine imexQ_initialize(this, lev)
    class(pf_imexQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m,n, nnodes

    this%npieces = 2

    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(this%QtilE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QtilI(nnodes-1,nnodes))  !  S-BE

    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)

    ! Implicit matrix
    if (this%use_LUq) then 
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

  end subroutine imexQ_initialize

  subroutine imexQ_destroy(this, lev)
    class(pf_imexQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
  end subroutine imexQ_destroy



  ! Compute SDC integral
  subroutine imexQ_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_imexQ_t), intent(inout) :: this
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
  end subroutine imexQ_integrate

    subroutine imexQ_residual(this, lev, dt)
    class(pf_imexQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: dt
    call pf_generic_residual(this, lev, dt)
  end subroutine imexQ_residual
  
  subroutine imexQ_evaluate(this, lev, t, m)
    class(pf_imexQ_t),  intent(inout) :: this
    real(pfdp),        intent(in   ) :: t
    integer,           intent(in   ) :: m
    class(pf_level_t), intent(inout) :: lev
    if (this%explicit) &
       call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1),1)
    if (this%implicit) &
         call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2),2)
  end subroutine imexQ_evaluate

  subroutine imexQ_evaluate_all(this, lev, t)
    class(pf_imexQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t(:)
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine imexQ_evaluate_all

end module pf_mod_imexQ
