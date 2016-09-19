!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

module pf_mod_implicitQ
  use pf_mod_dtype
  use pf_mod_utils
  implicit none
  integer,parameter,private :: npieces = 1
  logical,save,private :: use_LUq_ = .false.

  interface
     subroutine pf_f2eval_p(y, t, level, ctx, f2)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, f2, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2eval_p
  end interface

  interface
     subroutine pf_f2comp_p(y, t, dt, rhs, level, ctx, f2)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, rhs, f2, ctx
       real(pfdp),     intent(in)        :: t, dt
       integer(c_int), intent(in)        :: level
     end subroutine pf_f2comp_p
  end interface

  type :: pf_implicitQ_t
     procedure(pf_f2eval_p),  pointer, nopass :: f2eval
     procedure(pf_f2comp_p),  pointer, nopass :: f2comp

     real(pfdp), ALLOCATABLE :: QdiffI(:,:)
     real(pfdp), ALLOCATABLE :: QtilI(:,:)
  end type pf_implicitQ_t

contains

  ! Perform one SDC sweep on level Lev and set qend appropriately.
  subroutine implicitQ_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:Lev%nnodes-1)
    type(c_ptr) :: rhs

    type(pf_implicitQ_t), pointer :: imp

    call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0_pfdp)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*imp%QdiffI(m,n), Lev%F(n,1))
       end do
       if (associated(Lev%tauQ)) then
          call Lev%encap%axpy(Lev%S(m), 1.0_pfdp, Lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%copy(Lev%Q(1), Lev%q0)

    call imp%f2eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,1))

    call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%ctx)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%setval(rhs, 0.0_pfdp)
       do n = 1, m
          call Lev%encap%axpy(rhs, dt*imp%QtilI(m,n), Lev%F(n,1))  
       end do
       call Lev%encap%axpy(rhs, 1.0_pfdp,lev%S(m))
       call Lev%encap%axpy(rhs, 1.0_pfdp, Lev%Q(1))

       call imp%f2comp(Lev%Q(m+1), t, dt*imp%QtilI(m,m+1), rhs, Lev%level, Lev%ctx, Lev%F(m+1,1))    

    end do
       
    ! Put the last node value into qend
    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine implicitQ_sweep

  ! Evaluate function values
  subroutine implicitQ_evaluate(Lev, t, m)
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_implicitQ_t), pointer :: imp
    call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    call imp%f2eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,1))
  end subroutine implicitQ_evaluate

  ! Initialize matrix
  subroutine implicitQ_initialize(Lev)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)

    integer :: m,n,nnodes
    type(pf_implicitQ_t), pointer :: imp
    call c_f_pointer(Lev%sweeper%sweeperctx, imp)

    nnodes = Lev%nnodes
    allocate(imp%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(imp%QtilI(nnodes-1,nnodes))  !  S-BE

    imp%QtilI = 0.0_pfdp
    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       do n = 1,m
          imp%QtilI(m,n+1) =  dsdc(n)
       end do
    end do

    if(use_LUq_) then
      call myLUq(Lev%qmat,imp%QtilI,Nnodes,1)
    end if

    imp%QdiffI = Lev%qmat-imp%QtilI

!    print *,'QtilI',imp%QtilI
!    print *,'QdiffI',imp%QdiffI
!    print *,'Qmat',Lev%qmat

  end subroutine implicitQ_initialize


  ! Compute SDC integral
  subroutine implicitQ_integrate(Lev, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in)    :: Lev
    type(c_ptr),      intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),       intent(in)    :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, Lev%nnodes-1
       call Lev%encap%setval(fintSDC(n), 0.0_pfdp)
       do m = 1, Lev%nnodes
          do p = 1, npieces
             call Lev%encap%axpy(fintSDC(n), dt*Lev%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine implicitQ_integrate

  ! Create implicitQ sweeper
  subroutine pf_implicitQ_create(sweeper,f2eval,f2comp,use_LUq)
    type(pf_sweeper_t),intent(inout) :: sweeper
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp
    logical,intent(in),optional :: use_LUq

    type(pf_implicitQ_t), pointer :: imp

    if(present(use_LUq)) then
      use_LUq_ = use_LUq
    endif

    allocate(imp)
    imp%f2eval => f2eval
    imp%f2comp => f2comp

    sweeper%npieces = npieces
    sweeper%sweep      => implicitQ_sweep
    sweeper%evaluate   => implicitQ_evaluate
    sweeper%destroy   => pf_implicitQ_destroy
    sweeper%initialize => implicitQ_initialize
    sweeper%integrate  => implicitQ_integrate
    sweeper%evaluate_all => pf_generic_evaluate_all
    sweeper%residual     => pf_generic_residual

    sweeper%sweeperctx = c_loc(imp)
  end subroutine pf_implicitQ_create

  subroutine pf_implicitQ_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_implicitQ_t), pointer :: imp
    call c_f_pointer(sweeper%sweeperctx, imp)
    deallocate(imp%QdiffI)
    deallocate(imp%QtilI)
    deallocate(imp)
  end subroutine pf_implicitQ_destroy

end module pf_mod_implicitQ

