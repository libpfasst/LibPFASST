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
  use pf_mod_explicitQ, only: pf_f1eval_p
  use pf_mod_implicitQ, only: pf_f2eval_p, pf_f2comp_p
  implicit none

  integer, parameter, private :: npieces = 2

  type :: pf_imexQ_t
     procedure(pf_f1eval_p), pointer, nopass :: f1eval
     procedure(pf_f2eval_p), pointer, nopass :: f2eval
     procedure(pf_f2comp_p), pointer, nopass :: f2comp

     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
  end type pf_imexQ_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine imexQ_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:Lev%nnodes-1)
    type(c_ptr) :: rhs

    type(pf_imexQ_t), pointer :: imexQ

    call c_f_pointer(Lev%sweeper%sweeperctx, imexQ)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0_pfdp)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*imexQ%QdiffE(m,n), Lev%F(n,1))
          call Lev%encap%axpy(Lev%S(m), dt*imexQ%QdiffI(m,n), Lev%F(n,2))
       end do
       if (associated(Lev%tauQ)) then
          call Lev%encap%axpy(Lev%S(m), 1.0_pfdp, Lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%unpack(Lev%Q(1), Lev%q0)

    call imexQ%f1eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,1))
    call imexQ%f2eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,2))

    call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%ctx)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%setval(rhs, 0.0_pfdp)
       do n = 1, m
          call Lev%encap%axpy(rhs, dt*imexQ%QtilE(m,n), Lev%F(n,1))  
          call Lev%encap%axpy(rhs, dt*imexQ%QtilI(m,n), Lev%F(n,2))  
       end do


!       call Lev%encap%axpy(rhs, dtsdc(m), Lev%F(m,1))
       call Lev%encap%axpy(rhs, 1.0_pfdp, Lev%S(m))
       !  Add the starting value
       call Lev%encap%axpy(rhs,1.0_pfdp, Lev%Q(1))


!       call imexQ%f2comp(Lev%Q(m+1), t, dtsdc(m), rhs, Lev%level, Lev%ctx, Lev%F(m+1,2))
       call imexQ%f2comp(Lev%Q(m+1), t, dt*imexQ%QtilI(m,m+1), rhs, Lev%level, Lev%ctx, Lev%F(m+1,2))
       call imexQ%f1eval(Lev%Q(m+1), t, Lev%level, Lev%ctx, Lev%F(m+1,1))

    end do

    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine imexQ_sweep

  ! Evaluate function values
  subroutine imexQ_evaluate(Lev, t, m)
    real(pfdp),       intent(in   ) :: t
    integer,          intent(in   ) :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_imexQ_t), pointer :: imexQ
    call c_f_pointer(Lev%sweeper%sweeperctx, imexQ)

    call imexQ%f1eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,1))
    call imexQ%f2eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,2))
  end subroutine imexQ_evaluate

  ! Initialize matrices
  subroutine imexQ_initialize(Lev)
    type(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)
    integer    :: m,n, nnodes

    type(pf_imexQ_t), pointer :: imexQ
    call c_f_pointer(Lev%sweeper%sweeperctx, imexQ)

    nnodes = Lev%nnodes
    allocate(imexQ%QdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(imexQ%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(imexQ%QtilE(nnodes-1,nnodes))  !  S-FE
    allocate(imexQ%QtilI(nnodes-1,nnodes))  !  S-BE

    imexQ%QtilE = 0.0_pfdp
    imexQ%QtilI = 0.0_pfdp

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       do n = 1,m
          imexQ%QtilE(m,n)   =  dsdc(n)
          imexQ%QtilI(m,n+1) =  dsdc(n)
       end do
    end do

!    do m = 1,nnodes-1
!       print *,'row i of qmat', m,Lev%qmat(m,:)
!    end do
!    call myLUq(Lev%qmat,imexQ%QtilI,Nnodes,0)
    imexQ%QdiffE = Lev%qmat-imexQ%QtilE
    imexQ%QdiffI = Lev%qmat-imexQ%QtilI

  end subroutine imexQ_initialize

  ! Compute SDC integral
  subroutine imexQ_integrate(Lev, qSDC, fSDC, dt, fintSDC)
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
  end subroutine imexQ_integrate

  ! Create/destroy IMEXQ sweeper
  subroutine pf_imexQ_create(sweeper, f1eval, f2eval, f2comp)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp

    type(pf_imexQ_t), pointer :: imexQ

    allocate(imexQ)
    imexQ%f1eval => f1eval
    imexQ%f2eval => f2eval
    imexQ%f2comp => f2comp

    sweeper%npieces = npieces
    sweeper%sweep        => imexQ_sweep
    sweeper%evaluate     => imexQ_evaluate
    sweeper%initialize   => imexQ_initialize
    sweeper%integrate    => imexQ_integrate
    sweeper%destroy      => pf_imexQ_destroy
    sweeper%evaluate_all => pf_generic_evaluate_all
    sweeper%residual     => pf_generic_residual

    sweeper%sweeperctx = c_loc(imexQ)
  end subroutine pf_imexQ_create

  subroutine pf_imexQ_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_imexQ_t), pointer :: imexQ
    call c_f_pointer(sweeper%sweeperctx, imexQ)


    deallocate(imexQ%QdiffI)
    deallocate(imexQ%QdiffE)
    deallocate(imexQ%QtilI)
    deallocate(imexQ%QtilE)
    deallocate(imexQ)
  end subroutine pf_imexQ_destroy

end module pf_mod_imexQ

