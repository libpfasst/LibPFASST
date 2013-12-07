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

module pf_mod_imex
  use pf_mod_dtype
  use pf_mod_explicit, only: pf_f1eval_p
  use pf_mod_implicit, only: pf_f2eval_p, pf_f2comp_p
  implicit none

  integer, parameter, private :: npieces = 2

  type :: pf_imex_t
     procedure(pf_f1eval_p), pointer, nopass :: f1eval
     procedure(pf_f2eval_p), pointer, nopass :: f2eval
     procedure(pf_f2comp_p), pointer, nopass :: f2comp

     !  Matrices
     real(pfdp), ALLOCATABLE :: SdiffE(:,:)
     real(pfdp), ALLOCATABLE :: SdiffI(:,:)

  end type pf_imex_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine imex_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:Lev%nnodes-1)
    type(c_ptr) :: rhs

    type(pf_imex_t), pointer :: imex

    call c_f_pointer(Lev%sweeper%sweeperctx, imex)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0d0)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*imex%SdiffE(m,n), Lev%F(n,1))
          call Lev%encap%axpy(Lev%S(m), dt*imex%SdiffI(m,n), Lev%F(n,2))
       end do
       if (associated(Lev%tau)) then
          call Lev%encap%axpy(Lev%S(m), 1.0d0, Lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%unpack(Lev%Q(1), Lev%q0)

    call imex%f1eval(Lev%Q(1), t0, Lev%level, Lev%levelctx, Lev%F(1,1))
    call imex%f2eval(Lev%Q(1), t0, Lev%level, Lev%levelctx, Lev%F(1,2))

    call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%levelctx, Lev%encap%encapctx)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%copy(rhs, Lev%Q(m))
       call Lev%encap%axpy(rhs, dtsdc(m), Lev%F(m,1))
       call Lev%encap%axpy(rhs, 1.0d0, Lev%S(m))

       call imex%f2comp(Lev%Q(m+1), t, dtsdc(m), rhs, Lev%level, Lev%levelctx, Lev%F(m+1,2))
       call imex%f1eval(Lev%Q(m+1), t, Lev%level, Lev%levelctx, Lev%F(m+1,1))
    end do

    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine imex_sweep

  ! Evaluate function values
  subroutine imex_evaluate(Lev, t, m)
    real(pfdp),       intent(in   ) :: t
    integer,          intent(in   ) :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_imex_t), pointer :: imex
    call c_f_pointer(Lev%sweeper%sweeperctx, imex)

    call imex%f1eval(Lev%Q(m), t, Lev%level, Lev%levelctx, Lev%F(m,1))
    call imex%f2eval(Lev%Q(m), t, Lev%level, Lev%levelctx, Lev%F(m,2))
  end subroutine imex_evaluate

  ! Initialize matrices
  subroutine imex_initialize(Lev)
    type(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)
    integer    :: m, nnodes

    type(pf_imex_t), pointer :: imex
    call c_f_pointer(Lev%sweeper%sweeperctx, imex)

    nnodes = Lev%nnodes
    allocate(imex%SdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(imex%SdiffI(nnodes-1,nnodes))  !  S-BE

    imex%SdiffE = Lev%s0mat
    imex%SdiffI = Lev%s0mat

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       imex%SdiffE(m,m)   = imex%SdiffE(m,m)   - dsdc(m)
       imex%SdiffI(m,m+1) = imex%SdiffI(m,m+1) - dsdc(m)
    end do
  end subroutine imex_initialize

  ! Compute SDC integral
  subroutine imex_integrate(Lev, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in)    :: Lev
    type(c_ptr),      intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),       intent(in)    :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, Lev%nnodes-1
       call Lev%encap%setval(fintSDC(n), 0.0d0)
       do m = 1, Lev%nnodes
          do p = 1, npieces
             call Lev%encap%axpy(fintSDC(n), dt*Lev%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine imex_integrate

  ! Create/destroy IMEX sweeper
  subroutine pf_imex_create(sweeper, f1eval, f2eval, f2comp)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp

    type(pf_imex_t), pointer :: imex

    allocate(imex)
    imex%f1eval => f1eval
    imex%f2eval => f2eval
    imex%f2comp => f2comp

    sweeper%npieces = npieces
    sweeper%sweep      => imex_sweep
    sweeper%evaluate   => imex_evaluate
    sweeper%initialize => imex_initialize
    sweeper%integrate  => imex_integrate
    sweeper%destroy    => pf_imex_destroy

    sweeper%sweeperctx = c_loc(imex)
  end subroutine pf_imex_create

  subroutine pf_imex_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_imex_t), pointer :: imex
    call c_f_pointer(sweeper%sweeperctx, imex)

    deallocate(imex%SdiffI)
    deallocate(imex%SdiffE)
    deallocate(imex)
  end subroutine pf_imex_destroy

end module pf_mod_imex

