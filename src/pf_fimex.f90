!
! Copyright (C) 2014 Matthew Emmett.
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

module pf_mod_fimex
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_explicit, only: pf_f1eval_p
  use pf_mod_implicit, only: pf_f2eval_p, pf_f2comp_p
  implicit none

  integer, parameter, private :: npieces = 2

  interface
     subroutine pf_force_p(y, t, level, ctx, f)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, f, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_force_p
  end interface

  type :: pf_fimex_t
     procedure(pf_f1eval_p), pointer, nopass :: f1eval
     procedure(pf_f2eval_p), pointer, nopass :: f2eval
     procedure(pf_f2comp_p), pointer, nopass :: f2comp
     procedure(pf_force_p),  pointer, nopass :: force

     real(pfdp), allocatable :: SdiffE(:,:)
     real(pfdp), allocatable :: SdiffI(:,:)

     real(pfdp) :: last_forcing_time

     type(c_ptr), pointer :: F(:) ! forcing terms

  end type pf_fimex_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine fimex_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:Lev%nnodes-1)
    type(c_ptr) :: rhs

    type(pf_fimex_t), pointer :: fimex

    call c_f_pointer(Lev%sweeper%sweeperctx, fimex)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute forcing if necessary
    if (fimex%last_forcing_time /= t0) then
       t = t0
       dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
       do m = 1, Lev%nnodes
          if (.not. c_associated(fimex%F(m))) then
             call Lev%encap%create(fimex%F(m), Lev%level, SDC_KIND_SOL_FEVAL, &
                  Lev%nvars, Lev%shape, Lev%ctx)
          end if
          call fimex%force(Lev%Q(m), t, Lev%level, Lev%ctx, fimex%F(m))

          if (m < Lev%nnodes) t = t + dtsdc(m)
       end do
       fimex%last_forcing_time = t0
    end if

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0d0)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*Lev%s0mat(m,n), fimex%F(n))
          call Lev%encap%axpy(Lev%S(m), dt*fimex%SdiffE(m,n), Lev%F(n,1))
          call Lev%encap%axpy(Lev%S(m), dt*fimex%SdiffI(m,n), Lev%F(n,2))
       end do
       if (associated(Lev%tau)) then
          call Lev%encap%axpy(Lev%S(m), 1.0d0, Lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%copy(Lev%Q(1), Lev%q0)

    call fimex%f1eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,1))
    call fimex%f2eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,2))

    call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%ctx)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%copy(rhs, Lev%Q(m))
       call Lev%encap%axpy(rhs, dtsdc(m), Lev%F(m,1))
       call Lev%encap%axpy(rhs, 1.0d0, Lev%S(m))

       call fimex%f2comp(Lev%Q(m+1), t, dtsdc(m), rhs, Lev%level, Lev%ctx, Lev%F(m+1,2))
       call fimex%f1eval(Lev%Q(m+1), t, Lev%level, Lev%ctx, Lev%F(m+1,1))
    end do

    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine fimex_sweep

  ! Evaluate function values
  subroutine fimex_evaluate(Lev, t, m)
    real(pfdp),       intent(in   ) :: t
    integer,          intent(in   ) :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_fimex_t), pointer :: fimex
    call c_f_pointer(Lev%sweeper%sweeperctx, fimex)

    call fimex%f1eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,1))
    call fimex%f2eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,2))
  end subroutine fimex_evaluate

  ! Initialize matrices
  subroutine fimex_initialize(Lev)
    type(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)
    integer    :: m, nnodes

    type(pf_fimex_t), pointer :: fimex
    call c_f_pointer(Lev%sweeper%sweeperctx, fimex)

    nnodes = Lev%nnodes
    allocate(fimex%SdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(fimex%SdiffI(nnodes-1,nnodes))  !  S-BE

    fimex%SdiffE = Lev%s0mat
    fimex%SdiffI = Lev%s0mat

    fimex%last_forcing_time = -1
    allocate(fimex%F(nnodes))
    do m = 1, nnodes
       fimex%F(m) = c_null_ptr
    end do

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       fimex%SdiffE(m,m)   = fimex%SdiffE(m,m)   - dsdc(m)
       fimex%SdiffI(m,m+1) = fimex%SdiffI(m,m+1) - dsdc(m)
    end do
  end subroutine fimex_initialize

  ! Compute SDC integral
  subroutine fimex_integrate(Lev, qSDC, fSDC, dt, fintSDC)
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
  end subroutine fimex_integrate

  ! Create/destroy FIMEX sweeper
  subroutine pf_fimex_create(sweeper, f1eval, f2eval, f2comp, force)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp
    procedure(pf_force_p)  :: force

    type(pf_fimex_t), pointer :: fimex

    allocate(fimex)
    fimex%f1eval => f1eval
    fimex%f2eval => f2eval
    fimex%f2comp => f2comp
    fimex%force  => force

    sweeper%npieces = npieces
    sweeper%sweep        => fimex_sweep
    sweeper%evaluate     => fimex_evaluate
    sweeper%initialize   => fimex_initialize
    sweeper%integrate    => fimex_integrate
    sweeper%destroy      => pf_fimex_destroy
    sweeper%evaluate_all => pf_generic_evaluate_all
    sweeper%residual     => pf_generic_residual

    sweeper%sweeperctx = c_loc(fimex)
  end subroutine pf_fimex_create

  subroutine pf_fimex_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_fimex_t), pointer :: fimex
    call c_f_pointer(sweeper%sweeperctx, fimex)

    deallocate(fimex%SdiffI)
    deallocate(fimex%SdiffE)
    deallocate(fimex)
  end subroutine pf_fimex_destroy

end module pf_mod_fimex

