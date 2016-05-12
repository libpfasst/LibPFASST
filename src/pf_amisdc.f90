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
!  Asynchronous misdc
module pf_mod_amisdc
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_explicitQ, only: pf_f1eval_p
  use pf_mod_implicitQ, only: pf_f2eval_p, pf_f2comp_p
  implicit none

  integer, parameter, private :: npieces = 3

  interface
     subroutine pf_f3eval_p(y, t, level, ctx, f3)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, f3, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f3eval_p
  end interface

  interface
     subroutine pf_f3comp_p(y, t, dt, rhs, level, ctx, f3)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: y, rhs, f3, ctx
       real(pfdp),     intent(in)        :: t, dt
       integer(c_int), intent(in)        :: level
     end subroutine pf_f3comp_p
  end interface
  type :: pf_amisdc_t
     procedure(pf_f1eval_p), pointer, nopass :: f1eval
     procedure(pf_f2eval_p), pointer, nopass :: f2eval
     procedure(pf_f2comp_p), pointer, nopass :: f2comp
     procedure(pf_f3eval_p), pointer, nopass :: f3eval
     procedure(pf_f3comp_p), pointer, nopass :: f3comp

     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
  end type pf_amisdc_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine amisdc_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:Lev%nnodes-1)
    type(c_ptr) :: rhs
    type(c_ptr), pointer :: S3(:)

    type(pf_amisdc_t), pointer :: amisdc

    call c_f_pointer(Lev%sweeper%sweeperctx, amisdc)

    call start_timer(pf, TLEVEL+Lev%level-1)
    allocate(S3(Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       call Lev%encap%create(S3(m),Lev%level,SDC_KIND_SOL_FEVAL,Lev%nvars,Lev%shape,Lev%ctx)
    end do

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%encap%setval(Lev%S(m), 0.0_pfdp)
       call Lev%encap%setval(S3(m), 0.0d0)
       do n = 1, Lev%nnodes
          call Lev%encap%axpy(Lev%S(m), dt*amisdc%QdiffE(m,n), Lev%F(n,1))
          call Lev%encap%axpy(Lev%S(m), dt*amisdc%QdiffI(m,n), Lev%F(n,2))
          call Lev%encap%axpy(Lev%S(m), dt*Lev%qmat(m,n), Lev%F(n,3))
          call Lev%encap%axpy(S3(m), dt*amisdc%QtilI(m,n), Lev%F(n,3))
          !  Note we have to leave off the -dt*Qtil here and put it in after f2comp
       end do
       if (associated(Lev%tauQ)) then
          call Lev%encap%axpy(Lev%S(m), 1.0_pfdp, Lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call Lev%encap%copy(Lev%Q(1), Lev%q0)

    call amisdc%f1eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,1))
    call amisdc%f2eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,2))
    call amisdc%f3eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,3))

    call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%ctx)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%encap%setval(rhs, 0.0_pfdp)
       do n = 1, m
          call Lev%encap%axpy(rhs, dt*amisdc%QtilE(m,n), Lev%F(n,1))  
          call Lev%encap%axpy(rhs, dt*amisdc%QtilI(m,n), Lev%F(n,2))  
       end do
       !  Add the tau term
       call Lev%encap%axpy(rhs, 1.0_pfdp, Lev%S(m))
       !  Add the starting value
       call Lev%encap%axpy(rhs,1.0_pfdp, Lev%Q(1))

       call amisdc%f2comp(Lev%Q(m+1), t, dt*amisdc%QtilI(m,m+1), rhs, Lev%level, Lev%ctx, Lev%F(m+1,2))

       !  Now we need to do the final subtraction for the f3 piece
       call Lev%encap%copy(rhs, Lev%Q(m+1))       
       do n = 1, m
          call Lev%encap%axpy(rhs, dt*amisdc%QtilI(m,n), Lev%F(n,3))  
       end do

       call Lev%encap%axpy(rhs, -1.0_pfdp, S3(m))

       call amisdc%f3comp(Lev%Q(m+1), t, dt*amisdc%QtilI(m,m+1), rhs, Lev%level, Lev%ctx, Lev%F(m+1,3))
       call amisdc%f1eval(Lev%Q(m+1), t, Lev%level, Lev%ctx, Lev%F(m+1,1))
       call amisdc%f2eval(Lev%Q(m+1), t, Lev%level, Lev%ctx, Lev%F(m+1,2))
    end do

    call Lev%encap%copy(Lev%qend, Lev%Q(Lev%nnodes))

    ! done
    call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine amisdc_sweep

  ! Evaluate function values
  subroutine amisdc_evaluate(Lev, t, m)
    real(pfdp),       intent(in   ) :: t
    integer,          intent(in   ) :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_amisdc_t), pointer :: amisdc
    call c_f_pointer(Lev%sweeper%sweeperctx, amisdc)

    call amisdc%f1eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,1))
    call amisdc%f2eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,2))
    call amisdc%f3eval(Lev%Q(m), t, Lev%level, Lev%ctx, Lev%F(m,3))
  end subroutine amisdc_evaluate

  ! Initialize matrices
  subroutine amisdc_initialize(Lev)
    type(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)
    integer    :: m,n, nnodes

    type(pf_amisdc_t), pointer :: amisdc
    call c_f_pointer(Lev%sweeper%sweeperctx, amisdc)

    nnodes = Lev%nnodes
    allocate(amisdc%QdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(amisdc%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(amisdc%QtilE(nnodes-1,nnodes))  !  S-FE
    allocate(amisdc%QtilI(nnodes-1,nnodes))  !  S-BE

    amisdc%QtilE = 0.0_pfdp
    amisdc%QtilI = 0.0_pfdp

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       do n = 1,m
          amisdc%QtilE(m,n)   =  dsdc(n)
          amisdc%QtilI(m,n+1) =  dsdc(n)
       end do
    end do

!    do m = 1,nnodes-1
!       print *,'row i of qmat', m,Lev%qmat(m,:)
!    end do
!    call myLUq(Lev%qmat,amisdc%QtilI,Nnodes,0)
    amisdc%QdiffE = Lev%qmat-amisdc%QtilE
    amisdc%QdiffI = Lev%qmat-amisdc%QtilI

  end subroutine amisdc_initialize

  ! Compute SDC integral
  subroutine amisdc_integrate(Lev, qSDC, fSDC, dt, fintSDC)
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
  end subroutine amisdc_integrate

  ! Create/destroy IMEXQ sweeper
  subroutine pf_amisdc_create(sweeper, f1eval, f2eval, f2comp,f3eval,f3comp)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval
    procedure(pf_f2eval_p) :: f2eval
    procedure(pf_f2comp_p) :: f2comp
    procedure(pf_f3eval_p) :: f3eval
    procedure(pf_f3comp_p) :: f3comp

    type(pf_amisdc_t), pointer :: amisdc

    allocate(amisdc)
    amisdc%f1eval => f1eval
    amisdc%f2eval => f2eval
    amisdc%f2comp => f2comp
    amisdc%f3eval => f3eval
    amisdc%f3comp => f3comp

    sweeper%npieces = npieces
    sweeper%sweep        => amisdc_sweep
    sweeper%evaluate     => amisdc_evaluate
    sweeper%initialize   => amisdc_initialize
    sweeper%integrate    => amisdc_integrate
    sweeper%destroy      => pf_amisdc_destroy
    sweeper%evaluate_all => pf_generic_evaluate_all
    sweeper%residual     => pf_generic_residual

    sweeper%sweeperctx = c_loc(amisdc)
  end subroutine pf_amisdc_create

  subroutine pf_amisdc_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_amisdc_t), pointer :: amisdc
    call c_f_pointer(sweeper%sweeperctx, amisdc)


    deallocate(amisdc%QdiffI)
    deallocate(amisdc%QdiffE)
    deallocate(amisdc%QtilI)
    deallocate(amisdc%QtilE)
    deallocate(amisdc)
  end subroutine pf_amisdc_destroy

end module pf_mod_amisdc

