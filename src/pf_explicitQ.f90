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

module pf_mod_explicitQ
  use pf_mod_dtype
  use pf_mod_utils
  implicit none
  integer, parameter, private :: npieces = 1

  interface
     subroutine pf_f1eval_p(y, t, level, levelctx, f1)
       import c_ptr, c_int, pfdp, pf_context_t
       type(c_ptr),    intent(in), value :: y, f1!, levelctx
       class(pf_context_t), intent(in) :: levelctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_f1eval_p
  end interface

  type :: pf_explicitQ_t
     procedure(pf_f1eval_p),   pointer, nopass :: f1eval
     real(pfdp), ALLOCATABLE :: QdiffE(:,:)
     real(pfdp), ALLOCATABLE :: QtilE(:,:)
  end type pf_explicitQ_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine explicitQ_sweep(pf, Lev, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: Lev

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:Lev%nnodes-1)
    type(pf_explicitQ_t), pointer :: exp

    call c_f_pointer(Lev%sweeper%sweeperctx, exp)

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%S(m)%setval(0.0_pfdp)
       do n = 1, Lev%nnodes
          call Lev%S(m)%axpy(dt*exp%QdiffE(m,n), Lev%F(n,1))
       end do
       if (associated(Lev%tau)) then
          call Lev%S(m)%axpy(1.0_pfdp, Lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call Lev%Q(1)%unpack(Lev%q0)

    call exp%f1eval(Lev%Q(1), t0, Lev%level, Lev%levelctx, Lev%F(1,1))

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call Lev%Q(m+1)%copy(Lev%Q(1))
       do n = 1, m
          call Lev%Q(m+1)%axpy(dt*exp%QtilE(m,n), Lev%F(n,1))
       end do

!       call Lev%encap%axpy(Lev%Q(m+1), dtsdc(m), Lev%F(m,1))
       call Lev%Q(m+1)%axpy(1.0_pfdp, Lev%S(m))

       call exp%f1eval(Lev%Q(m+1), t, Lev%level, Lev%levelctx, Lev%F(m+1,1))
    end do

    call Lev%qend%copy(Lev%Q(Lev%nnodes))

    ! done
    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine explicitQ_sweep

  ! Evaluate function values
  subroutine explicitQ_evaluate(Lev, t, m)
    use pf_mod_dtype

    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: Lev

    type(pf_explicitQ_t), pointer :: exp

    call c_f_pointer(Lev%sweeper%sweeperctx, exp)

    call exp%f1eval(Lev%Q(m), t, Lev%level, Lev%levelctx, Lev%F(m,1))
  end subroutine explicitQ_evaluate

  ! Initialize matrix
  subroutine explicitQ_initialize(Lev)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: Lev
    real(pfdp) :: dsdc(Lev%nnodes-1)

    integer :: m,n,nnodes
    type(pf_explicitQ_t), pointer :: exp
    call c_f_pointer(Lev%sweeper%sweeperctx, exp)

    nnodes = Lev%nnodes
    allocate(exp%QdiffE(nnodes-1,nnodes))  ! S-FE
    allocate(exp%QtilE(nnodes-1,nnodes))  ! S-FE

    exp%QtilE = 0.0_pfdp

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       do n = 1,m
!          exp%QtilE(m,n)   =  dsdc(n)
       end do
    end do
    !  Or do the LU trick
!    call pf_LUexp(Lev%qmat,exp%QtilE,nnodes)

    exp%QdiffE = Lev%qmat-exp%QtilE

  end subroutine explicitQ_initialize


  ! Compute SDC integral
  subroutine explicitQ_integrate(Lev, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in)    :: Lev
    type(c_ptr),      intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),       intent(in)    :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, Lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, Lev%nnodes
          do p = 1, npieces
             call fintSDC(n)%axpy(dt*Lev%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine explicitQ_integrate

  ! Create explicitQ sweeper
  subroutine pf_explicitQ_create(sweeper, f1eval)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_f1eval_p) :: f1eval

    type(pf_explicitQ_t), pointer :: exp

    allocate(exp)
    exp%f1eval => f1eval

    sweeper%npieces = npieces
    sweeper%sweep      => explicitQ_sweep
    sweeper%evaluate   => explicitQ_evaluate
    sweeper%initialize => explicitQ_initialize
    sweeper%integrate  => explicitQ_integrate
    sweeper%destroy    => pf_explicitQ_destroy
    sweeper%evaluate_all => pf_generic_evaluate_all
    sweeper%residual     => pf_generic_residual

    sweeper%sweeperctx = c_loc(exp)
  end subroutine pf_explicitQ_create

  subroutine pf_explicitQ_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_explicitQ_t), pointer :: explicitQ
    call c_f_pointer(sweeper%sweeperctx, explicitQ)
    deallocate(explicitQ%QtilE)
    deallocate(explicitQ%QdiffE)

    deallocate(explicitQ)
  end subroutine pf_explicitQ_destroy

  subroutine pf_LUexp(A,U,Nnodes)
    real(pfdp),       intent(in)    :: A(Nnodes,Nnodes)
    real(pfdp),     intent(inout)   :: U(Nnodes,Nnodes)
    integer,        intent (in)     :: Nnodes
    ! Return the transpose of U from LU decomposition of
    !  an explicit integration matrix       without pivoting
    integer :: i,j
    real(pfdp) :: c

    U=transpose(A)
    do i = 1,Nnodes-1
       if (U(i,i+1) /= 0.0) then
          do j=i+1,Nnodes
             c = U(j,i+1)/U(i,i+1)
             U(j,i:Nnodes)=U(j,i:Nnodes)-c*U(i,i:Nnodes)
          end do
       end if
    end do
    U=transpose(U)

    !  Now scale the rows of U to match the sum of A
    do j=1,Nnodes
       c = sum(U(j,:))
       if (c /=  0.0) then
          U(j,:)=U(j,:)*sum(A(j,:))/c
!          print *,c,sum(A(j,:))/c
       end if
    end do

  end subroutine pf_LUexp

end module pf_mod_explicitQ
