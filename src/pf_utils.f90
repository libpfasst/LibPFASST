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

module pf_mod_utils
  use pf_mod_dtype
  use pf_mod_timer
  implicit none
contains

  !
  ! Build time interpolation matrix.
  !
  subroutine pf_time_interpolation_matrix(nodesF, nnodesF, nodesG, nnodesG, tmat)
    integer,    intent(in)  :: nnodesF, nnodesG
    real(pfdp), intent(in)  :: nodesF(0:nnodesF-1), nodesG(0:nnodesG-1)
    real(pfdp), intent(out) :: tmat(0:nnodesF-1,0:nnodesG-1)

    integer    :: i, j, k
    real(pfdp) :: xi, num, den

    do i = 0, nnodesF-1
       xi = nodesF(i)

       do j = 0, nnodesG-1
          den = 1.0_pfdp
          num = 1.0_pfdp

          do k = 0, nnodesG-1
             if (k == j) cycle
             den = den * (nodesG(j) - nodesG(k))
             num = num * (xi        - nodesG(k))
          end do

          tmat(i, j) = num/den
       end do
    end do
  end subroutine pf_time_interpolation_matrix


  !
  ! Spread initial condition.
  !
  subroutine spreadq0(Lev, t0)
    type(pf_level_t), intent(inout) :: Lev
    real(pfdp),       intent(in)    :: t0

    integer :: m, p

    call Lev%encap%unpack(Lev%Q(1), Lev%q0)
    call Lev%sweeper%evaluate(Lev, t0, 1)

    do m = 2, Lev%nnodes
       call Lev%encap%copy(Lev%Q(m), Lev%Q(1))
       do p = 1, Lev%sweeper%npieces
          call Lev%encap%copy(Lev%F(m,p), Lev%F(1,p))
       end do
    end do
  end subroutine spreadq0


  !
  ! Save current Q and F.
  !
  subroutine save(Lev)
    type(pf_level_t), intent(inout) :: Lev

    integer :: m, p

    if (Lev%Finterp) then
       !  Changed by MM Dec. 21, 2013 to save solutions too
       if (associated(Lev%pF)) then
          do m = 1, Lev%nnodes
             do p = 1,size(Lev%F(1,:))
                call Lev%encap%copy(Lev%pF(m,p), Lev%F(m,p))
             end do
             call Lev%encap%copy(Lev%pQ(m), Lev%Q(m))
          end do
!          call Lev%encap%copy(Lev%pQ(1), Lev%Q(1))
       end if
    else
       if (associated(Lev%pQ)) then
          do m = 1, Lev%nnodes
             call Lev%encap%copy(Lev%pQ(m), Lev%Q(m))
          end do
       end if
    end if
  end subroutine save


  !
  ! Compute full residual
  !
  ! During the process of computing the residual we compute the '0 to
  ! node' integral and store it in I.  This is used later when doing
  ! restriction (see restrict_time_space_fas).
  !
  subroutine pf_residual(pf, Lev, dt)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: Lev
    real(pfdp),        intent(in)    :: dt

    real(pfdp) :: norms(Lev%nnodes-1)
    integer :: m

!    if (pf%nlevels == 1 .and. pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) return
!   I think we often want the residual for diagnostics.  Maybe need flag to turn this off
!   for efficiency?

    call start_timer(pf, TRESIDUAL)

    call Lev%sweeper%residual(Lev, dt)

    ! compute max residual norm
    do m = 1, Lev%nnodes-1
       norms(m) = Lev%encap%norm(Lev%R(m))
    end do
    Lev%residual = maxval(abs(norms))
!    Lev%residual = norms(Lev%nnodes-1)

    call end_timer(pf, TRESIDUAL)

  end subroutine pf_residual


  !
  ! Apply an interpolation matrix (tmat or rmat) to src.
  !
  subroutine pf_apply_mat(dst, a, mat, src, encap, zero)
    type(c_ptr),       intent(inout) :: dst(:)
    real(pfdp),        intent(in)    :: a, mat(:, :)
    type(c_ptr),       intent(in)    :: src(:)
    type(pf_encap_t),  intent(in)    :: encap
    logical,           intent(in), optional :: zero

    logical :: lzero
    integer :: n, m, i, j

    lzero = .true.; if (present(zero)) lzero = zero

    n = size(mat, dim=1)
    m = size(mat, dim=2)

    ! XXX: test for nan's in matrices...

    do i = 1, n
       if (lzero) call encap%setval(dst(i), 0.0d0)
       do j = 1, m
          call encap%axpy(dst(i), a * mat(i, j), src(j))
       end do
    end do
  end subroutine pf_apply_mat

  !
  ! Generic residual
  !
  subroutine pf_generic_residual(Lev, dt)
    type(pf_level_t),  intent(inout) :: Lev
    real(pfdp),        intent(in)    :: dt

    integer :: m, n

    call Lev%sweeper%integrate(Lev, Lev%Q, Lev%F, dt, Lev%I)
!MMQ    do m = 2, Lev%nnodes-1
!       call Lev%encap%axpy(Lev%I(m), 1.0_pfdp, Lev%I(m-1))
!    end do

    ! add tau (which is 'node to node')
    if (associated(Lev%tau)) then
       do m = 1, Lev%nnodes-1
!  MMQ        do n = 1, m
!             call Lev%encap%axpy(Lev%I(m), 1.0_pfdp, Lev%tau(n))
!          end do
          call Lev%encap%axpy(Lev%I(m), 1.0_pfdp, Lev%tauQ(m))
       end do
    end if

    ! subtract out Q
    do m = 1, Lev%nnodes-1
       call Lev%encap%copy(Lev%R(m), Lev%I(m))
       call Lev%encap%axpy(Lev%R(m), 1.0_pfdp, Lev%Q(1))
       call Lev%encap%axpy(Lev%R(m), -1.0_pfdp, Lev%Q(m+1))
    end do

  end subroutine pf_generic_residual

  !
  ! Generic evaluate all
  !
  subroutine pf_generic_evaluate_all(Lev, t)
    type(pf_level_t),  intent(inout) :: Lev
    real(pfdp),        intent(in)    :: t(:)

    integer :: m
    do m = 1, Lev%nnodes
       call Lev%sweeper%evaluate(Lev, t(m), m)
    end do
  end subroutine pf_generic_evaluate_all


  subroutine pf_myLUexp(A,L,U,Nnodes)
    real(pfdp),       intent(in)    :: A(Nnodes,Nnodes)
    real(pfdp),      intent(inout)  :: L(Nnodes,Nnodes)
    real(pfdp),     intent(inout)   :: U(Nnodes,Nnodes)
    integer,        intent (in)     :: Nnodes
    ! Return the LU decomposition of an explicit integration matrix
    !   without pivoting
    integer :: i,j
    real(pfdp) :: c
    L = 0.0_pfdp
    U = 0.0_pfdp

    do i = 1,Nnodes-1
       L(i,i) = 1.0_pfdp
    end do
    U=transpose(A)
    do i = 1,Nnodes-1
       if (U(i,i+1) /= 0.0) then
          do j=i+1,Nnodes
             c = U(j,i+1)/U(i,i+1)
             U(j,i:Nnodes)=U(j,i:Nnodes)-c*U(i,i:Nnodes)
             L(j,:)=L(j,:)-c*L(i,:)
          end do
       end if
    end do

    !  Now scale the columns of U to match the sum of A
    do j=1,Nnodes
       c = sum(U(:,j))
       if (c /=  0.0) then
          U(:,j)=U(:,j)*sum(A(j,:))/c
       end if
    end do
  end subroutine pf_myLUexp

end module pf_mod_utils
