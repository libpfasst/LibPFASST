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
  subroutine spreadq0(F, t0)
    type(pf_level_t), intent(inout) :: F
    real(pfdp),       intent(in)    :: t0

    integer :: m, p

    call F%encap%unpack(F%Q(1), F%q0)
    call F%sweeper%evaluate(F, t0, 1)

    do m = 2, F%nnodes
       call F%encap%copy(F%Q(m), F%Q(1))
       do p = 1, F%sweeper%npieces
          call F%encap%copy(F%F(m,p), F%F(1,p))
       end do
    end do
  end subroutine spreadq0


  !
  ! Save current Q and F.
  !
  subroutine save(F)
    type(pf_level_t), intent(inout) :: F

    integer :: m, p

    if (F%Finterp) then
       if (associated(F%pF)) then
          do m = 1, F%nnodes
             do p = 1,size(F%F(1,:))
                call F%encap%copy(F%pF(m,p), F%F(m,p))
             end do
          end do
          call F%encap%copy(F%pQ(1), F%Q(1))
       end if
    else
       if (associated(F%pQ)) then
          do m = 1, F%nnodes
             call F%encap%copy(F%pQ(m), F%Q(m))
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
  subroutine pf_residual(F, dt)
    type(pf_level_t), intent(inout) :: F
    real(pfdp),       intent(in)    :: dt

    real(pfdp) :: norms(F%nnodes-1)
    integer :: m, n

    call F%sweeper%integrate(F, F%Q, F%F, dt, F%I)
    do m = 2, F%nnodes-1
       call F%encap%axpy(F%I(m), 1.0_pfdp, F%I(m-1))
    end do

    ! add tau (which is 'node to node')
    if (associated(F%tau)) then
       do m = 1, F%nnodes-1
          do n = 1, m
             call F%encap%axpy(F%I(m), 1.0_pfdp, F%tau(n))
          end do
       end do
    end if


    ! subtract out Q
    do m = 1, F%nnodes-1
       call F%encap%copy(F%R(m), F%Q(1))
       call F%encap%axpy(F%R(m),  1.0_pfdp, F%I(m))
       call F%encap%axpy(F%R(m), -1.0_pfdp, F%Q(m+1))
    end do

    
    ! compute max residual norm
    do m = 1, F%nnodes-1
       norms(m) = F%encap%norm(F%R(m))
    end do

    F%residual = maxval(abs(norms))

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

end module pf_mod_utils
