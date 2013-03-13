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

  ! Build the time interpolation matrix
  subroutine pf_time_interpolation_matrix(nodesF, nnodesF, nodesG, nnodesG, tmat)
    implicit none
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

  ! Spread initial condition
  subroutine spreadq0(F, t0)
    type(pf_level_t), intent(inout) :: F
    real(pfdp),       intent(in)    :: t0

    integer :: m, p

    call F%encap%unpack(F%qSDC(1), F%q0)
    call F%sweeper%evaluate(F, t0, 1)

    do m = 2, F%nnodes
       call F%encap%copy(F%qSDC(m), F%qSDC(1))
       do p = 1, F%sweeper%npieces
          call F%encap%copy(F%fSDC(m,p), F%fSDC(1,p))
       end do
    end do
  end subroutine spreadq0

  ! Save current qSDC and fSDC
  subroutine save(F)
    type(pf_level_t), intent(inout) :: F

    integer :: m,p

    if (F%Finterp) then
       if (associated(F%pfSDC)) then
          do m = 1, F%nnodes
             do p = 1,size(F%fSDC(1,:))
                call F%encap%copy(F%pfSDC(m,p), F%fSDC(m,p))
             end do
          end do
          call F%encap%copy(F%pSDC(1), F%qSDC(1))
       end if
    else
       if (associated(F%pSDC)) then
          do m = 1, F%nnodes
             call F%encap%copy(F%pSDC(m), F%qSDC(m))
          end do
       end if
    end if
  end subroutine save

  ! Compute residual (generic but probably inefficient)
  subroutine pf_residual(F, dt, residual)
    type(pf_level_t), intent(inout) :: F
    real(pfdp),       intent(in)    :: dt
    type(c_ptr),      intent(in)    :: residual

    type(c_ptr) :: fintSDC(F%nnodes-1)
    integer :: n

    do n = 1, F%nnodes-1
       call F%encap%create(fintSDC(n), F%level, .true., F%nvars, F%shape, F%ctx)
    end do

    ! integrate and compute residual
    call F%sweeper%integrate(F, F%fSDC, dt, fintSDC)

    call F%encap%copy(residual, F%qSDC(1))
    do n = 1, F%nnodes-1
       call F%encap%axpy(residual, 1.0_pfdp, fintSDC(n))
    end do
    call F%encap%axpy(residual, -1.0_pfdp, F%qSDC(F%nnodes))

    do n = 1, F%nnodes-1
       call F%encap%destroy(fintSDC(n))
    end do
  end subroutine pf_residual

end module pf_mod_utils
