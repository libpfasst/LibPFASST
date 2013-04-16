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

module pf_mod_quadrature
  use pf_mod_dtype
  use sdc_mod_quadrature
  implicit none
contains

  subroutine pf_quadrature(qtype, nnodes, nnodes0, nodes, nflags, smat, qmat)
    integer,    intent(in)  :: qtype, nnodes, nnodes0
    real(pfdp), intent(out) :: nodes(nnodes)
    real(pfdp), intent(out) :: smat(nnodes-1,nnodes), qmat(nnodes-1,nnodes)
    integer,    intent(out) :: nflags(nnodes)

    real(c_long_double) :: qnodes0(nnodes0), qnodes(nnodes)
    integer(c_int)      :: flags0(nnodes0)

    integer :: refine

    if (qtype > SDC_PROPER_NODES) then

       ! nodes are given by proper quadrature rules

       call sdc_qnodes(qnodes, nflags, qtype-SDC_PROPER_NODES, nnodes)

    else

       ! nodes are given by refining the finest set of nodes

       refine = (nnodes0 - 1) / (nnodes - 1)

       call sdc_qnodes(qnodes0, flags0, qtype, nnodes0)

       qnodes = qnodes0(::refine)
       nodes  = real(qnodes, pfdp)
       nflags = flags0(::refine)

    end if

    call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    if (all(nodes == 0.0d0)) then
       stop 'ERROR: pf_quadrature: invalid SDC nnodes.'
    end if

  end subroutine pf_quadrature

end module pf_mod_quadrature
