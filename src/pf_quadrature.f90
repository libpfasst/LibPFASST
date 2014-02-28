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

    real(c_long_double) :: qnodes0(nnodes0), qnodes(nnodes), dt
    real(pfdp)          :: qmat0(nnodes0-1,nnodes0), smat0(nnodes0-1,nnodes0)
    integer(c_int)      :: flags0(nnodes0)

    integer :: i, r, refine

    qmat = 0
    smat = 0

    if (qtype > SDC_COMPOSITE_NODES) then

       ! nodes are given by repeating the coarsest set of nodes.  note
       ! that in this case nnodes0 corresponds to the coarsest number
       ! of nodes.

       refine = (nnodes - 1) / (nnodes0 - 1)

       call sdc_qnodes(qnodes0, flags0, qtype-SDC_COMPOSITE_NODES, nnodes0)
       call sdc_qmats(qmat0, smat0, qnodes0, qnodes0, flags0, nnodes0, nnodes0)

       dt = 1.q0 / refine
       do i = 1, refine
          r = (i-1)*(nnodes0-1)+1
          qnodes(r:r+nnodes0) = dt * ((i-1) + qnodes0)
          smat(r:r+nnodes0,r:r+nnodes0-1) = smat0 / refine
       end do

       nodes = real(qnodes, pfdp)

    else if (qtype < SDC_COMPOSITE_NODES .and. qtype > SDC_PROPER_NODES) then

       ! nodes are given by proper quadrature rules

       call sdc_qnodes(qnodes, nflags, qtype-SDC_PROPER_NODES, nnodes)
       nodes = real(qnodes, pfdp)

       call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    else

       ! nodes are given by refining the finest set of nodes.  note
       ! that in this case nnodes0 corresponds to the finest number of
       ! nodes.

       refine = (nnodes0 - 1) / (nnodes - 1)

       call sdc_qnodes(qnodes0, flags0, qtype, nnodes0)

       qnodes = qnodes0(::refine)
       nodes  = real(qnodes, pfdp)
       nflags = flags0(::refine)

       call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    end if


    if (all(nodes == 0.0d0)) then
       stop 'ERROR: pf_quadrature: invalid SDC nnodes.'
    end if

  end subroutine pf_quadrature

end module pf_mod_quadrature
