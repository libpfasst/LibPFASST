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

module pf_mod_quadrature
  use pf_mod_dtype
  use sdc_mod_quadrature
  implicit none
contains

  subroutine pf_quadrature(qtype_in, nnodes, nnodes0, nodes, nflags, smat, qmat,LUmat,FEmat,BEmat)
    integer,    intent(in)  :: qtype_in, nnodes, nnodes0
    real(pfdp), intent(out) :: nodes(nnodes)
    real(pfdp), intent(out) :: smat(nnodes-1,nnodes), qmat(nnodes-1,nnodes)
    real(pfdp), intent(out) :: LUmat(nnodes-1,nnodes), FEmat(nnodes-1,nnodes), BEmat(nnodes-1,nnodes)
    integer,    intent(out) :: nflags(nnodes)

    real(c_long_double) :: qnodes0(nnodes0), qnodes(nnodes), dt
    real(pfdp)          :: qmat0(nnodes0-1,nnodes0), smat0(nnodes0-1,nnodes0)
    integer(c_int)      :: flags0(nnodes0)

    integer :: qtype, i, r, refine,m,n
    logical :: composite, proper, no_left
    real(pfdp) :: dsdc(nnodes-1)

    proper    = btest(qtype_in, 8)
    composite = btest(qtype_in, 9)
    no_left   = btest(qtype_in, 10)

    qmat = 0
    smat = 0
    flags0 = 0
    nflags = 0

    qtype = qtype_in
    if (proper)    qtype = qtype - SDC_PROPER_NODES
    if (composite) qtype = qtype - SDC_COMPOSITE_NODES
    if (no_left)   qtype = qtype - SDC_NO_LEFT


    if (composite) then

       ! nodes are given by repeating the coarsest set of nodes.  note
       ! that in this case nnodes0 corresponds to the coarsest number
       ! of nodes.

       refine = (nnodes - 1) / (nnodes0 - 1)

       call sdc_qnodes(qnodes0, flags0, qtype, nnodes0)
       call sdc_qmats(qmat0, smat0, qnodes0, qnodes0, flags0, nnodes0, nnodes0)

       dt = 1.q0 / refine
       do i = 1, refine
          r = (i-1)*(nnodes0-1)+1
          qnodes(r:r+nnodes0) = dt * ((i-1) + qnodes0)
          smat(r:r+nnodes0,r:r+nnodes0-1) = smat0 / refine
       end do

       nodes = real(qnodes, pfdp)

    else if (proper) then

       ! nodes are given by proper quadrature rules

       call sdc_qnodes(qnodes, nflags, qtype, nnodes)
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

       if (no_left) nflags(1) = 0

       call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    end if


    if (all(nodes == 0.0d0)) then
       stop 'ERROR: pf_quadrature: invalid SDC nnodes.'
    end if

    ! make the approximate integration matrices
    BEmat = 0.0_pfdp
    FEmat = 0.0_pfdp
    dsdc = qnodes(2:nnodes) - qnodes(1:nnodes-1)
    do m = 1, nnodes-1
       do n = 1,m
          BEmat(m,n+1) =  dsdc(n)
       end do
       do n = 1,m
          FEmat(m,n) =  dsdc(n)
       end do
    end do
    ! Get the LU
    call myLUq(qmat,LUmat,nnodes,1)
!    print *,LUmat

  end subroutine pf_quadrature
  subroutine pf_myLUexp(A,L,U,Nnodes,scaleLU)
    real(pfdp),       intent(in)    :: A(Nnodes,Nnodes)
    real(pfdp),      intent(inout)  :: L(Nnodes,Nnodes)
    real(pfdp),     intent(inout)   :: U(Nnodes,Nnodes)
    integer,        intent (in)     :: Nnodes
    integer,        intent (in)     :: scaleLU
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

    U=transpose(U)
    !  Now scale the columns of U to match the sum of A
    if (scaleLU .eq. 1) then
       do j=1,Nnodes
          c = sum(U(j,:))
          if (c /=  0.0) then
             U(j,:)=U(j,:)*sum(A(j,:))/c
          end if
       end do
    end if

  end subroutine pf_myLUexp
  subroutine myLUq(Q,Qtil,Nnodes,fillq)
    real(pfdp),       intent(in)    :: Q(Nnodes-1,Nnodes)
    real(pfdp),     intent(inout)   :: Qtil(Nnodes-1,Nnodes)
    integer,        intent (in)     :: Nnodes
    integer,        intent (in)     :: fillq

    ! Return the Qtil=U^T where U is the LU decomposition of Q without pivoting
    ! if fillq is positive, then the first row of Qtil is filled to make
    ! the matrix consistent

    integer :: i,j,N
    real(pfdp) :: c
    real(pfdp)  :: U(Nnodes-1,Nnodes-1)
    real(pfdp)  :: L(Nnodes-1,Nnodes-1)
    L = 0.0_pfdp
    U = 0.0_pfdp
    N = Nnodes-1
    U=transpose(Q(1:Nnodes-1,2:Nnodes))

    do i = 1,N
       if (U(i,i) /= 0.0) then
          do j=i+1,N
             c = U(j,i)/U(i,i)
             U(j,i:N)=U(j,i:N)-c*U(i,i:N)
             L(j,i)=c
          end do
       end if
       L(i,i) = 1.0_pfdp
    end do

    !  Check
    if (maxval(abs(matmul(L,U)-transpose(Q(1:Nnodes-1,2:Nnodes)))) > 1.0d-12)  print *,'LU error'
    
    Qtil = 0.0_pfdp
    Qtil(1:Nnodes-1,2:Nnodes)=transpose(U)
    !  Now scale the columns of U to match the sum of A
    if (fillq .eq. 1) then
       do i=1,Nnodes-1
          Qtil(i,1)=sum(Q(i,1:Nnodes))-sum(Qtil(i,2:Nnodes))
       end do
    end if

  end subroutine myLUq

end module pf_mod_quadrature
