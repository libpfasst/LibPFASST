!! Quadrature matrices and accompanying routines
!
! This file is part of LIBPFASST.
!
!
!> Module to create quadrature matrices and accompanying routines
module pf_mod_quadrature
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  integer,  parameter :: qp = c_long_double   
  integer,  parameter :: dp = c_double
  real(qp), parameter :: eps = 1.0e-23_qp
  
  private :: qsort_partition

  interface poly_eval
     module procedure poly_eval
     module procedure poly_eval_complex
  end interface


  
contains
  !>  
  subroutine pf_init_sdcmats(SDCmats,qtype,nnodes,nnodes0,nflags)
    type(pf_sdcmats_t), intent(inout) :: SDCmats
     integer, intent(in) :: nnodes,nnodes0
     integer, intent(in) :: qtype
     integer,    intent(inout) :: nflags(nnodes)
     integer :: ierr

     allocate(SDCmats%qmat(nnodes-1,nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error qmat"
     allocate(SDCmats%qmatFE(nnodes-1,nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error qmatFE"
     allocate(SDCmats%qmatBE(nnodes-1,nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error qmatBE"
     allocate(SDCmats%qmatTrap(nnodes-1,nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error qmatBE"
     allocate(SDCmats%qmatVer(nnodes-1,nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error qmatBE"
     allocate(SDCmats%qmatLU(nnodes-1,nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error qmatLU"
     allocate(SDCmats%s0mat(nnodes-1,nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error s0mat"
     allocate(SDCmats%qnodes(nnodes),stat=ierr)
     if (ierr /= 0) stop "allocate error qnodes"

     call pf_quadrature(qtype, nnodes, nnodes0, &
          SDCmats%qnodes, nflags, SDCmats%s0mat, SDCmats%qmat,SDCmats%qmatFE,SDCmats%qmatBE)

     call myLUq(SDCmats%qmat,SDCmats%qmatLU,nnodes,0)
     SDCmats%qmatTrap=0.5_pfdp*(SDCmats%qmatFE+SDCmats%qmatBE)

  end subroutine pf_init_sdcmats
  
  subroutine pf_destroy_sdcmats(SDCmats)
    type(pf_sdcmats_t), intent(inout) :: SDCmats


     deallocate(SDCmats%qmat)  
     deallocate(SDCmats%qmatFE)  
     deallocate(SDCmats%qmatBE)  
     deallocate(SDCmats%qmatLU)  
     deallocate(SDCmats%qmatTrap)  
     deallocate(SDCmats%qmatVer)  
     deallocate(SDCmats%s0mat)  
     deallocate(SDCmats%qnodes)       
  end subroutine pf_destroy_sdcmats
  
    
  !>  Routine to compute the LU decomposition of spectral integration matrix
  subroutine myLUq(Q,QLU,Nnodes,fillq)
    integer,        intent (in)     :: Nnodes
    real(pfdp),     intent(in)      :: Q(Nnodes-1,Nnodes)
    real(pfdp),     intent(inout)   :: QLU(Nnodes-1,Nnodes)
    integer,        intent (in)     :: fillq

    ! Return the QLU=U^T where U is the LU decomposition of Q without pivoting
    ! if fillq is positive, then the first row of QLU is filled to make
    ! the matrix consistent

    integer :: i,j,N
    real(pfdp) :: c
    real(pfdp)  :: U(Nnodes-1,Nnodes-1)
    real(pfdp)  :: L(Nnodes-1,Nnodes-1)
    real(pfdp)  :: LUerror
    
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
    LUerror = maxval(abs(matmul(L,U)-transpose(Q(1:Nnodes-1,2:Nnodes))))
    if (LUerror .gt. 1e-14)  then
       print *,'error in LU too high',LUerror
       stop
    end if

    QLU = 0.0_pfdp
    QLU(1:Nnodes-1,2:Nnodes)=transpose(U)
    !  Now scale the columns of U to match the sum of A
    if (fillq .eq. 1) then
       do j=1,Nnodes-1
          QLU(j,1)=sum(Q(j,1:Nnodes))-sum(U(j,1:Nnodes-1))
       end do
    end if

  end subroutine myLUq

  !>  Subroutine to create quadrature matrices
  subroutine pf_quadrature(qtype_in, nnodes, nnodes0, nodes, nflags, smat, qmat,qmatFE,qmatBE)
    integer,    intent(in)  :: qtype_in, nnodes, nnodes0
    real(pfdp), intent(out) :: nodes(nnodes)
    real(pfdp), intent(out) :: smat(nnodes-1,nnodes), qmat(nnodes-1,nnodes)
    real(pfdp), intent(out) :: qmatFE(nnodes-1,nnodes), qmatBE(nnodes-1,nnodes)
    integer,    intent(out) :: nflags(nnodes)

    real(pfqp) :: qnodes0(nnodes0), qnodes(nnodes), dt
    real(pfdp)          :: qmat0(nnodes0-1,nnodes0), smat0(nnodes0-1,nnodes0)
    integer             :: flags0(nnodes0)

    integer :: qtype, i, r, refine,m,n
    logical :: composite, proper, no_left

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


    !  Make forward and backward Euler matrices
    qmatFE=0.0_pfdp
    qmatBE=0.0_pfdp
    do m = 1, nnodes-1
       do n = 1,m
          qmatBE(m,n+1) =  qnodes(n+1)-qnodes(n)
       end do
    end do
    ! Explicit matrix
    do m = 1, nnodes-1
       do n = 1,m
          qmatFE(m,n)   =  qnodes(n+1)-qnodes(n)
       end do
    end do


    if (all(nodes == 0.0d0)) then
       stop 'ERROR: pf_quadrature: invalid SDC nnodes.'
    end if

  end subroutine pf_quadrature

  !>  Function to decide if the restriction of the nodes is pointwise, e.g. coarse nodes are every other fine node
  logical function not_proper(flags, node)
    integer , intent(in) :: flags(:)
    integer,        intent(in) :: node

    not_proper = .not. btest(flags(node), 0)
  end function not_proper



  !> Subroutine to compute high precision quadrature nodes.
  subroutine sdc_qnodes(qnodes, flags, qtype, nnodes) 
    integer ,       intent(in), value  :: nnodes          !!  Number of nodes
    integer ,       intent(in), value  :: qtype           !!  Type of nodes (see pf_dtype)
    real(pfqp),  intent(out)        :: qnodes(nnodes)  !!  The computed quadrature nodes
    integer ,       intent(out)        :: flags(nnodes)   !!

    integer :: j, degree
    real(pfqp), allocatable :: roots(:)
    real(pfqp), allocatable :: coeffs(:), coeffs2(:)

    real(pfqp), parameter :: pi = 3.141592653589793115997963468544185161590576171875_pfdp

    flags = 0

    select case(qtype)

    case (SDC_GAUSS_LEGENDRE)

       degree = nnodes-2
       allocate(roots(degree))
       allocate(coeffs(degree+1))

       call poly_legendre(coeffs, degree)
       call poly_roots(roots, coeffs, degree)

       qnodes(1) = 0.0_pfqp
       qnodes(2:nnodes-1) = 0.5_pfqp * (1.0_pfqp + roots)
       qnodes(nnodes) = 1.0_pfqp

       deallocate(coeffs)
       deallocate(roots)

       do j = 2, nnodes-1
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_GAUSS_LOBATTO)
       
       degree = nnodes - 1
       allocate(roots(degree-1))
       allocate(coeffs(degree+1))

       call poly_legendre(coeffs, degree)
       call poly_diff(coeffs, degree)
       call poly_roots(roots, coeffs(:degree), degree-1)

       qnodes(1)          = 0.0_pfdp
       qnodes(2:nnodes-1) = 0.5_pfdp * (1.0_pfdp + roots)
       qnodes(nnodes)     = 1.0_pfdp

       deallocate(coeffs)
       deallocate(roots)

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_GAUSS_RADAU)

       degree = nnodes - 1
       allocate(roots(degree))
       allocate(coeffs(degree+1))
       allocate(coeffs2(degree))


       call poly_legendre(coeffs, degree)
       call poly_legendre(coeffs2, degree-1)
       coeffs(:degree) = coeffs(:degree) + coeffs2
       call poly_roots(roots, coeffs, degree)

       qnodes(1)      = 0.0_pfdp
       do j = 2, nnodes-1
          qnodes(j) = 0.5_pfdp * (1.0_pfdp - roots(nnodes+1-j))
       end do
       qnodes(nnodes) = 1.0_pfdp

       deallocate(coeffs2)
       deallocate(coeffs)
       deallocate(roots)

       do j = 2, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_CLENSHAW_CURTIS)

       do j = 0, nnodes-1
          qnodes(j+1) = 0.5_pfqp * (1.0_pfdp - cos(j * pi / (nnodes-1)))
       end do

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_UNIFORM)

       do j = 0, nnodes-1
          qnodes(j+1) = j * (1.0_pfqp / (nnodes-1))
       end do

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case default
       print *,'qtype = ',qtype
       stop "ERROR: Invalid qtype in sdc_quadrature.f90."

    end select

  end subroutine sdc_qnodes

  !>  Subroutine to compute the quadrature matrices 
  subroutine sdc_qmats(qmat, smat, dst, src, flags, ndst, nsrc)
    integer ,  intent(in), value  :: ndst   !!  Number of destination points
    integer ,   intent(in), value  :: nsrc   !!  Number of source points
    real(pfqp), intent(in)  :: dst(ndst)     !!  Destination points
    real(pfqp), intent(in)  :: src(nsrc)     !!  Source points
    real(pfdp),      intent(out) :: qmat(ndst-1, nsrc)  !!  O to dst quadrature weights
    real(pfdp),      intent(out) :: smat(ndst-1, nsrc)  !! dst(m) to dst(m+1) quadrature weights
    integer ,      intent(in)  :: flags(nsrc)     

    integer  :: i, j, m
    real(pfqp) :: q, s, den, p(0:nsrc)

    qmat = 0.0_pfqp
    smat = 0.0_pfqp

    ! construct qmat and smat
    do i = 1, nsrc

       if (not_proper(flags, i)) cycle

       ! construct interpolating polynomial coefficients
       p    = 0.0_pfdp
       p(0) = 1.0_pfdp
       do m = 1, nsrc
          if (not_proper(flags, m) .or. m == i) cycle
          p = eoshift(p, -1) - src(m) * p
       end do

       den = poly_eval(p, nsrc, src(i))

       call poly_int(p, nsrc)

       ! evaluate integrals
       do j = 2, ndst
          q = poly_eval(p, nsrc, dst(j)) - poly_eval(p, nsrc,   0.0_pfqp)
          s = poly_eval(p, nsrc, dst(j)) - poly_eval(p, nsrc, dst(j-1))

          qmat(j-1,i) = real(q / den, pfqp)
          smat(j-1,i) = real(s / den, pfqp)
       end do
    end do

  end subroutine sdc_qmats


  !> Polynomial manipulation routines.
  !!
  !! A polynomial p
  !!
  !!   p(x) = a_n x^n + ... + a_2 x^2 + a_1 x + a_0
  !!
  !! is stored as a Fortran array p(0:n) according to
  !!
  !!   p = [ a_0, a_1, ..., a_n ].
  !!
  
  !> Function to evaluate real polynomial
  real(pfqp) function poly_eval(p, n, x) result(v) 
    integer, intent(in), value :: n
    real(pfqp),       intent(in)        :: p(0:n), x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function

  !> Function to evaluate complex polynomial
  complex(pfdp) function poly_eval_complex(p, n, x) result(v)
    integer, intent(in), value :: n
    real(pfqp),       intent(in)        :: p(0:n)
    complex(pfqp),    intent(in)        :: x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function



  !> Subroutine to differentiate polynomial (in place)
  subroutine poly_diff(p, n) 
    integer, intent(in),   value :: n
    real(pfqp),       intent(inout) :: p(0:n)

    integer  :: j
    real(pfdp) :: pp(0:n)

    pp = 0.0_pfqp

    do j = 1, n
       pp(j-1) = j * p(j)
    end do

    p = pp
  end subroutine poly_diff


  !> Subroutine to integrate polynomial (in place)
  subroutine poly_int(p, n) 
    integer, intent(in),   value :: n
    real(pfqp),       intent(inout) :: p(0:n)

    integer  :: j
    real(pfqp) :: pp(0:n)

    pp = 0.0_pfqp

    do j = 0, n-1
       pp(j+1) = p(j) / (j+1)
    end do

    p = pp
  end subroutine poly_int


  !> Subroutine to compute Legendre polynomial coefficients using Bonnet's recursion formula.
  subroutine poly_legendre(p, n) 
    integer, intent(in), value :: n
    real(pfqp),       intent(out)       :: p(0:n)

    real(pfqp), dimension(0:n) :: p0, p1, p2
    integer :: j, m

    if (n == 0) then
       p = [ 1.0_pfqp ]
       return
    end if

    if (n == 1) then
       p = [ 0.0_pfqp, 1.0_pfqp ]
       return
    end if

    p0 = 0.0_pfqp;
    p1 = 0.0_pfqp;
    p2 = 0.0_pfqp

    p0(0) = 1.0_pfqp
    p1(1) = 1.0_pfqp

    ! (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
    do m = 1, n-1
       do j = 1, n
          p2(j) = ( (2*m + 1) * p1(j-1) - m * p0(j) ) / (m + 1)
       end do
       p2(0) = - m * p0(0) / (m + 1)

       p0 = p1
       p1 = p2
    end do

    p = p2
  end subroutine poly_legendre

  !> Subroutine to compute polynomial roots using the Durand-Kerner algorithm.
  !! The roots are assumed to be real.
  subroutine poly_roots(roots, p0, n) 
    integer,  intent(in), value   :: n
    real(pfqp),        intent(out)  :: roots(n)
    real(pfqp),        intent(in)   :: p0(0:n)

    integer     :: i, j, k
    complex(pfqp) :: num, den, z0(n), z1(n)
    real(pfqp)    :: p(0:n)
    
    real(pfqp) ::  eps 

    eps = epsilon(1.0_pfqp)*100.0_pfqp
    p = p0 / p0(n)

    ! initial guess
    do i = 1, n
       z0(i) = (0.4_pfqp, 0.9_pfqp)**i
    end do

    ! durand-kerner-weierstrass iterations
    z1 = z0
    do k = 1, 100
       do i = 1, n

          ! evaluate poly at z0(i)
          num = poly_eval_complex(p, n, z0(i))

          ! evaluate denominator
          den = 1.0_pfqp
          do j = 1, n
             if (j == i) cycle
             den = den * (z0(i) - z0(j))
          end do

          ! update
          z0(i) = z0(i) - num / den
       end do

       ! converged?
       if (sum(abs(z0 - z1)) < eps) exit

       z1 = z0
    end do

    roots = real(z0)
    where (abs(roots) < eps) roots = 0.0_pfqp
    call qsort(roots)

  end subroutine poly_roots

  !> Subroutine to sort (inplace) using the quick sort algorithm.
  !> Adapted from http://www.fortran.com/qsort_c.f95.
  recursive subroutine qsort(a)
    real(pfqp), intent(inout) :: a(:)
    integer :: iq

    if (size(a) > 1) then
       call qsort_partition(a, iq)
       call qsort(a(:iq-1))
       call qsort(a(iq:))
    end if
  end subroutine qsort

  subroutine qsort_partition(a, marker)
    real(pfqp), intent(inout) :: a(:)
    integer,  intent(out)   :: marker

    integer  :: i, j
    real(pfqp) :: temp, x

    x = a(1)
    i = 0
    j = size(a) + 1

    do
       j = j-1
       do
          if (a(j) <= x) exit
          j = j-1
       end do

       i = i+1
       do
          if (a(i) >= x) exit
          i = i+1
       end do

       if (i < j) then
          temp = a(i)
          a(i) = a(j)
          a(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
  end subroutine qsort_partition

end module pf_mod_quadrature






