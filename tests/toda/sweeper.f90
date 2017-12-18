!-------------------------------------------------------------------------------
! Copyright (c) 2017, Brandon Krull.  All rights reserved.
!--------------------------------------------------------------------------------
! MODULE: sweeper
! !> @author
!> Brandon Krull, Berkeley Lab
!
! Description:
!> This module contains sweeper and related functionality
module sweeper
  use pf_mod_dtype
  use pf_mod_magnus_picard
  use factory

  implicit none

  real(pfdp), parameter :: &
       pi = 3.141592653589793_pfdp, &
       two_pi = 6.2831853071795862_pfdp

  complex(pfdp), parameter :: &
       z0 = (0.0_pfdp, 0.0_pfdp), &
       z1 = (1.0_pfdp, 0.0_pfdp), &
       zm1 = (-1.0_pfdp, 0.0_pfdp), &
       zi = (0.0_pfdp, 1.0_pfdp), &
       zmi = (0.0_pfdp, -1.0_pfdp)

  external :: zgemm

  type, extends(pf_user_level_t) :: magpicard_context
   contains
     procedure :: restrict => restrict
     procedure :: interpolate => interpolate
  end type magpicard_context

  type, extends(pf_magpicard_t) :: magpicard_sweeper_t
     integer :: dim, exp_iterations
     logical :: debug
     complex(pfdp), allocatable :: commutator(:,:)
   contains
     procedure :: f_eval => compute_B
     procedure :: compute_omega
     procedure :: compute_time_ev_ops
     procedure :: propagate_solution
     procedure :: destroy => destroy_magpicard_sweeper
  end type magpicard_sweeper_t

contains

  function cast_as_magpicard_sweeper(sweeper) result(magpicard)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(magpicard_sweeper_t), pointer :: magpicard

    select type(sweeper)
    type is (magpicard_sweeper_t)
       magpicard => sweeper
    class default
       print*, 'invalid sweeper class'
       stop
    end select

  end function cast_as_magpicard_sweeper

  subroutine initialize_magpicard_sweeper(this, level, qtype, debug, shape, nvars)
    use probin, only: nprob, basis, molecule, magnus_order, nparticles, dt
    class(pf_sweeper_t), intent(inout) :: this
    integer, intent(in) :: level, qtype
    logical, intent(in) :: debug
    integer, intent(inout) :: shape(:)
    integer, intent(inout) :: nvars

    class(magpicard_sweeper_t), pointer :: magpicard !< context data containing integrals, etc

    ! print*, 'Creating context for basis ', basis(level)
    magpicard => cast_as_magpicard_sweeper(this)

    magpicard%qtype = qtype
    magpicard%dt = dt
    magpicard%magnus_order = magnus_order(level)
    magpicard%exp_iterations = 0
    magpicard%debug = debug
    magpicard%dim = nparticles

    allocate(magpicard%commutator(nparticles, nparticles))

    magpicard%commutator = z0

    ! inouts necessary for the base class structure
    shape(:) = magpicard%dim
    nvars = magpicard%dim * magpicard%dim * 2

    nullify(magpicard)
  end subroutine initialize_magpicard_sweeper

  !> Set initial condition.
  subroutine initial(L)
    class(pf_encap_t), intent(inout) :: L
    class(zndarray), pointer :: L_p

    L_p => cast_as_zndarray(L)

    call exact(L_p, 0.0_pfdp)

    nullify(L_p)
  end subroutine initial

  !> Set exact solution
  !! Exact solution mat(u(t)) = mat(U(t,t0))*mat(u(t0))
  !! U(t,t0) = e^(-iH(t-t0)),
  !! project U into basis of eigenstates of H
  !! aka diagonalize H, then
  !! U(t,t0) = sum_j(|j>e^(-ie_j(t-t0))<j|)
  subroutine exact(L, t)
    use probin, only: alpha, beta, vab, toda_periodic

    type(zndarray), intent(inout) :: L
    real(pfdp), intent(in)  :: t

    complex(pfdp), allocatable :: time_ev_op(:,:), time_ev_op_dagger(:,:), f(:,:)
    complex(pfdp), allocatable :: u(:,:), uinit(:,:) !> definition of initial state
    real(pfdp), allocatable :: q(:), p(:)
    real(pfdp) :: tol=1d-15, alpha_toda
    integer :: dim, i

    dim = L%dim
    allocate(time_ev_op(dim, dim), &
          time_ev_op_dagger(dim, dim), &
          f(dim, dim), &
          u(dim, dim), &
          uinit(dim, dim), q(dim), p(dim))

    u = z0
    uinit = z0
    time_ev_op = z0
    time_ev_op_dagger = z0
    f = z0
    q = 0.0_pfdp
    p = 0.0_pfdp

       ! See Zanna thesis on "On the Numerical Solution of Isospectral Flows"
    if (toda_periodic .eqv. .true.) then
       q = 0
    else
       do i = 1, dim
           q(i) = i - dim/2 - 1 ! everyone assumes a unique position on the number line
       enddo
    endif

    do i = 1, dim
       if (i <= dim/2) then
          p(i) = 4
       elseif (i == dim/2+1) then
          p(i) = 0
       else
          p(i) = -4
       endif
       u(i,i) = p(i) ! everyone gets a different initial momentum
    enddo

    do i = 1, dim-1
       alpha_toda = exp(-0.5_pfdp * (q(i+1) - q(i)))
       if (i == 1 .and. toda_periodic .eqv. .true.) then
          u(i, dim) = alpha_toda
          u(dim, i) = alpha_toda
       endif
       u(i, i+1) = alpha_toda
       u(i+1, i) = alpha_toda
    enddo

    u(:,:) = 0.5_pfdp * u(:,:)
    L%array = u

    deallocate(time_ev_op, time_ev_op_dagger, f, u, uinit, q, p)

  end subroutine exact

  subroutine compute_B(this, y, t, level, f)
      use probin, only: toda_periodic

      class(magpicard_sweeper_t), intent(inout) :: this
      class(pf_encap_t), intent(inout) :: y ! prev solution
      class(pf_encap_t), intent(inout) :: f ! output RHS
      real(pfdp), intent(in) :: t
      integer, intent(in) :: level

      type(zndarray), pointer :: L, B
      integer :: i

      L => cast_as_zndarray(y)
      B => cast_as_zndarray(f)

      do i = 1, this%dim
         B%array(i,i) = 0.0_pfdp
      enddo

      do i = 1, this%dim-1
         B%array(i, i+1) = -1.0_pfdp * L%array(i, i+1)
         B%array(i+1, i) = L%array(i, i+1)
      enddo

      if (toda_periodic .eqv. .true.) then
         B%array(1, this%dim) = L%array(1, this%dim)
         B%array(this%dim, 1) = -1.0_pfdp * L%array(this%dim, 1)
      endif

      nullify(L, B)

  end subroutine compute_B

 !> Compute the matrix $\Omega$ to be exponentiated from AO-Fock matrix
 !! \Omega = \sum_i omega_i
 !! omega coming in already has omega_1 on it (from magpicard_integrate (lev%S(m)))
 !! this subroutine really just adds additional omega_i
 !! \f[\Omega_1(t+dt,t) = -i\int_{t}^{t+dt} F(\tau)d\tau\f]
 !! this can be computed using any number of approximate numerical methods
 !! Trapezoid/midpoint rule for low-order accuracy, Simpson's for higher-order
 subroutine compute_omega(this, omega, integrals, f, nodes, qmat, dt, this_node, coefs)
   class(magpicard_sweeper_t), intent(inout) :: this
   class(pf_encap_t), intent(inout) :: omega, integrals(:), f(:,:)
   real(pfdp), intent(in) :: coefs(:,:), nodes(:), qmat(:,:), dt
   integer, intent(in) :: this_node

   class(zndarray), pointer :: omega_p, ints
   integer :: dim

   dim = this%dim
   omega_p => cast_as_zndarray(omega)
   ints => cast_as_zndarray(integrals(this_node))

   omega_p%array = ints%array
   if (this%magnus_order > 1) then
      call add_single_commutator_terms(omega_p%array, f, coefs(:,1), dim, this%commutator)

      if (this%magnus_order > 2) then
         call add_double_commutator_terms(omega_p%array, f, coefs(:,2), dim, this%commutator)
         call add_triple_commutator_terms(omega_p%array, f, nodes, this_node, qmat, dt, &
              coefs(1,3), dim, this%commutator)
      endif
   endif


   nullify(omega_p, ints)
 end subroutine compute_omega

 subroutine add_single_commutator_terms(omega, f, coefs, dim, commutator)
   class(pf_encap_t), intent(inout) :: f(:,:)
   complex(pfdp), intent(inout) :: omega(dim, dim), commutator(dim,dim)
   real(pfdp), intent(in) :: coefs(:)
   integer, intent(in) :: dim

   class(zndarray), pointer :: f1, f2, f3
   complex(pfdp), allocatable :: tmp(:,:)
   integer :: i, nnodes

   ! print*, 'omega before single term = ', (real(omega(i,i+2), pfdp), i=1,dim-5)
   allocate(tmp(dim,dim))
   tmp = z0
   nnodes = size(f)
   if (nnodes > 3) then !hopefully this only corresponds to 3gauss nodes aka 5 total nodes
      f1 => cast_as_zndarray(f(2,1))
      f2 => cast_as_zndarray(f(3,1))
      f3 => cast_as_zndarray(f(4,1))
   else
      f1 => cast_as_zndarray(f(1,1))
      f2 => cast_as_zndarray(f(2,1))
      f3 => cast_as_zndarray(f(3,1))
   endif

   call compute_commutator(f1%array, f2%array, dim, commutator)
   tmp = tmp + coefs(1) * commutator

   call compute_commutator(f1%array, f3%array, dim, commutator)
   tmp = tmp + coefs(2) * commutator

   call compute_commutator(f2%array, f3%array, dim, commutator)
   tmp = tmp + coefs(3) * commutator
   omega = omega + tmp

   ! print*, 'single commutator term = ', (real(tmp(i,i+2), pfdp), i=1,dim-5)
   ! print*, 'omega after single term = ', (real(omega(i,i+2), pfdp), i=1,dim-5)
   print*, 'coefs = ', coefs
   deallocate(tmp)
   nullify(f1, f2, f3)
 end subroutine add_single_commutator_terms

 subroutine add_double_commutator_terms(omega, f, coef, dim, commutator)
   class(pf_encap_t), intent(inout) :: f(:,:)
   complex(pfdp), intent(inout) :: omega(dim,dim), commutator(dim,dim)
   real(pfdp), intent(in) :: coef(:)
   integer, intent(in) :: dim

   class(zndarray), pointer :: f1, f2, f3
   complex(pfdp), allocatable :: tmp(:,:,:), tmp_comm(:,:)
   integer :: i, nnodes

   allocate(tmp(4,dim,dim), tmp_comm(dim,dim))

   tmp = z0
   tmp_comm = z0
   nnodes = size(f)

   do i = 1, 3
      f1 => cast_as_zndarray(f(i+1,1))
      tmp(1,:,:) = tmp(1,:,:) + coef(i)   * f1%array
      tmp(2,:,:) = tmp(2,:,:) + coef(i+3) * f1%array
      tmp(3,:,:) = tmp(3,:,:) + coef(i+6) * f1%array
   end do

   f1 => cast_as_zndarray(f(2,1))
   f2 => cast_as_zndarray(f(3,1))
   f3 => cast_as_zndarray(f(4,1))
   call compute_commutator(f1%array, f2%array, dim, tmp_comm)
   call compute_commutator(tmp(1,:,:), tmp_comm, dim, commutator)
   tmp(4,:,:) = tmp(4,:,:) + commutator

   call compute_commutator(f1%array, f3%array, dim, tmp_comm)
   call compute_commutator(tmp(2,:,:), tmp_comm, dim, commutator)
   tmp(4,:,:) = tmp(4,:,:) + commutator

   call compute_commutator(f2%array, f3%array, dim, tmp_comm)
   call compute_commutator(tmp(3,:,:), tmp_comm, dim, commutator)
   tmp(4,:,:) = tmp(4,:,:) + commutator

   omega = omega + tmp(4,:,:)

   print*, 'double commutator term = ', (real(tmp(4,i,i+1)), i=1,dim-5)
   print*, 'omega after double term = ', (real(omega(i,i+1)), i=1,dim-5)
   ! print*, 'double commutator term = ', real(tmp(4,:,:))
   ! print*, 'omega after adding double terms', omega(1:3,11)
   deallocate(tmp, tmp_comm)
   nullify(f1, f2, f3)
 end subroutine add_double_commutator_terms

 subroutine add_triple_commutator_terms(omega, f, nodes, this_node, qmat, dt, coef, dim, commutator)
   class(pf_encap_t), intent(inout) :: f(:,:)
   complex(pfdp), intent(inout) :: omega(dim,dim), commutator(dim,dim)
   real(pfdp), intent(in) :: coef, nodes(:), qmat(:,:), dt
   integer, intent(in) :: dim, this_node

   class(zndarray), pointer :: f1
   complex(pfdp), allocatable :: tmp(:,:), a(:,:,:)
   real(pfdp) :: time_scaler
   integer :: i, j, nnodes

   allocate(tmp(dim,dim), a(3,dim,dim))
   tmp = z0
   a = z0

   nnodes = size(nodes)

   do j = 2, nnodes-1
       f1 => cast_as_zndarray(f(j,1))
       time_scaler = (nodes(j)-0.5_pfdp)
       tmp = dt * qmat(this_node,j) * f1%array
       a(1,:,:) = a(1,:,:) + tmp
       a(2,:,:) = a(2,:,:) + tmp * time_scaler
       a(3,:,:) = a(3,:,:) + tmp * time_scaler**2
   enddo

   print*, 'triple a0 = ', (real(a(1,i,i+1)), i=1,dim-5)
   print*, 'triple a1 = ', (real(a(2,i,i+1)), i=1,dim-5)

   tmp = a(2,:,:)
   call compute_commutator(a(1,:,:), tmp, dim, commutator)
   a(3,:,:) = commutator
   call compute_commutator(a(1,:,:), a(3,:,:), dim, commutator)
   a(3,:,:) = commutator
   call compute_commutator(a(1,:,:), a(3,:,:), dim, commutator)
   ! do i = 1, 3
   !    call compute_commutator(a(1,:,:), tmp, dim, commutator)
   !    tmp = commutator
   ! end do

   omega = omega + commutator / 60.d0

   ! print*, 'triple commutator coef = ', coef
   ! print*, 'triple commutator term = ', (real(commutator(i,i+2)), i=1,dim-5)
   ! print*, 'omega after triple term = ', (real(omega(i,i+2)), i=1,dim-5)
   ! print*, 'triple commutator term = ', (coef(1)*commutator(i,i+1), i=1,nnodes-1)
   ! print*, 'triple commutator term = ', real(tmp(4,:,:))
   ! print*, 'omega after adding triple terms', omega(1:3,11)
   deallocate(tmp, a)
   nullify(f1)
 end subroutine add_triple_commutator_terms

 subroutine compute_commutator(a, b, dim, commutator)
   complex(pfdp), intent(in) :: a(dim,dim), b(dim,dim)
   integer, intent(in) :: dim
   complex(pfdp), intent(inout) :: commutator(dim,dim)

   call zgemm('n', 'n', dim, dim, dim, &
        z1, b, dim, &
        a, dim, &
        z0, commutator, dim)

   call zgemm('n', 'n', dim, dim, dim, &
        z1, a, dim, &
        b, dim, &
        zm1, commutator, dim)
 end subroutine compute_commutator

 !> Compute matrix exponential
 !! u = exp{omega}
 !! u_dagger = exp{-omega}
 subroutine compute_time_ev_ops(this, time_ev_op, omega, level)
   use probin, only: exptol
   class(magpicard_sweeper_t), intent(inout) :: this
   class(pf_encap_t), intent(inout) :: time_ev_op, omega
   integer, intent(in) :: level

   class(zndarray), pointer :: omega_p
   class(zndarray), pointer :: time_ev_op_p

   integer :: dim

   omega_p => cast_as_zndarray(omega)
   time_ev_op_p => cast_as_zndarray(time_ev_op)

   dim = omega_p%dim
   if(this%debug) then
      print*, '----------Inside compute_time_ev_ops----------'
      print*, 'dim omega= ', omega_p%dim
      ! print*, 'shape omega = ', shape(omega_p%array)
      print*, 'omega = ', omega_p%array
   end if

   call c8mat_expm1(dim, omega_p%array, time_ev_op_p%array)

   if (time_ev_op_p%norm()*0.0 /= 0.0) then
      print*, omega_p%array
      print*, time_ev_op_p%array
      stop
   endif

   nullify(omega_p, time_ev_op_p)
 end subroutine compute_time_ev_ops

 !> Computes the P_t = U*P_t0*U^dagger
 subroutine propagate_solution(this, sol_t0, sol_tn, u)
   use probin, only: nprob
   class(magpicard_sweeper_t), intent(inout) :: this
   class(pf_encap_t), intent(inout) :: sol_t0
   class(pf_encap_t), intent(inout) :: sol_tn
   class(pf_encap_t), intent(inout) :: u !< Time-evolution operator
   integer :: dim !< size of dimensions of P, U
   class(zndarray), pointer :: sol_t0_p, sol_tn_p, u_p
   complex(pfdp), allocatable :: tmp(:,:)

   sol_t0_p => cast_as_zndarray(sol_t0)
   sol_tn_p => cast_as_zndarray(sol_tn)
   u_p => cast_as_zndarray(u)

   dim = sol_t0_p%dim
   allocate(tmp(dim, dim))

   if (nprob < 10) then
      call zgemm('n', 'n', dim, dim, dim, &
           z1, u_p%array, dim, &
           sol_t0_p%array, dim, &
           z0, tmp, dim)

      call zgemm('n', 'c', dim, dim, dim, &
           z1, tmp, dim, &
           u_p%array, dim, &
           z0, sol_tn_p%array, dim)
   else
      call zgemm('n', 'n', dim, dim, dim, &
           z1, u_p%array, dim, &
           sol_t0_p%array, dim, &
           z0, sol_tn_p%array, dim)
   endif

   deallocate(tmp)
   nullify(sol_t0_p, sol_tn_p, u_p)
 end subroutine propagate_solution

 function compute_inf_norm(matrix, n) result(norm)
   integer, intent(in) :: n
   complex(pfdp), intent(in) :: matrix(n,n)

   integer :: i, j
   real(pfdp) :: norm, tmp(n*n)

   i = 0
   j = 0
   norm = 0d0
   tmp = 0d0

   do j = 1, n
      do i = 1, n
         tmp(i) = tmp(i) + abs(matrix(i,j))
      enddo
   enddo

   do i = 1, n
      norm = max(norm, tmp(i))
   enddo

 end function compute_inf_norm

 !> array of ctx data deallocation
 subroutine destroy_magpicard_sweeper(this, lev)
   class(magpicard_sweeper_t), intent(inout) :: this
   class(pf_level_t),   intent(inout) :: lev
   integer :: io

   call this%magpicard_destroy(lev)

 end subroutine destroy_magpicard_sweeper

 subroutine restrict(this, levelF, levelG, qF, qG, t)
   class(magpicard_context), intent(inout) :: this
   class(pf_level_t), intent(inout) :: levelF
   class(pf_level_t), intent(inout) :: levelG
   class(pf_encap_t), intent(inout) :: qF, qG
   real(pfdp),        intent(in   ) :: t

   class(zndarray), pointer :: f, g
   f => cast_as_zndarray(qF)
   g => cast_as_zndarray(qG)

   g%array = f%array
 end subroutine restrict

 subroutine interpolate(this, levelF, levelG, qF, qG, t)
   class(magpicard_context), intent(inout) :: this
   class(pf_level_t), intent(inout) :: levelF
   class(pf_level_t), intent(inout) :: levelG
   class(pf_encap_t), intent(inout) :: qF, qG
   real(pfdp),        intent(in   ) :: t

   class(zndarray), pointer :: f, g
   f => cast_as_zndarray(qF)
   g => cast_as_zndarray(qG)

   f%array = g%array
 end subroutine interpolate

 subroutine initialize_as_identity_real(matrix)
   real(pfdp), intent(inout) :: matrix(:,:)

   integer :: i, dim, shp(2)

   shp = shape(matrix)
   dim = shp(1)

   matrix = 0.0_pfdp
   forall (i=1:dim) matrix(i,i) = 1.0_pfdp
 end subroutine initialize_as_identity_real

 subroutine initialize_as_identity(zmatrix)
   complex(pfdp), intent(inout) :: zmatrix(:,:)

   integer :: i, dim, shp(2)

   shp = shape(zmatrix)
   dim = shp(1)

   zmatrix = z0
   forall (i=1:dim) zmatrix(i,i) = z1
 end subroutine initialize_as_identity
end module sweeper
