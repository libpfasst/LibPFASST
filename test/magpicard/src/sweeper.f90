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
  use pf_mod_zndarray
  use utils

  implicit none

  real(pfdp), parameter :: &
       pi = 3.141592653589793_pfdp, &
       two_pi = 6.2831853071795862_pfdp

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
     integer :: indices(3,2)
   contains
     procedure :: f_eval => compute_B
     procedure :: compute_single_commutators
     procedure :: compute_omega
     ! procedure :: compute_time_ev_ops
     procedure :: propagate_solution
     procedure :: destroy
     procedure :: initialize
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

  subroutine initialize_magpicard_sweeper(this,pf, level_index)
    use probin, only: nprob, basis, molecule, magnus_order, nparticles, dt
    class(pf_sweeper_t), intent(inout) :: this
   type(pf_pfasst_t),  target, intent(inout) :: pf
    integer, intent(in) :: level_index


    class(magpicard_sweeper_t), pointer :: magpicard !< context data containing integrals, etc
    integer :: coefs=9

    magpicard => cast_as_magpicard_sweeper(this)

    magpicard%qtype = pf%qtype
    magpicard%dt = dt
    magpicard%magnus_order = magnus_order(level_index)
    magpicard%exp_iterations = 0
    magpicard%debug = pf%debug
    magpicard%dim = nparticles

    ! inouts necessary for the base class structure

    allocate(magpicard%commutator(nparticles, nparticles), &
         magpicard%commutators(magpicard%dim, magpicard%dim, 9))

    magpicard%commutators(:,:,:) = z0
    magpicard%commutator = z0

    magpicard%indices(1,1) = 1
    magpicard%indices(1,2) = 2

    magpicard%indices(2,1) = 1
    magpicard%indices(2,2) = 3

    magpicard%indices(3,1) = 2
    magpicard%indices(3,2) = 3

    nullify(magpicard)
  end subroutine initialize_magpicard_sweeper

    subroutine initialize(this,pf, level_index)
    use probin, only: nprob, basis, molecule, magnus_order, nparticles, dt
    class(magpicard_sweeper_t), intent(inout) :: this
   type(pf_pfasst_t),  target, intent(inout) :: pf
    integer, intent(in) :: level_index


    integer :: ncoefs=9

    this%qtype = pf%qtype
    this%dt = dt
    this%magnus_order = magnus_order(level_index)
    this%exp_iterations = 0
    this%debug = pf%debug
    this%dim = nparticles

    ! inouts necessary for the base class structure
    call this%magpicard_initialize(pf,level_index)

    allocate(this%commutator(nparticles, nparticles), &
         this%commutators(this%dim, this%dim, ncoefs))

    this%commutators(:,:,:) = z0
    this%commutator = z0

    this%indices(1,1) = 1
    this%indices(1,2) = 2

    this%indices(2,1) = 1
    this%indices(2,2) = 3

    this%indices(3,1) = 2
    this%indices(3,2) = 3

  end subroutine initialize


  subroutine compute_B(this, y, t, level, f)
    use probin, only: Nprob

    
    class(magpicard_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: y ! prev solution
    class(pf_encap_t), intent(inout) :: f ! output RHS
    real(pfdp), intent(in) :: t
    integer, intent(in) :: level
    
    type(zndarray), pointer :: L
    complex(pfdp),      pointer :: L_array(:,:), B_array(:,:)
    integer :: i,j,n,m,dhalf,Nmat
    real(pfdp) :: xi,xj,xn,xm,cst

    L => cast_as_zndarray(y)
    L_array=>get_array2d(y)
    B_array=>get_array2d(f)    
    Nmat = L%shape(1)  !  Assume square matrix
    
    if (nprob .eq. 1) then
       call compute_F_toda(L_array,B_array,Nmat,t,level)
    else
       call compute_Facke(L_array,B_array,Nmat,t,level)
    endif
    nullify(L)
    
  end subroutine compute_B
 

  subroutine compute_single_commutators(this, f)
    class(magpicard_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: f(:,:)

    class(zndarray), pointer :: f1, f2
    complex(pfdp),      pointer :: f1_array(:,:), f2_array(:,:)    

    integer :: i, j, k, nnodes, node_offset


    if (this%qtype == 1) then
       node_offset = 0
    else
       node_offset = 1
    endif

    !$omp parallel do private(i, j, k, f1,f2)
    do i = 1, 3
       j = this%indices(i, 1) + node_offset
       k = this%indices(i, 2) + node_offset
        f1 => cast_as_zndarray(f(j,1))
       f2 => cast_as_zndarray(f(k,1))
       f1_array =>get_array2d(f1)
       f2_array =>get_array2d(f2)     
       call compute_commutator(f1_array, f2_array, this%dim, this%commutators(:,:,i))
    enddo
    !$omp end parallel do

    !$omp barrier
    nullify(f1, f2)


  end subroutine compute_single_commutators

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
   integer ::  dim
   complex(pfdp),      pointer :: omega_array(:,:), ints_array(:,:)
   
   dim = this%dim


   omega_array=>get_array2d(omega)
   ints_array=>get_array2d(integrals(this_node))    

   omega_array = ints_array

   if (this%magnus_order > 1) then
      call add_single_commutator_terms(this, omega_array, coefs(:,1), dim)
   endif

   if (this%magnus_order > 2 .and. this%qtype == 5) then
      call add_double_commutator_terms(this, omega_array, f, coefs(:,2), dim)
      call add_triple_commutator_terms(this, omega_array, f, nodes, this_node,&
           qmat, dt, coefs(1,3), dim)
   endif

 end subroutine compute_omega

 subroutine add_single_commutator_terms(this, omega, coefs, N)
   class(magpicard_sweeper_t), intent(inout) :: this
   complex(pfdp), intent(inout) :: omega(:,:)
   integer, intent(in) :: N
   real(pfdp), intent(in) :: coefs(:)

   integer :: i
   complex(pfdp) :: tmp(N,N)

   !$omp parallel do reduction(+:omega) private(i, tmp)
   do i = 1, 3
      tmp = coefs(i) * this%commutators(:,:,i)
      omega = omega + tmp
   end do
   !$omp end parallel do

 end subroutine add_single_commutator_terms

 subroutine add_double_commutator_terms(this, omega, f, coefs, N)
   class(magpicard_sweeper_t), intent(inout) :: this
   class(pf_encap_t), intent(inout) :: f(:,:)
   complex(pfdp), intent(inout) :: omega(:,:)
   real(pfdp), intent(in) :: coefs(:)
   integer, intent(in) :: N

   complex(pfdp),      pointer :: f1_array(:,:)
   integer :: i, j, k, nnodes, node_offset
   complex(pfdp) :: tmp(N,N,3)

   node_offset = 0
   nnodes = size(f)
   if (nnodes > 3) node_offset = 1
   tmp = z0

   do i = 1, 3
      j = i + node_offset
      f1_array=>get_array2d(f(j,1))      
      
      tmp(:,:,1) = tmp(:,:,1) + coefs(i)   * f1_array
      tmp(:,:,2) = tmp(:,:,2) + coefs(i+3) * f1_array
      tmp(:,:,3) = tmp(:,:,3) + coefs(i+6) * f1_array
   end do

   !$omp parallel do reduction(+:omega) private(i, j, k)
   do i = 1, 3
      j = this%indices(i, 1) + node_offset
      k = this%indices(i, 2) + node_offset

      call compute_commutator(tmp(:,:,i), this%commutators(:,:,i), this%dim, this%commutators(:,:,i+3))
      omega = omega + this%commutators(:,:,i+3)
   end do
   !$omp end parallel do

 end subroutine add_double_commutator_terms

 subroutine add_triple_commutator_terms(this, omega, f, nodes, this_node, qmat, dt, coef, N)
   class(magpicard_sweeper_t), intent(inout) :: this
   class(pf_encap_t), intent(inout) :: f(:,:)
   complex(pfdp), intent(inout) :: omega(:,:)
   real(pfdp), intent(in) :: coef, nodes(:), qmat(:,:), dt
   integer, intent(in) :: N, this_node

   class(zndarray), pointer :: f1
   complex(pfdp),      pointer :: f1_array(:,:)
   
   complex(pfdp), allocatable :: tmp(:,:), a(:,:,:)
   real(pfdp) :: time_scaler
   integer :: i, j, nnodes

   nnodes = size(nodes)
   this%commutator = z0

   allocate(tmp(N,N), a(N,N,3))
   tmp = z0
   a = z0

   do i = 2, nnodes-1
      f1_array =>get_array2d(f(i,1))
       time_scaler = (nodes(i)-0.5_pfdp)
       tmp = dt * qmat(this_node, i) * f1_array
       a(:,:,1) = a(:,:,1) + tmp
       a(:,:,2) = a(:,:,2) + tmp * time_scaler
       a(:,:,3) = a(:,:,3) + tmp * time_scaler**2
   enddo

   tmp = a(:,:,2)
   call compute_commutator(a(:,:,1), tmp, this%dim, this%commutator)
   a(:,:,3) = this%commutator
   call compute_commutator(a(:,:,1), a(:,:,3), this%dim, this%commutator)
   a(:,:,3) = this%commutator
   call compute_commutator(a(:,:,1), a(:,:,3), this%dim, this%commutator)

   omega = omega + this%commutator / 60.d0

   deallocate(tmp, a)
 end subroutine add_triple_commutator_terms

 subroutine compute_commutator(a, b, N, output)
   complex(pfdp), intent(in) :: a(:,:), b(:,:)
   integer, intent(in) :: N
   complex(pfdp), intent(inout) :: output(:,:)

   call zgemm('n', 'n', N, N, N, &
        z1, b, N, &
        a, N, &
        z0, output, N) ! output is zeroed here

   call zgemm('n', 'n', N, N, N, &
        z1, a, N, &
        b, N, &
        zm1, output, N)
 end subroutine compute_commutator


 !> Computes the P_t = U*P_t0*U^dagger
 subroutine propagate_solution(this, sol_t0, sol_tn, omega, level)
   use probin, only: nprob, exptol
   class(magpicard_sweeper_t), intent(inout) :: this
   class(pf_encap_t), intent(inout) :: sol_t0
   class(pf_encap_t), intent(inout) :: sol_tn
   class(pf_encap_t), intent(inout) :: omega !< Time-evolution operator
   integer, intent(in) :: level
   integer :: dim !< size of dimensions of P, U
   class(zndarray), pointer :: sol_t0_p ! , sol_tn_p, omega_p
   complex(pfdp), allocatable :: tmp(:,:), time_ev_op(:,:)
   complex(pfdp),      pointer :: sol_t0_array(:,:), sol_tn_array(:,:),omega_array(:,:)
   
   sol_t0_p =>cast_as_zndarray(sol_t0)
   sol_t0_array =>get_array2d(sol_t0)
   sol_tn_array =>get_array2d(sol_tn)       
   omega_array =>get_array2d(omega)       

   dim = sol_t0_p%shape(1)  !  Assumes square matrix
   allocate(tmp(dim, dim), time_ev_op(dim, dim))

   time_ev_op = cmplx(0.0, 0.0, pfdp)
   time_ev_op = compute_matrix_exp(omega_array, dim, exptol(level))

   if (nprob < 10) then
      call zgemm('n', 'n', dim, dim, dim, &
           z1, time_ev_op, dim, &
           sol_t0_array, dim, &
           z0, tmp, dim)

      call zgemm('n', 'c', dim, dim, dim, &
           z1, tmp, dim, &
           time_ev_op, dim, &
           z0, sol_tn_array, dim)
   else
      call zgemm('n', 'n', dim, dim, dim, &
           z1, time_ev_op, dim, &
           sol_t0_array, dim, &
           z0, sol_tn_array, dim)
   endif

   deallocate(tmp, time_ev_op)

 end subroutine propagate_solution

function compute_matrix_exp(matrix_in, dim, tol) result(matexp)
   ! sum and square method
   integer, intent(in) :: dim
   complex(pfdp), intent(in) :: matrix_in(dim, dim)
   real(pfdp), intent(in) :: tol
   complex(pfdp) :: matexp(dim, dim)

   integer, parameter :: MAX_TERMS = 3
   integer, parameter :: maxk = 1000
   integer, parameter :: max_mscale = 16

   complex(pfdp) :: zinvk, zscale, matrix(dim,dim), prev(dim, dim), next(dim, dim)
   real(pfdp) :: invk, ratio, norm, scale_val
   integer :: i, ik, im, nterms, mscale
   logical :: converged

   matexp = z0
   norm = compute_inf_norm(matrix_in, dim)
   ratio = log(2.5_pfdp*norm) / log(2.0_pfdp)
   mscale = max(int(ratio), 0)



   scale_val = 1.0_pfdp/(2.0_pfdp**mscale)
   zscale = cmplx(scale_val)

   matrix = zscale * matrix_in
   call initialize_as_identity(prev)
   call initialize_as_identity(matexp)

   ik = 1
   nterms = 0
   next = z0
   converged = .false.

   do while (.not. converged)
      zinvk = z1 / cmplx(ik)
      call zgemm('n', 'n', dim, dim, dim, &
           zinvk, prev, dim, &
           matrix, dim, &
           z0, next, dim)
      matexp = matexp + next

      norm = compute_inf_norm(next, dim)
      if(norm < tol) nterms = nterms + 1
      if(nterms >= MAX_TERMS) converged = .true.

      prev = next

      ik = ik + 1
   end do

   next = matexp

   do im = 1, mscale
      prev = z0
      if (maxval(abs(matexp))*0.0 /= 0.0) then
         do ik = 1, 92
            next = matmul(next, next)
            if (maxval(abs(next))*0.0 /= 0.0) stop
         end do
      endif

      call zgemm('n', 'n', dim, dim, dim, &
           z1, matexp, dim, &
           matexp, dim, &
           z0, prev, dim)

      matexp = prev
   enddo
 end function compute_matrix_exp

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


 !> Destroy sweeper (called from base sweeper destroy from pf_pfasst_destroy)
  subroutine destroy(this,pf,level_index)
   class(magpicard_sweeper_t), intent(inout) :: this
   type(pf_pfasst_t),  target, intent(inout) :: pf
   integer,              intent(in)    :: level_index

   !>  Destroy the magpicard sweeper
   call this%magpicard_destroy(pf,level_index)

   !>  Destroy local variables
   deallocate(this%commutators)

 end subroutine destroy

 
  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(magpicard_context), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
   integer,           intent(in   ), optional :: flags

   class(zndarray), pointer :: f, g
   f => cast_as_zndarray(f_vec)
   g => cast_as_zndarray(c_vec)

   f%flatarray = g%flatarray
 end subroutine interpolate


 subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
   class(magpicard_context), intent(inout) :: this
   class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
   class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
   real(pfdp),        intent(in   ) :: t
   integer,           intent(in   ), optional :: flags

   class(zndarray), pointer :: f, g
   f => cast_as_zndarray(f_vec)
   g => cast_as_zndarray(c_vec)

   g%flatarray = f%flatarray
 end subroutine restrict

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
