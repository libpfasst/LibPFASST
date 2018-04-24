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
  use pf_mod_imk
  use factory
  use utils

  implicit none

  real(pfdp), parameter :: &
       pi = 3.141592653589793_pfdp, &
       two_pi = 6.2831853071795862_pfdp


  external :: zgemm

  type, extends(pf_user_level_t) :: imk_context
   contains
     procedure :: restrict => restrict
     procedure :: interpolate => interpolate
  end type imk_context

  type, extends(pf_imk_t) :: imk_sweeper_t
     integer :: dim
     complex(pfdp), allocatable :: commutator(:,:)
   contains
     procedure :: f_eval
     procedure :: dexpinv
     procedure :: propagate => propagate_solution
     procedure :: destroy => destroy_imk_sweeper
  end type imk_sweeper_t

contains

  function cast_as_imk_sweeper(sweeper) result(imk)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(imk_sweeper_t), pointer :: imk

    select type(sweeper)
    type is (imk_sweeper_t)
       imk => sweeper
    class default
       print*, 'invalid sweeper class'
       stop
    end select

  end function cast_as_imk_sweeper

  subroutine initialize_imk_sweeper(this, level, debug, use_sdc, qtype, nterms)
    use probin, only: nparticles, dt
    class(pf_sweeper_t), intent(inout) :: this
    integer, intent(in) :: level, qtype, nterms
    logical, intent(in) :: debug, use_sdc

    class(imk_sweeper_t), pointer :: imk !< context data containing integrals, etc

    imk => cast_as_imk_sweeper(this)

    imk%qtype = qtype
    imk%nterms = nterms
    imk%debug = debug
    imk%dim = nparticles
    imk%use_sdc = use_sdc

    allocate(imk%commutator(nparticles, nparticles))

    imk%commutator = z0

    nullify(imk)
  end subroutine initialize_imk_sweeper

  subroutine f_eval(this, y, t, level, f)
    use probin, only: toda_periodic

    class(imk_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: y ! prev solution
    class(pf_encap_t), intent(inout) :: f ! output RHS
    real(pfdp), intent(in) :: t
    integer, intent(in) :: level

    type(zndarray), pointer :: y_p, A_p
    integer :: i

    y_p => cast_as_zndarray(y)
    A_p => cast_as_zndarray(f)

    do i = 1, this%dim
       A_p%array(i,i) = 0.0_pfdp
    enddo

    do i = 1, this%dim-1
       A_p%array(i, i+1) = -1.0_pfdp * y_p%y(i, i+1)
       A_p%array(i+1, i) = y_p%y(i, i+1)
    enddo

    if (toda_periodic .eqv. .true.) then
       A_p%array(1, this%dim) = y_p%y(1, this%dim)
       A_p%array(this%dim, 1) = -1.0_pfdp * y_p%y(this%dim, 1)
    endif

    nullify(y_p, A_p)

  end subroutine f_eval

  subroutine dexpinv(this, a, omega, f)
      class(imk_sweeper_t), intent(inout) :: this
      class(pf_encap_t), intent(inout) :: a, f, omega

      type(zndarray), pointer :: a_p, omega_p, f_p
      integer :: i
      real(pfdp) :: factor, cc
      complex(pfdp), allocatable :: D(:,:), C(:,:)

      a_p => cast_as_zndarray(a)
      f_p => cast_as_zndarray(f)
      omega_p => cast_as_zndarray(omega)

      allocate(D(this%dim, this%dim), C(this%dim, this%dim))

      D = a_p%array
      C = a_p%array
      factor = 1.0_pfdp
      do i = 1, this%nterms
         call compute_commutator(omega_p%array, C, this%dim, this%commutator)
         factor = factor / real(i, pfdp)

         if (this%bernoullis(i) .ne. 0.0) then
            cc = this%bernoullis(i) * factor
            D = D + cc * this%commutator
         endif

         C = this%commutator
      end do

      f_p%array = D
      deallocate(D, C)

      nullify(a_p, f_p, omega_p)

  end subroutine dexpinv

  subroutine compute_commutator(a, b, dim, output)
    complex(pfdp), intent(in) :: a(dim,dim), b(dim,dim)
    integer, intent(in) :: dim
    complex(pfdp), intent(inout) :: output(dim,dim)

    call zgemm('n', 'n', dim, dim, dim, &
         z1, b, dim, &
         a, dim, &
         z0, output, dim) ! output is zeroed here

    call zgemm('n', 'n', dim, dim, dim, &
         z1, a, dim, &
         b, dim, &
         zm1, output, dim)
  end subroutine compute_commutator

  subroutine propagate_solution(this, q0, q)
    class(imk_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: q0
    class(pf_encap_t), intent(inout) :: q
    integer :: dim, nprob=1 !< size of dimensions of P, U
    class(zndarray), pointer :: q0_p, q_p
    complex(pfdp), allocatable :: tmp(:,:), time_ev_op(:,:)
    real(pfdp) :: exptol=1.d-15

    q0_p => cast_as_zndarray(q0)
    q_p => cast_as_zndarray(q)

    dim = q0_p%dim
    allocate(tmp(dim, dim), time_ev_op(dim, dim))

    time_ev_op = cmplx(0.0, 0.0, pfdp)
    time_ev_op = compute_matrix_exp(q_p%array, dim, exptol)

    if (nprob < 10) then
       call zgemm('n', 'n', dim, dim, dim, &
            z1, time_ev_op, dim, &
            q0_p%y, dim, &
            z0, tmp, dim)

       call zgemm('n', 'c', dim, dim, dim, &
            z1, tmp, dim, &
            time_ev_op, dim, &
            z0, q_p%y, dim)
    else
       call zgemm('n', 'n', dim, dim, dim, &
            z1, time_ev_op, dim, &
            q0_p%y, dim, &
            z0, q_p%y, dim)
    endif

    deallocate(tmp, time_ev_op)
    nullify(q0_p, q_p)
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
     norm = compute_inf_norm(matrix, dim)
     ratio = log(2.5_pfdp*norm) / log(2.0_pfdp)
     mscale = max(int(ratio), 0)

     ! print*, 'norm=', norm, 'ratio=', ratio, 'mscale=', mscale

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
              print*, 'ik=', ik
              print*, next
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

 !> array of ctx data deallocation
 subroutine destroy_imk_sweeper(this, lev)
   class(imk_sweeper_t), intent(inout) :: this
   class(pf_level_t),   intent(inout) :: lev
   integer :: io

   deallocate(this%commutator)
   call this%imk_destroy(lev)

 end subroutine destroy_imk_sweeper

 subroutine restrict(this, levelF, levelG, qF, qG, t,flags)
   class(imk_context), intent(inout) :: this
   class(pf_level_t), intent(inout) :: levelF
   class(pf_level_t), intent(inout) :: levelG
   class(pf_encap_t), intent(inout) :: qF, qG
   real(pfdp),        intent(in   ) :: t
   integer, intent(in), optional :: flags
   class(zndarray), pointer :: f, g
   f => cast_as_zndarray(qF)
   g => cast_as_zndarray(qG)

   g%array = f%array
   g%y = f%y
 end subroutine restrict

 subroutine interpolate(this, levelF, levelG, qF, qG, t,flags)
   class(imk_context), intent(inout) :: this
   class(pf_level_t), intent(inout) :: levelF
   class(pf_level_t), intent(inout) :: levelG
   class(pf_encap_t), intent(inout) :: qF, qG
   real(pfdp),        intent(in   ) :: t
   integer, intent(in), optional :: flags

   class(zndarray), pointer :: f, g
   f => cast_as_zndarray(qF)
   g => cast_as_zndarray(qG)

   f%array = g%array
   f%y = g%y
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
