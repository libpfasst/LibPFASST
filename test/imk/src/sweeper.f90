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
  use mod_zmkpair
  use utils

  implicit none

  external :: zgemm

  type, extends(pf_imk_t) :: imk_sweeper_t
     integer :: dim
     complex(pfdp), allocatable :: commutator(:,:)
   contains
     procedure :: f_eval
     procedure :: dexpinv
     procedure :: propagate => propagate_solution
     procedure :: commutator_p
     procedure :: initialize
     procedure :: destroy
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

  subroutine initialize(this,pf, level_index)
    use probin, only: nparticles, dt,  use_sdc, rk, mkrk,  nterms
    class(imk_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index

!    integer, intent(in) :: level, qtype, nterms
!    logical, intent(in) :: debug, use_sdc, rk, mkrk

 !   class(imk_sweeper_t), pointer :: imk !< context data containing integrals, etc

 !   imk => cast_as_imk_sweeper(this)
    !  Call the imk sweeper initialize


    this%qtype = pf%qtype
    this%nterms = nterms(level_index)
    this%debug = pf%debug
    this%dim = nparticles
    this%use_sdc = use_sdc
    this%rk = rk
    this%mkrk = mkrk

    print *,'calling sweeper initialize',rk,mkrk,use_sdc,nparticles
    call this%imk_initialize(pf,level_index)    

    allocate(this%commutator(nparticles, nparticles))
    this%commutator = z0

!    nullify(imk)
  end subroutine initialize

  subroutine f_eval(this, y, t, level, f)
    use probin, only: toda_periodic

    class(imk_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: y ! prev solution
    class(pf_encap_t), intent(inout) :: f ! output RHS
    real(pfdp), intent(in) :: t
    integer, intent(in) :: level

    type(zmkpair), pointer :: y_p, A_p
    integer :: i

    y_p => cast_as_zmkpair(y)
    A_p => cast_as_zmkpair(f)

    if (nprob .eq. 1) then
       call compute_F_toda(y_p%y,A_p%array,this%dim,t,level)
    else
       call compute_Facke(y_p%y,A_p%array,this%dim,t,level)
    endif


    nullify(y_p, A_p)

  end subroutine f_eval

  subroutine dexpinv(this, a, omega, f)
      class(imk_sweeper_t), intent(inout) :: this
      class(pf_encap_t), intent(inout) :: a, f, omega

      type(zmkpair), pointer :: a_p, omega_p, f_p
      integer :: i
      real(pfdp) :: factor, cc
      complex(pfdp), allocatable :: D(:,:), C(:,:)

      a_p => cast_as_zmkpair(a)
      f_p => cast_as_zmkpair(f)
      omega_p => cast_as_zmkpair(omega)

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

  subroutine commutator_p(this, a, b, out, flags)
    ! this interface routine is just a wrapper to the actual
    ! compute commutator routine and is ONLY used by imk_sweeper_t%rk_step()
    ! for rk_step() it requires the matmul of B and Y
    class(imk_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: a, b, out
    integer, intent(in), optional :: flags

    type(zmkpair), pointer :: a_p, b_p, out_p
    integer :: dim

    a_p => cast_as_zmkpair(a)
    b_p => cast_as_zmkpair(b)
    out_p => cast_as_zmkpair(out)

    dim = a_p%dim

    if (present(flags)) then
       if (flags .eq. 1) b_p%array = b_p%y
       if (flags .eq. 2) then
          b_p%y = b_p%array
          return
       end if
    endif

    call compute_commutator(a_p%array, b_p%array, dim, out_p%array)

  end subroutine commutator_p

  subroutine compute_commutator(a, b, Nmat, output)
    complex(pfdp), intent(in) :: a(Nmat,Nmat), b(Nmat,Nmat)
    integer, intent(in) :: Nmat
    complex(pfdp), intent(inout) :: output(Nmat,Nmat)

    call zgemm('n', 'n', Nmat, Nmat, Nmat, &
         z1, b, Nmat, &
         a, Nmat, &
         z0, output, Nmat) ! output is zeroed here

    call zgemm('n', 'n', Nmat, Nmat, Nmat, &
         z1, a, Nmat, &
         b, Nmat, &
         zm1, output, Nmat)
  end subroutine compute_commutator

  subroutine propagate_solution(this, q0, q)
    class(imk_sweeper_t), intent(inout) :: this
    class(pf_encap_t), intent(inout) :: q0
    class(pf_encap_t), intent(inout) :: q
    integer :: i, dim, nprob=1 !< size of dimensions of P, U
    class(zmkpair), pointer :: q0_p, q_p
    complex(pfdp), allocatable :: tmp(:,:), time_ev_op(:,:)
    real(pfdp) :: exptol=1.0d-20

    q0_p => cast_as_zmkpair(q0)
    q_p => cast_as_zmkpair(q)

    dim = q0_p%dim
    allocate(tmp(dim, dim), time_ev_op(dim, dim))

    if (this%debug) then
      print*, 'expm input'
      call q_p%eprint()
    endif

    tmp = cmplx(0.0, 0.0, pfdp)
    time_ev_op = cmplx(0.0, 0.0, pfdp)
    time_ev_op = compute_matrix_exp(q_p%array, dim, exptol)

    if (this%debug) then
      print*, 'solution before propagate'
      call q0_p%eprint()

      print*, 'time_ev_ops'
      do i = 1, dim
         print*, time_ev_op(i,i)
      end do
    endif

    call zgemm('n', 'n', dim, dim, dim, &
        z1, time_ev_op, dim, &
        q0_p%y, dim, &
        z0, tmp, dim)

    if (nprob < 10) then
       call zgemm('n', 'c', dim, dim, dim, &
            z1, tmp, dim, &
            time_ev_op, dim, &
            z0, q_p%y, dim)
    endif

    if (this%debug) then
      print*, 'solution after propagate'
      print*, 'this omega was used to compute time_ev_op'
      call q_p%eprint()

      print*, 'asymmetry in solution = ', q_p%y(1,1) + q_p%y(dim,dim)
      !if (real(q_p%y(1,1) + q_p%y(dim,dim)) > 0.001) stop
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

     matexp=z0
     norm = compute_inf_norm(matrix_in, dim)
     if (norm .eq. 0.0d0) then
        call initialize_as_identity(matexp,dim)
        return
     else

        ratio = log(2.5_pfdp*norm) / log(2.0_pfdp)
        mscale = max(int(ratio), 0)
        
        ! print*, 'norm=', norm, 'ratio=', ratio, 'mscale=', mscale
        
        scale_val = 1.0_pfdp/(2.0_pfdp**mscale)
        zscale = cmplx(scale_val)
        
        matrix = zscale * matrix_in
        call initialize_as_identity(prev,dim)
        call initialize_as_identity(matexp,dim)
        
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
     endif
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

 !> Destroy sweeper (bypasses base sweeper destroy)
 subroutine destroy(this,pf,level_index)
   class(imk_sweeper_t), intent(inout) :: this
   type(pf_pfasst_t),  target, intent(inout) :: pf
   integer,              intent(in)    :: level_index
   
!   class(imk_sweeper_t), pointer :: sweeper
   !>  Call base sweeper destroy
   print *,'calling base sweeper destroy'
   call this%imk_destroy(pf,level_index)
   
   deallocate(this%commutator)

 end subroutine destroy
 


 subroutine initialize_as_identity_real(matrix,dim)
   integer, intent(in)  :: dim
   real(pfdp), intent(inout) :: matrix(dim,dim)

   integer :: i 

   matrix = 0.0_pfdp
   forall (i=1:dim) matrix(i,i) = 1.0_pfdp
 end subroutine initialize_as_identity_real

 subroutine initialize_as_identity(zmatrix,dim)
   integer, intent(in)  :: dim
   complex(pfdp), intent(inout) :: zmatrix(dim,dim)

   integer :: i 

   zmatrix = z0
   forall (i=1:dim) zmatrix(i,i) = z1
 end subroutine initialize_as_identity
end module sweeper
