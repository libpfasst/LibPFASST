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

module pf_mod_implicitQ
  use pf_mod_implicit
  implicit none
  
  type, extends(pf_implicit_t), abstract :: pf_implicitQ_t
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
     logical                 :: use_LUq_ = .false.
   contains
     procedure :: sweep      => implicitQ_sweep
     procedure :: initialize => implicitQ_initialize
     procedure :: integrate  => implicitQ_integrate
     procedure :: destroy    => implicitQ_destroy
     procedure :: implicitQ_destroy
  end type pf_implicitQ_t

contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine implicitQ_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer

    class(pf_implicitQ_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    class(pf_level_t), intent(inout) :: lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:lev%nnodes-1)
    class(pf_encap_t), allocatable :: rhs

    call start_timer(pf, TLEVEL+lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1
       call lev%S(m)%setval(0.0_pfdp)
       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,1))
       end do
       if (allocated(lev%tauQ)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%unpack(lev%q0)

    call this%f2eval(lev%Q(1), t0, lev%level, lev%F(1,1))

    call lev%ulevel%factory%create0(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       call rhs%setval(0.0_pfdp)
       do n = 1, m
          call rhs%axpy(dt*this%QtilI(m,n), lev%F(n,1))
       end do
       !  Add the tau term
       call rhs%axpy(1.0_pfdp, lev%S(m))
       !  Add the starting value
       call rhs%axpy(1.0_pfdp, lev%Q(1))

       call this%f2comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), rhs, lev%level,lev%F(m+1,1))

    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    call lev%ulevel%factory%destroy0(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call end_timer(pf, TLEVEL+lev%level-1)
  end subroutine implicitQ_sweep

  ! Initialize matrices
  subroutine implicitQ_initialize(this, lev)
    class(pf_implicitQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m,n, nnodes

    this%npieces = 1

    nnodes = lev%nnodes
    allocate(this%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(this%QtilI(nnodes-1,nnodes))  !  S-BE

    this%QtilI = 0.0_pfdp

    if (this%use_LUq_) then 
       ! Get the LU
       call myLUq(lev%qmat,lev%LUmat,lev%nnodes,1)
       this%QtilI = lev%LUmat
    else 
       dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
       do m = 1, nnodes-1
          do n = 1,m
             this%QtilI(m,n+1) =  dsdc(n)
          end do
       end do
    end if

    this%QdiffI = lev%qmat-this%QtilI

  end subroutine implicitQ_initialize

  
  ! Destroy the matrices
  subroutine implicitQ_destroy(this, lev)
    class(pf_implicitQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffI)
    deallocate(this%QtilI)
  end subroutine implicitQ_destroy


  ! Compute SDC integral
  subroutine implicitQ_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_implicitQ_t), intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in   ) :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine implicitQ_integrate

end module pf_mod_implicitQ
