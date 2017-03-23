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

module pf_mod_imexQ
  use pf_mod_imex
  implicit none

  type, extends(pf_imex_t), abstract :: pf_imexQ_t
     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
   contains
     procedure :: sweep      => imexQ_sweep
     procedure :: initialize => imexQ_initialize
     procedure :: integrate  => imexQ_integrate
     procedure :: destroy    => imexQ_destroy
     procedure :: imexQ_destroy
  end type pf_imexQ_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine imexQ_sweep(this, pf, Lev, t0, dt)
    use pf_mod_timer

    class(pf_imexQ_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    class(pf_level_t), intent(inout) :: lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:Lev%nnodes-1)
    class(pf_encap_t), allocatable :: rhs

    call start_timer(pf, TLEVEL+lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1
       call lev%S(m)%setval(0.0_pfdp)
       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
          call lev%S(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2))
       end do
       if (allocated(lev%tauQ)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%unpack(lev%q0)

    call this%f1eval(lev%Q(1), t0, lev%level, lev%F(1,1))
    call this%f2eval(lev%Q(1), t0, lev%level, lev%F(1,2))

    call lev%ulevel%factory%create0(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       call rhs%setval(0.0_pfdp)
       do n = 1, m
          call rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1))
          call rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2))
       end do
       !  Add the tau term
       call rhs%axpy(1.0_pfdp, lev%S(m))
       !  Add the starting value
       call rhs%axpy(1.0_pfdp, lev%Q(1))

       call this%f2comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), rhs, lev%level,lev%F(m+1,2))
       call this%f1eval(lev%Q(m+1), t, lev%level, lev%F(m+1,1))

    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    call lev%ulevel%factory%destroy0(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    call end_timer(pf, TLEVEL+lev%level-1)
  end subroutine imexQ_sweep

  ! Initialize matrices
  subroutine imexQ_initialize(this, lev)
    class(pf_imexQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m,n, nnodes

    this%npieces = 2

    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(this%QtilE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QtilI(nnodes-1,nnodes))  !  S-BE

    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       do n = 1,m
          this%QtilE(m,n)   =  dsdc(n)
          this%QtilI(m,n+1) =  dsdc(n)
       end do
    end do

    this%QdiffE = lev%qmat-this%QtilE
    this%QdiffI = lev%qmat-this%QtilI

  end subroutine imexQ_initialize

  subroutine imexQ_destroy(this, lev)
    class(pf_imexQ_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
  end subroutine imexQ_destroy


  ! Compute SDC integral
  subroutine imexQ_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_imexQ_t), intent(inout) :: this
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
  end subroutine imexQ_integrate

end module pf_mod_imexQ
