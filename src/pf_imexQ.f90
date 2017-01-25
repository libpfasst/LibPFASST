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
  end type pf_imexQ_t

contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine imexQ_sweep(this, pf, Lev, t0, dt)
    use pf_mod_timer

    class(pf_imexQ_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    class(pf_level_t), intent(inout) :: Lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:Lev%nnodes-1)
    class(pf_encap_t), allocatable :: rhs

    call start_timer(pf, TLEVEL+Lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, Lev%nnodes-1
       call Lev%S(m)%setval(0.0_pfdp)
       do n = 1, Lev%nnodes
          call Lev%S(m)%axpy(dt*this%QdiffE(m,n), Lev%F(n,1))
          call Lev%S(m)%axpy(dt*this%QdiffI(m,n), Lev%F(n,2))
       end do
       if (allocated(Lev%tauQ)) then
          call Lev%S(m)%axpy(1.0_pfdp, Lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call Lev%Q(1)%unpack(Lev%q0)

    call this%f1eval(Lev%Q(1), t0, Lev%level, Lev%F(1,1))
    call this%f2eval(Lev%Q(1), t0, Lev%level, Lev%F(1,2))

    call Lev%ulevel%factory%create0(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape)

    t = t0
    dtsdc = dt * (Lev%nodes(2:Lev%nnodes) - Lev%nodes(1:Lev%nnodes-1))
    do m = 1, Lev%nnodes-1
       t = t + dtsdc(m)

       call rhs%setval(0.0_pfdp)
       do n = 1, m
          call rhs%axpy(dt*this%QtilE(m,n), Lev%F(n,1))
          call rhs%axpy(dt*this%QtilI(m,n), Lev%F(n,2))
       end do


!       call rhs%axpy(dtsdc(m), Lev%F(m,1))
       call rhs%axpy(1.0_pfdp, Lev%S(m))
       !  Add the starting value
       call rhs%axpy(1.0_pfdp, Lev%Q(1))


!       call this%f2comp(Lev%Q(m+1), t, dtsdc(m), rhs, Lev%level, Lev%levelctx, Lev%F(m+1,2))
       call this%f2comp(Lev%Q(m+1), t, dt*this%QtilI(m,m+1), rhs, Lev%level,Lev%F(m+1,2))
       call this%f1eval(Lev%Q(m+1), t, Lev%level, Lev%F(m+1,1))

    end do

    call Lev%qend%copy(Lev%Q(Lev%nnodes))

    ! done
    !call Lev%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+Lev%level-1)
  end subroutine imexQ_sweep

  ! Initialize matrices
  subroutine imexQ_initialize(this, Lev)
    class(pf_imexQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: Lev

    real(pfdp) :: dsdc(Lev%nnodes-1)
    integer    :: m,n, nnodes

    this%npieces = 2

    nnodes = Lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QdiffI(nnodes-1,nnodes))  !  S-BE
    allocate(this%QtilE(nnodes-1,nnodes))  !  S-FE
    allocate(this%QtilI(nnodes-1,nnodes))  !  S-BE

    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    dsdc = Lev%nodes(2:nnodes) - Lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       do n = 1,m
          this%QtilE(m,n)   =  dsdc(n)
          this%QtilI(m,n+1) =  dsdc(n)
       end do
    end do

    ! do m = 1,nnodes-1
    !    print *,'row i of qmat', m,Lev%qmat(m,:)
    ! end do
    call myLUq(Lev%qmat,this%QtilI,Nnodes,0)
    this%QdiffE = Lev%qmat-this%QtilE
    this%QdiffI = Lev%qmat-this%QtilI

  end subroutine imexQ_initialize

  ! Compute SDC integral
  subroutine imexQ_integrate(this, Lev, qSDC, fSDC, dt, fintSDC)
    class(pf_imexQ_t), intent(inout) :: this
    class(pf_level_t), intent(in   ) :: Lev
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in   ) :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, Lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, Lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*Lev%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine imexQ_integrate

end module pf_mod_imexQ
