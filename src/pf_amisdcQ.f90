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

module pf_mod_amisdcQ
  use pf_mod_amisdc
  implicit none

  type, extends(pf_amisdc_t), abstract :: pf_amisdcQ_t
     real(pfdp), allocatable :: QdiffE(:,:)
     real(pfdp), allocatable :: QdiffI(:,:)
     real(pfdp), allocatable :: QtilE(:,:)
     real(pfdp), allocatable :: QtilI(:,:)
   contains 
     procedure :: sweep        => amisdcQ_sweep
     procedure :: initialize   => amisdcQ_initialize
     procedure :: integrate    => amisdcQ_integrate
  end type pf_amisdcQ_t

contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine amisdcQ_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer
    class(pf_amisdcQ_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout) :: pf
    real(pfdp),          intent(in)    :: dt, t0
    class(pf_level_t),   intent(inout) :: lev

    integer                        :: m, n
    real(pfdp)                     :: t
    real(pfdp)                     :: dtsdc(1:lev%nnodes-1)
    class(pf_encap_t), allocatable :: rhsA, rhsB, QA, QB
    class(pf_encap_t), allocatable :: S2(:), S3(:)

    call start_timer(pf, TLEVEL+lev%level-1)
    
    call lev%ulevel%factory%create1(S2,lev%nnodes-1,lev%level,SDC_KIND_SOL_FEVAL,lev%nvars,lev%shape)
    call lev%ulevel%factory%create1(S3,lev%nnodes-1,lev%level,SDC_KIND_SOL_FEVAL,lev%nvars,lev%shape)
    
    
    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1

       call lev%S(m)%setval(0.0_pfdp)
       call S2(m)%setval(0.0d0)
       call S3(m)%setval(0.0d0)

       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
          call lev%S(m)%axpy(dt*lev%qmat(m,n),    lev%F(n,2))
          call lev%S(m)%axpy(dt*lev%qmat(m,n),    lev%F(n,3))
          if (m > n-1) then
             call S2(m)%axpy(dt*this%QtilI(m,n),     lev%F(n,2))
             call S3(m)%axpy(dt*this%QtilI(m,n),     lev%F(n,3))
          else
             call S2(m)%axpy(2.0_pfdp*dt*this%QtilI(m,n),     lev%F(n,2))
             call S3(m)%axpy(2.0_pfdp*dt*this%QtilI(m,n),     lev%F(n,3))
          end if 
       end do
       if (allocated(lev%tauQ)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tauQ(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%unpack(lev%q0)

    call this%f1eval(lev%Q(1), t0, lev%level, lev%F(1,1))
    call this%f2eval(lev%Q(1), t0, lev%level, lev%F(1,2))
    call this%f3eval(lev%Q(1), t0, lev%level, lev%F(1,3))

    call lev%ulevel%factory%create0(rhsA, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%create0(rhsB, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%create0(QA,   lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%create0(QB,   lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)
             
       call rhsA%setval(0.0_pfdp)
       ! First compute the explicit part of the right-hand side
       do n = 1, m
          call rhsA%axpy(dt*this%QtilE(m,n), lev%F(n,1))  
       end do
       call rhsA%axpy(1.0_pfdp, lev%S(m))
       call rhsA%axpy(1.0_pfdp, lev%Q(1))
       
       ! Save the right-hand side with only the explicit contribution
       call rhsB%copy(rhsA)

       ! Add the first implicit part to the right-hand side and solve for the first asynchronous update
       do n = 1, m
          call rhsA%axpy(dt*this%QtilI(m,n), lev%F(n,2))  
       end do
       call rhsA%axpy(-1.0_pfdp, S2(m))  

       call this%f2comp(QA, t, 2.0_pfdp*dt*this%QtilI(m,m+1), rhsA, lev%level, lev%F(m+1,2))

       ! Add the second implicit part to the right-hand side and solve for the second asynchronous update
       do n = 1, m
          call rhsB%axpy(dt*this%QtilI(m,n), lev%F(n,3))  
       end do
       call rhsB%axpy(-1.0_pfdp, S3(m))  

       call this%f3comp(QB, t, 2.0_pfdp*dt*this%QtilI(m,m+1), rhsB, lev%level, lev%F(m+1,3))

       ! Now we average the two asynchronous updates
       call lev%Q(m+1)%setval(0.0_pfdp)
       call lev%Q(m+1)%axpy(0.5_pfdp, QA)
       call lev%Q(m+1)%axpy(0.5_pfdp, QB)

       ! Evaluate the three right-hand sides with the updated variables
       call this%f1eval(lev%Q(m+1), t, lev%level, lev%F(m+1,1))
       call this%f2eval(lev%Q(m+1), t, lev%level, lev%F(m+1,2))
       call this%f3eval(lev%Q(m+1), t, lev%level, lev%F(m+1,3))
    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    call end_timer(pf, TLEVEL+lev%level-1)

  end subroutine amisdcQ_sweep
    
  ! Initialize matrices
  subroutine amisdcQ_initialize(this, lev)
    class(pf_amisdcQ_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m, n, nnodes

    this%npieces = 3

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
  end subroutine amisdcQ_initialize

  ! Compute SDC integral
  subroutine amisdcQ_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_amisdcQ_t),  intent(inout) :: this
    class(pf_level_t),  intent(in)      :: lev
    class(pf_encap_t), intent(in)       :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in)       :: dt
    class(pf_encap_t), intent(inout)    :: fintSDC(:)

    integer :: n, m, p
    
    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,p))
          end do
       end do
    end do  
  end subroutine amisdcQ_integrate
 
end module pf_mod_amisdcQ
