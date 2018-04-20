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
!> Module with useful subroutines that don't  fit in other modules
module pf_mod_utils
  use pf_mod_dtype
  use pf_mod_timer
  implicit none
contains



  !
  !> Compute full residual at each node and measure it's size
  subroutine pf_residual(pf, lev, dt, flag)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: lev
    real(pfdp),        intent(in)    :: dt
    integer, optional, intent(in)    :: flag

    real(pfdp) :: norms(lev%nnodes-1)
    integer :: m, which
    
    which = 1
    if(present(flag)) which = flag

    call start_timer(pf, TRESIDUAL)

    call lev%ulevel%sweeper%residual(lev, dt, which)

    ! compute max residual norm
    do m = 1, lev%nnodes-1
       norms(m) = lev%R(m)%norm(which)
!       print *, 'norm(', m, ') = ', norms(m)
    end do
!    lev%residual = maxval(abs(norms))
    lev%residual = norms(lev%nnodes-1)

    call end_timer(pf, TRESIDUAL)

  end subroutine pf_residual

  !
  !> Generic residual
  !! Each sweeper can define its own residual, or use this generic one
  subroutine pf_generic_residual(this, lev, dt, flags)
    class(pf_sweeper_t), intent(in)  :: this
    class(pf_level_t),  intent(inout) :: lev
    real(pfdp),        intent(in)    :: dt
    integer,  intent(in), optional  :: flags

    integer :: m, which
    
    which = 1
    if(present(flags)) which = flags

    !>  Compute the integral of F
    call lev%ulevel%sweeper%integrate(lev, lev%Q, lev%F, dt, lev%I, which)

    !> add tau if it exists
    if (allocated(lev%tauQ)) then
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), which)
       end do
    end if

    !> subtract out the solution value
    do m = 1, lev%nnodes-1
       if( (which .eq. 0) .or. (which .eq. 1) ) then
         call lev%R(m)%copy(lev%I(m), 1)
         call lev%R(m)%axpy(1.0_pfdp, lev%Q(1), 1)
         call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1), 1)
       end if
       if( (which .eq. 0) .or. (which .eq. 2) ) then
         call lev%R(m)%copy(lev%I(m), 2)
         call lev%R(m)%axpy(1.0_pfdp, lev%Q(lev%nnodes), 2)
         call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m), 2)
       end if
    end do

  end subroutine pf_generic_residual


  !
  !> Generic evaluate all
  !! Each sweeper can define its own evaluate_all or use this generic one
  subroutine pf_generic_evaluate_all(this, lev, t, flags, step)
    class(pf_sweeper_t), intent(in)  :: this
    class(pf_level_t),  intent(inout) :: lev
    real(pfdp),        intent(in)    :: t(:)
    integer, optional, intent(in)    :: flags, step

    integer :: m, which, mystep
        
    which = 1
    if(present(flags)) which = flags
    
    mystep = 1
    if(present(step)) mystep = step
    
    do m = 1, lev%nnodes
       call lev%ulevel%sweeper%evaluate(lev, t(m), m, which, mystep)
    end do
  end subroutine pf_generic_evaluate_all

  
  !> Generic routine to spread initial conditions
  !! Each sweeper can define its own spreadq0 or use this generic one
  subroutine pf_generic_spreadq0(this,lev, t0)
    class(pf_sweeper_t), intent(in)  :: this
    class(pf_level_t), intent(inout) :: lev  !<  Level on which to spread
    real(pfdp),       intent(in)    :: t0    !<  time at beginning of interval

    integer :: m, p

    !  Stick initial condition into first node slot
    call lev%Q(1)%copy(lev%q0)

    !  Evaluate F at first spot
    call lev%ulevel%sweeper%evaluate(lev, t0, 1)

    ! Spread F and solution to all nodes
    do m = 2, lev%nnodes
       call lev%Q(m)%copy(lev%Q(1))
       do p = 1, lev%ulevel%sweeper%npieces
         call lev%F(m,p)%copy(lev%F(1,p))
       end do
    end do
  end subroutine pf_generic_spreadq0

  !>  Routine to compute the LU decomposition of spectral integration matrix
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
    if (maxval(abs(matmul(L,U)-transpose(Q(1:Nnodes-1,2:Nnodes)))) .gt. 1e-14) then
       print *,'LU error in pf_utils'
    endif 
    
    Qtil = 0.0_pfdp
    Qtil(1:Nnodes-1,2:Nnodes)=transpose(U)
    !  Now scale the columns of U to match the sum of A
    if (fillq .eq. 1) then
       do j=1,Nnodes-1
          Qtil(j,1)=sum(Q(j,1:Nnodes))-sum(U(j,1:Nnodes-1))
       end do
    end if

  end subroutine myLUq


end module pf_mod_utils
