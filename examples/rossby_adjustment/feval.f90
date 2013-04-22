!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use iso_c_binding
  use encap_array1d
  use spatialdiscretization, only : RHS, RHS_coarse, ReadInitialvalue

  implicit none

  type :: weno
     double precision, allocatable :: Y(:)
  end type weno
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine feval_create_workspace(ctx, nvars)
    type(c_ptr), intent(out) :: ctx
    integer,     intent(in)  :: nvars

    type(weno), pointer :: work

    allocate(work)
    allocate(work%Y(nvars))
    work%Y = DBLE(1.0)
    ctx = c_loc(work)

  end subroutine feval_create_workspace

  subroutine feval_destroy_workspace(ctx)
    type(c_ptr), intent(in) :: ctx

  end subroutine feval_destroy_workspace

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    type(array1d), intent(inout) :: q0

    CALL ReadInitialValue(q0%array)
    
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, nvars, yex)
    real(pfdp), intent(in)  :: t
    integer,    intent(in)  :: nvars
    real(pfdp), intent(out) :: yex(nvars)

    yex = DBLE(2.0)

  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine eval_f1(yptr, t, level, ctx, f1ptr)
    type(c_ptr), intent(in), value :: yptr, f1ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level

    type(weno), pointer :: y, fy, c
    real(pfdp), pointer :: u(:), fu(:)
    
! yptr points  to y, y points to weno datatype
! f1ptr points to fy, fy points to weno datatype

    call c_f_pointer(yptr, y)
    call c_f_pointer(f1ptr, fy)
    call c_f_pointer(ctx, c)
    
    u => array(yptr)
    fu => array(f1ptr)

    CALL RHS(u, fu, t, 'F', 'F')

  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(yptr, t, level, ctx, f2ptr)
    type(c_ptr), intent(in), value :: yptr, f2ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level

  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine comp_f2(yptr, t, dt, rhsptr, level, ctx, f2ptr)
    type(c_ptr), intent(in), value :: yptr, rhsptr, f2ptr, ctx
    real(pfdp),  intent(in)        :: t, dt
    integer,     intent(in)        :: level

  end subroutine comp_f2

end module feval
