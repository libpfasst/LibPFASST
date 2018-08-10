!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 1-D advection/diffusion example.
!>     u_t + v*u_x = nu*u_xx
module feval
  use pf_mod_dtype
  use encap
  use pf_mod_imexQ


  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: my_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type my_level_t

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imexQ_t) :: my_sweeper_t

   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves 

  end type my_sweeper_t

contains

  !>  Helper function to return sweeper pointer
  function as_my_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(my_sweeper_t), pointer :: r
    select type(sweeper)
    type is (my_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_my_sweeper

  !>  Routine to set up sweeper variables and operators
  subroutine sweeper_setup(sweeper, grid_shape)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: grid_shape(1)

    class(my_sweeper_t), pointer :: this
    this => as_my_sweeper(sweeper)

    !>  Set variables for explicit and implicit parts
    this%implicit=.TRUE.
    this%explicit=.TRUE.

  end subroutine sweeper_setup

  !>  destroy the sweeper type
  subroutine destroy(this, lev)
    class(my_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev

    call this%imexQ_destroy(lev)

  end subroutine destroy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    use probin, only:  lam1, lam2
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    class(scalar_encap), pointer :: y_encap, f_encap
    y_encap => cast_as_scalar(y)
    f_encap => cast_as_scalar(f)

    ! Compute the function values
    select case (piece)
    case (1)  ! Explicit piece
       f_encap%y = lam1*y_encap%y
    case (2)  ! Implicit piece
       f_encap%y = lam2*y_encap%y
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select

  end subroutine f_eval

  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  lam1, lam2
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    class(scalar_encap), pointer :: y_encap, f_encap, rhs_encap
    y_encap => cast_as_scalar(y)
    f_encap => cast_as_scalar(f)
    rhs_encap => cast_as_scalar(rhs)

    
    if (piece == 2) then

       !  Do the solve
       y_encap%y =  rhs_encap%y/(1.0_pfdp - dtq*lam2)

       !  The function is easy to derive  (equivalent to lam2*y)
       f_encap%y = (y_encap%y - rhs_encap%y) / dtq
    else
       print *,'Bad piece in f_comp ',piece
       call exit(0)
    end if
  end subroutine f_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  These are the transfer functions that must be  provided for the level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF,qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags


    class(scalar_encap), pointer :: y_encap_f, y_encap_c

    y_encap_f => cast_as_scalar(qF)
    y_encap_c => cast_as_scalar(qG)

    y_encap_f%y=y_encap_c%y
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, levelf, levelG, qF, qG, t, flags)
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelf  !<  fine level
    class(pf_level_t), intent(inout) :: levelG  !<  coarse level
    class(pf_encap_t), intent(inout) :: qF    !<  fine solution
    class(pf_encap_t), intent(inout) :: qG    !<  coarse solution
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags

    class(scalar_encap), pointer :: y_encap_f, y_encap_c

    y_encap_f => cast_as_scalar(qF)
    y_encap_c => cast_as_scalar(qG)

    y_encap_c%y = y_encap_f%y
  end subroutine restrict


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial(y_0)
    type(scalar_encap), intent(inout) :: y_0
    call exact(0.0_pfdp, y_0%y)
  end subroutine initial

  !> Routine to return the exact solution
  subroutine exact(t, yex)
    use probin, only: lam1,lam2
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex

    yex=exp((lam1+lam2)*t)

  end subroutine exact


end module feval
