!
! This file is part of LIBPFASST.
!
!> Sweeper and RHS specification for Dahlquist example.
!>     u_t = lam1*u + lam2*u
module pf_my_sweeper
  use encap
  use pf_mod_imex_sweeper

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: my_sweeper_t

   contains

     procedure :: f_eval    !  Computes the explicit rhs terms
     procedure :: f_comp    !  Does implicit solves 

     procedure :: initialize !  Overwrites imex sweeper initialize
     procedure :: destroy    !  Overwrites imex sweeper destroy

  end type my_sweeper_t

contains


  !>  Routine to set up sweeper variables and operators
  subroutine initialize(this, pf,level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    integer, intent(in) :: level_index

    !>  Call the imex sweeper initialization
    call this%imex_initialize(pf,level_index)

    !>  Set variables for explicit and implicit parts (just to show you can)
    this%implicit=.TRUE.
    this%explicit=.TRUE.

  end subroutine initialize

  !>  destroy the sweeper type
  subroutine destroy(this, pf,level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    integer, intent(in) :: level_index

    !>  Call the imex sweeper destroy
    call this%imex_destroy(pf,level_index)

    !  Nothing to do 

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines to help out
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


end module pf_my_sweeper

