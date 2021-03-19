!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS specification for Dahlquist example.
!>     u_t = lam1*u + lam2*u
module pf_my_stepper
  use pf_mod_ndarray
  use pf_mod_rkstepper

  !>  extend the imex stepper type with stuff we need to compute rhs
  type, extends(pf_ark_stepper_t) :: my_stepper_t

   contains

     procedure :: f_eval    !  Computes the explicit rhs terms
     procedure :: f_comp    !  Does implicit solves 

  end type my_stepper_t

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    use probin, only:  lam1, lam2
    class(my_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    
    real(pfdp),      pointer :: yvec(:), fvec(:)

    !> Grab the arrays from the encapsulation
    yvec  => get_array1d(y)
    fvec => get_array1d(f)

    ! Compute the function values
    select case (piece)
    case (1)  ! Explicit piece
       fvec = lam1*yvec
    case (2)  ! Implicit piece
       fvec = lam2*yvec
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select

  end subroutine f_eval

  !> Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  lam1, lam2
    class(my_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece

    real(pfdp),      pointer :: yvec(:), rhsvec(:), fvec(:)
    
    if (piece == 2) then
       !> Grab the arrays from the encapsulation       
       yvec  => get_array1d(y)
       rhsvec => get_array1d(rhs)
       fvec => get_array1d(f)

       !  Do the solve
       yvec =  rhsvec/(1.0_pfdp - dtq*lam2)

       !  The function is easy to derive  (equivalent to lam2*yvec)
       fvec = (yvec - rhsvec) / dtq
    else
       print *,'Bad piece in f_comp ',piece
       call exit(0)
    end if
  end subroutine f_comp

end module pf_my_stepper

