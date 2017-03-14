!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_implicit
  implicit none

  real(pfdp), parameter :: &
       a      = 1.0_pfdp,  &  
       b      = 2.0_pfdp,  &  
       c      = 3.0_pfdp    

  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_implicit_t) :: ad_sweeper_t
   contains
     procedure :: f2eval
     procedure :: f2comp
!     final :: destroy0, destroy1
  end type ad_sweeper_t

contains

  function as_ad_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(ad_sweeper_t), pointer :: r
    select type(sweeper)
    type is (ad_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_ad_sweeper

  subroutine setup(sweeper, nvars)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: nvars

    class(ad_sweeper_t), pointer :: this
    integer     :: i

    this => as_ad_sweeper(sweeper)

  end subroutine setup

  subroutine destroy0(this)
    type(ad_sweeper_t), intent(inout) :: this

  end subroutine destroy0

  subroutine destroy1(this)
    type(ad_sweeper_t), intent(inout) :: this(:)
    integer :: i
    do i = 1, size(this)
       call destroy0(this(i))
    end do
  end subroutine destroy1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0
    call exact(0.0_pfdp, q0%flatarray)
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, yex)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    yex(1) = dexp((a+b+c)*t)
  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine f2eval(this, y, t, level, f2)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f2
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level

    real(pfdp), pointer :: yvec(:), f2vec(:)
    real(pfdp)          :: d             
    real(pfdp)          :: val

    d = (a+b+c)

    yvec  => array1(y)
    f2vec => array1(f2)

    f2vec(1) = d*yvec(1)
    
  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine f2comp(this, y, t, dt, rhs, level, f2)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: rhs
    class(pf_encap_t),   intent(inout) :: y, f2
    real(pfdp),          intent(in   ) :: t, dt
    integer,             intent(in   ) :: level

    real(pfdp), pointer :: yvec(:), rhsvec(:), f2vec(:)
    real(pfdp)          :: val
    real(pfdp)          :: d             

    d = (a+b+c)

    yvec   => array1(y)
    rhsvec => array1(rhs)
    f2vec  => array1(f2)
 
    yvec  = rhsvec / (1.0_pfdp - d*dt)
    f2vec = (yvec - rhsvec) / dt 

  end subroutine f2comp

  subroutine interpolate(this, levelF, levelG, qFp, qGp, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qFp, qGp
    real(pfdp),        intent(in   ) :: t

    ! do nothing

  end subroutine interpolate

  subroutine restrict(this, levelF, levelG, qFp, qGp, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qFp, qGp
    real(pfdp),        intent(in   ) :: t
    
    ! do nothing

  end subroutine restrict

end module feval
