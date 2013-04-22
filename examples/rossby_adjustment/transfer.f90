!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use iso_c_binding
  use encap_array1d
  use feval
  implicit none
contains

  subroutine interpolate(qFp, qGp, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

  
  end subroutine interpolate

  subroutine restrict(qFp, qGp, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

  end subroutine restrict
end module transfer
