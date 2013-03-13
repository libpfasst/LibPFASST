!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use iso_c_binding
  use encap
  implicit none
contains

  subroutine interpolate(qF, qG, levelF, ctxF, levelG, ctxG)
    use feval
    type(pf_encap_t), intent(inout) :: qF
    type(pf_encap_t), intent(in)    :: qG
    integer,          intent(in)    :: levelF, levelG
    type(c_ptr),      intent(in)    :: ctxF, ctxG

    complex(kind=8), pointer :: wkF(:), wkG(:)
    integer :: nvarF, nvarG, xrat
    
    nvarF = size(qF%array) 
    nvarG = size(qG%array)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       qF%array = qG%array
       return
    endif

    wkF => levels(levelF)%wk
    wkG => levels(levelG)%wk

    wkG = qG%array
    call fftw_execute_dft(levels(levelG)%ffft, wkG, wkG)
    wkG = wkG / nvarG

    wkF = 0.0d0
    wkF(1:nvarG/2) = wkG(1:nvarG/2)
    wkF(nvarF-nvarG/2+2:nvarF) = wkG(nvarG/2+2:nvarG)

    call fftw_execute_dft(levels(levelF)%ifft, wkF, wkF)

    qF%array = real(wkF)
  end subroutine interpolate

  subroutine restrict(qF, qG, levelF, ctxF, levelG, ctxG)
    integer,          intent(in)    :: levelF, levelG
    type(pf_encap_t), intent(in)    :: qF
    type(pf_encap_t), intent(inout) :: qG
    type(c_ptr),      intent(in)    :: ctxF, ctxG

    integer :: nvarF, nvarG, xrat
    
    nvarF = size(qF%array) 
    nvarG = size(qG%array)
    xrat  = nvarF / nvarG

    qG%array = qF%array(::xrat)
  end subroutine restrict
end module transfer
