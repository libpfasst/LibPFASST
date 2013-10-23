!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use iso_c_binding
  use feval
  implicit none
contains

  subroutine interpolate(qFp, qGp, levelF, levelctxF, levelG, levelctxG,t)
    type(c_ptr), intent(in), value :: qFp, qGp, levelctxF, levelctxG
    integer,     intent(in)        :: levelF, levelG
    real(pfdp),  intent(in) :: t

    type(ad_work_t), pointer :: workF, workG
    real(pfdp),      pointer :: qF(:), qG(:)
    complex(kind=8), pointer :: wkF(:), wkG(:)

    integer :: nvarF, nvarG, xrat

    call c_f_pointer(levelctxF, workF)
    call c_f_pointer(levelctxG, workG)

    qF => array1(qFp)
    qG => array1(qGp)

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       qF = qG
       return
    endif

    wkF => workF%wk
    wkG => workG%wk

    wkG = qG
    call fftw_execute_dft(workG%ffft, wkG, wkG)
    wkG = wkG / nvarG

    wkF = 0.0d0
    wkF(1:nvarG/2) = wkG(1:nvarG/2)
    wkF(nvarF-nvarG/2+2:nvarF) = wkG(nvarG/2+2:nvarG)

    call fftw_execute_dft(workF%ifft, wkF, wkF)

    qF = real(wkF)
  end subroutine interpolate

  subroutine restrict(qFp, qGp, levelF, levelctxF, levelG, levelctxG,t)
    type(c_ptr), intent(in), value :: qFp, qGp, levelctxF, levelctxG
    integer,     intent(in)        :: levelF, levelG
    real(pfdp),  intent(in) :: t

    real(pfdp), pointer :: qF(:), qG(:)

    integer :: nvarF, nvarG, xrat

    qF => array1(qFp)
    qG => array1(qGp)

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    qG = qF(::xrat)
  end subroutine restrict
end module transfer
