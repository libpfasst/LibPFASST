!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use feval
  use probin
  implicit none
contains

  subroutine interp1(qF, qG, levelctxF, levelctxG)
    type(c_ptr), intent(in), value :: levelctxF, levelctxG
    real(pfdp),  intent(inout) :: qF(:), qG(:)

    type(work1),   pointer :: workF, workG
    complex(pfdp), pointer :: wkF(:), wkG(:)
    integer                :: nvarF, nvarG, xrat

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       qF = qG
       return
    endif

    call c_f_pointer(levelctxF, workF)
    call c_f_pointer(levelctxG, workG)

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
  end subroutine interp1

  subroutine interp2(qF, qG, levelctxF, levelctxG)
    type(c_ptr), intent(in), value :: levelctxF, levelctxG
    real(pfdp),  intent(inout) :: qF(:,:), qG(:,:)

    type(work2),   pointer :: workF, workG
    complex(pfdp), pointer :: wkF(:,:), wkG(:,:)
    integer                :: nvarF, nvarG, xrat, i, j

    nvarF = size(qF, dim=1)
    nvarG = size(qG, dim=1)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       qF = qG
       return
    endif

    call c_f_pointer(levelctxF, workF)
    call c_f_pointer(levelctxG, workG)

    wkF => workF%wk
    wkG => workG%wk

    wkG = qG
    call fftw_execute_dft(workG%ffft, wkG, wkG)
    wkG = wkG / nvarG

    wkF = 0
    do j = 1, nvarG/2
       do i = 1, nvarG/2
          wkF(i,j)       = wkG(i,j)
       end do
       do i = 0, nvarG/2-2
          wkF(nvarF-i,j) = wkG(nvarG-i,j)
       end do
    end do
    do j = 0, nvarG/2-2
       do i = 1, nvarG/2
          wkF(i,nvarF-j)       = wkG(i,nvarG-j)
       end do
       do i = 0, nvarG/2-2
          wkF(nvarF-i,nvarF-j) = wkG(nvarG-i,nvarG-j)
       end do
    end do

    call fftw_execute_dft(workF%ifft, wkF, wkF)

    qF = real(wkF)
  end subroutine interp2


  subroutine interpolate(qFp, qGp, levelF, levelctxF, levelG, levelctxG,t)
    type(c_ptr), intent(in), value :: qFp, qGp, levelctxF, levelctxG
    integer,     intent(in)        :: levelF, levelG
    real(pfdp),  intent(in) :: t

    real(pfdp), pointer :: qF(:), qG(:), qF2(:,:), qG2(:,:)

    select case(dim)
    case(1)
       qF => array1(qFp)
       qG => array1(qGp)
       call interp1(qF, qG, levelctxF, levelctxG)
    case(2)
       qF2 => array2(qFp)
       qG2 => array2(qGp)
       call interp2(qF2, qG2, levelctxF, levelctxG)
    case(3)
       stop
    end select
  end subroutine interpolate

  subroutine restrict(qFp, qGp, levelF, levelctxF, levelG, levelctxG,t)
    type(c_ptr), intent(in), value :: qFp, qGp, levelctxF, levelctxG
    integer,     intent(in)        :: levelF, levelG
    real(pfdp),  intent(in) :: t

    real(pfdp), pointer :: qF(:), qG(:), qF2(:,:), qG2(:,:)

    integer :: nvarF, nvarG, xrat

    qF => array1(qFp)
    qG => array1(qGp)
    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG
    qG = qF(::xrat)
  end subroutine restrict
end module transfer
