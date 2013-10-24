!
! Copyright (c) 2013, Matthew Emmett.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use iso_c_binding
  use encap
  implicit none
contains

  subroutine interpolate(qFptr, qGptr, levelF, ctxF, levelG, ctxG)
    use feval
    type(c_ptr), intent(in), value :: qFptr, qGptr, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    type(carray4), pointer :: qF, qG
    integer :: i, j, k, ii, jj, kk, c, nG, nF

    call c_f_pointer(qFptr, qF)
    call c_f_pointer(qGptr, qG)

    nG = qG%shape(1)
    nF = qG%shape(1)

    !$omp parallel workshare
    qF%array = 0.0d0
    !$omp end parallel workshare

    !$omp parallel do private(i,j,k,ii,jj,kk,c)
    do k = 1, nG
       if (k <= nG/2) then
          kk = k
       else if (k > nG/2+1) then
          kk = nF + nG/2 - k
       else
          cycle
       end if

       do j = 1, nG
          if (j <= nG/2) then
             jj = j
          else if (j > nG/2+1) then
             jj = nF + nG/2 - j
          else
             cycle
          end if

          do i = 1, nG
             if (i <= nG/2) then
                ii = i
             else if (i > nG/2+1) then
                ii = nF + nG/2 - i
             else
                cycle
             end if

             do c = 1, 3
                qF%array(ii, jj, kk, c) = qG%array(i, j, k, c)
             end do
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine interpolate

  subroutine restrict(qFptr, qGptr, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qFptr, qGptr, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    integer :: i, j, k, ii, jj, kk, c, nG, nF
    type(carray4), pointer :: qF, qG

    call c_f_pointer(qFptr, qF)
    call c_f_pointer(qGptr, qG)

    nG = qG%shape(1)
    nF = qG%shape(1)

    !$omp parallel workshare
    qF%array = 0.0d0
    !$omp end parallel workshare

    !$omp parallel do private(i,j,k,ii,jj,kk,c)
    do k = 1, nG
       if (k <= nG/2) then
          kk = k
       else if (k > nG/2+1) then
          kk = nF + nG/2 - k
       else
          cycle
       end if

       do j = 1, nG
          if (j <= nG/2) then
             jj = j
          else if (j > nG/2+1) then
             jj = nF + nG/2 - j
          else
             cycle
          end if

          do i = 1, nG
             if (i <= nG/2) then
                ii = i
             else if (i > nG/2+1) then
                ii = nF + nG/2 - i
             else
                cycle
             end if

             do c = 1, 3
                qG%array(i, j, k, c) = qF%array(ii, jj, kk, c)
             end do
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine restrict
end module transfer
