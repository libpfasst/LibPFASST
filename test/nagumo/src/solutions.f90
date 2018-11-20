module solutions
  use pf_mod_dtype
  use pf_mod_ndarray_oc
  use probin, only : sizex
  implicit none

contains


  ! Set initial condition.
  subroutine initial(q0, t0, tend)
    type(ndarray_oc), intent(inout) :: q0
    real(pfdp),       intent(in)    :: t0, tend
    integer             :: nvars, i
    real(pfdp)          :: dx, x
      
    q0%yflatarray = 0.0_pfdp
    nvars = size(q0%yflatarray)
    dx = sizex/dble(nvars)
    do i = 1, nvars
       x = dble(i-1)*dx
        if( x >= 9.0 .and. x <= 11.0 ) q0%yflatarray(i) = 1.2*sqrt(3.0)
!        if( x < 10.0 ) q0%yflatarray(i) = -1.2*sqrt(3.0)
!        if( x > 10.0 ) q0%yflatarray(i) =  1.2*sqrt(3.0)
    end do

    !print *, 'initial', t0, tend
    !call exact_y(t0, size(q0%yflatarray), q0%yflatarray)
    q0%pflatarray = 0.0_pfdp ! for adjoint, all zero terminal condition
    !call exact_p(tend, size(q0%pflatarray), q0%pflatarray)

  end subroutine initial


end module solutions
