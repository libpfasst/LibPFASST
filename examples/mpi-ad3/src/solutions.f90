module solutions
  use pf_mod_dtype
  use pf_mod_ndarray
  use probin
  implicit none

contains

  ! Set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0
    integer             :: shp2(2)
    real(pfdp), pointer :: y2(:,:)

    ! select case(problem)
    ! case (PROB_AD)
       call exact(0.0_pfdp, size(q0%flatarray), q0%flatarray)
    ! case (PROB_HEAT)
    !    call gaussian(q0%flatarray)
    ! case (PROB_VB)
    !    call gaussian(q0%flatarray)
    ! case (PROB_KS)
    !    call sinusoidal(q0%flatarray)
    ! case (PROB_WAVE)
    !    shp2 = q0%shape
    !    call c_f_pointer(q0%aptr, y2, shp2)
    !    y2 = 0
    !    call gaussian(y2(:, 1))
    ! case (PROB_SHEAR)
    !    shp2 = q0%shape
    !    call c_f_pointer(q0%aptr, y2, shp2)
    !    call vortex_sheets(y2)
    ! case default
    !    stop "ERROR: Unknown problem type (initial)."
    ! end select
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gaussian(q)
    real(pfdp), intent(inout) :: q(:)
    integer :: nvars, ii, i
    double precision :: x

    nvars = size(q)

    q = 0.0d0
    do ii = -3, 3
       do i = 1, nvars
          x = Lx * (dble(i-nvars/2-1)/dble(nvars) + ii)
          q(i) = q(i) + dexp(-x**2/sigma)
       end do
    end do
  end subroutine gaussian

  subroutine sinusoidal(q)
    real(pfdp), intent(inout) :: q(:)
    integer :: nvars, i
    double precision :: x,dx

    nvars = size(q)

    do i = 1, nvars
       x = dble(i-1)/dble(nvars)
       q(i) = 0.1*dsin(3*two_pi*x) &
            + 0.2*dsin(4*two_pi*x) &
            + 0.3*dsin(7*two_pi*x)
    end do
  end subroutine sinusoidal

  subroutine exact_ode(t, nvars, yex)
    use probin, only:  nu, Lx, spatial_order, kfreq
    real(pfdp), intent(in)  :: t
    integer,    intent(in)  :: nvars
    real(pfdp), intent(out) :: yex(nvars)

    integer :: i
    real(pfdp) :: dx,x,rho

    yex = 0.0_pfdp

    dx = 1.0_pfdp/dble(nvars)
    rho = (2.0_pfdp-2.0_pfdp*cos(kfreq*pi*dx))/(dx*dx)
    do i = 1, nvars
       x = dx*dble(i-1) - t*v
       yex(i) = yex(i) + dsin(kfreq*pi*x)*dexp(-rho*nu*t)
    end do

  end subroutine exact_ode
  subroutine exact(t, nvars, yex)
    real(pfdp), intent(in)  :: t
    integer,    intent(in)  :: nvars
    real(pfdp), intent(out) :: yex(nvars)

    integer :: i, ii
    real(pfdp) :: x,dx

    yex = 0.0_pfdp

    dx = 1.0_pfdp/dble(nvars)
    do i = 1, nvars
       x = dx*dble(i-1) - t*v
       yex(i) = yex(i) + dsin(kfreq*pi*x)*dexp(-kfreq*kfreq*pi*pi*nu*t)
    end do

  end subroutine exact


end module solutions
