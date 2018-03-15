module utils
  use pf_mod_dtype
  use probin
  use factory

  implicit none

  complex(pfdp), parameter :: &
       z0 = (0.0_pfdp, 0.0_pfdp), &
       z1 = (1.0_pfdp, 0.0_pfdp), &
       zm1 = (-1.0_pfdp, 0.0_pfdp), &
       zi = (0.0_pfdp, 1.0_pfdp), &
       zmi = (0.0_pfdp, -1.0_pfdp)

 contains
  !> Set initial condition.
  subroutine initial(L)

    class(pf_encap_t), intent(inout) :: L
    class(zndarray), pointer :: L_p

    L_p => cast_as_zndarray(L)

    call exact(L_p, 0.0_pfdp)

    nullify(L_p)
  end subroutine initial

  !> Set exact solution
  !! Exact solution mat(u(t)) = mat(U(t,t0))*mat(u(t0))
  !! U(t,t0) = e^(-iH(t-t0)),
  !! project U into basis of eigenstates of H
  !! aka diagonalize H, then
  !! U(t,t0) = sum_j(|j>e^(-ie_j(t-t0))<j|)
  subroutine exact(L, t)
    use probin, only: alpha, beta, vab, toda_periodic

    type(zndarray), intent(inout) :: L
    real(pfdp), intent(in)  :: t

    complex(pfdp), allocatable :: time_ev_op(:,:), time_ev_op_dagger(:,:), f(:,:)
    complex(pfdp), allocatable :: u(:,:), uinit(:,:) !> definition of initial state
    real(pfdp), allocatable :: q(:), p(:)
    real(pfdp) :: tol=1d-15, alpha_toda
    integer :: dim, i

    dim = L%dim
    allocate(time_ev_op(dim, dim), &
          time_ev_op_dagger(dim, dim), &
          f(dim, dim), &
          u(dim, dim), &
          uinit(dim, dim), q(dim), p(dim))

    u = z0
    uinit = z0
    time_ev_op = z0
    time_ev_op_dagger = z0
    f = z0
    q = 0.0_pfdp
    p = 0.0_pfdp

       ! See Zanna thesis on "On the Numerical Solution of Isospectral Flows"
    if (toda_periodic .eqv. .true.) then
       q = 0
    else
       do i = 1, dim
           q(i) = i - dim/2 - 1 ! everyone assumes a unique position on the number line
       enddo
    endif

    do i = 1, dim
       if (i <= dim/2) then
          p(i) = 4
       elseif (i == dim/2+1) then
          p(i) = 0
       else
          p(i) = -4
       endif
       u(i,i) = p(i) ! everyone gets a different initial momentum
    enddo

    do i = 1, dim-1
       alpha_toda = exp(-0.5_pfdp * (q(i+1) - q(i)))
       if (i == 1 .and. toda_periodic .eqv. .true.) then
          u(i, dim) = alpha_toda
          u(dim, i) = alpha_toda
       endif
       u(i, i+1) = alpha_toda
       u(i+1, i) = alpha_toda
    enddo

    u(:,:) = 0.5_pfdp * u(:,:)
    L%array = u

    deallocate(time_ev_op, time_ev_op_dagger, f, u, uinit, q, p)

  end subroutine exact
end module utils
