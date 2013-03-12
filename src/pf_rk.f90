!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module pf_mod_rk
  use pf_mod_dtype
  implicit none
contains

  ! Run in serial using SSPRK3
  subroutine rk3_run(pf, q0, dt, tend, nsteps, qend)
    use pf_mod_timer
    use pf_mod_hooks
    use feval, only: eval_f1

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_encap_t),  intent(in)    :: q0
    double precision,  intent(in)    :: dt, tend
    integer,           intent(in), optional    :: nsteps
    type(pf_encap_t),  intent(inout), optional :: qend

    type(pf_level_t), pointer :: F
    type(pf_encap_t) :: y0F, y0newF, fF

    integer    :: i, step, nblocks
    real(pfdp) :: t0

    F => pf%levels(1)
    pf%state%dt = dt

    call start_timer(pf, TTOTAL)

    call create(y0F, F%level, .false., F%nvars, F%shape, F%ctx)
    call create(y0newF, F%level, .false., F%nvars, F%shape, F%ctx)
    call create(Ff, F%level, .true., F%nvars, F%shape, F%ctx)

    !!!! set initial conditions

    call copy(y0F, q0)

    if(present(nsteps)) then
       nblocks = nsteps
    else
       nblocks = int(ceiling(tend/dt))
    end if

    !!!! time step loop
    do i = 1, nblocks

       step = i-1
       t0   = step * dt

       pf%state%step = step
       pf%state%t0   = t0

       ! y0newF = y0F + (1.0d0/3) * dt * fF
       call eval_f1(y0F, t0, F%level, F%ctx, fF)
       call copy(y0newF, y0F)
       call axpy(y0newF, (1.0d0/3) * dt, fF)

       ! y0newF = y0F + (1.0d0/2) * dt * fF
       call eval_f1(y0newF, t0, F%level, F%ctx, fF)
       call copy(y0newF, y0F)
       call axpy(y0newF, (1.0d0/2) * dt, fF)

       ! y0F = y0F + dt*fF
       call eval_f1(y0newF, t0, F%level, F%ctx, fF)
       call axpy(y0F, dt, fF)

       call copy(F%qend, y0F)
       call call_hooks(pf, F%level, PF_POST_STEP)

    end do

    call destroy(y0F)
    call destroy(y0newF)
    call destroy(fF)

    call end_timer(pf, TTOTAL)

    if (present(qend)) then
       F => pf%levels(1)
       call copy(qend, F%qend)
    end if

  end subroutine rk3_run

  subroutine rk4_run(pf, q0, dt, tend, nsteps, qend)
    use pf_mod_timer
    use pf_mod_hooks
    use feval, only: eval_f1

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_encap_t),  intent(in)    :: q0
    double precision,  intent(in)    :: dt, tend
    integer,           intent(in), optional    :: nsteps
    type(pf_encap_t),  intent(inout), optional :: qend

    integer    :: i, step, nblocks
    real(pfdp) :: t0

    type(pf_level_t), pointer :: F
    type(pf_encap_t) :: y0F, y0newF, f0, f1, f2, f3


    F => pf%levels(1)
    pf%state%dt = dt

    call start_timer(pf, TTOTAL)

    call create(y0F, F%level, .false., F%nvars, F%shape, F%ctx)
    call create(y0newF, F%level, .false., F%nvars, F%shape, F%ctx)
    call create(f0, F%level, .true., F%nvars, F%shape, F%ctx)
    call create(f1, F%level, .true., F%nvars, F%shape, F%ctx)
    call create(f2, F%level, .true., F%nvars, F%shape, F%ctx)
    call create(f3, F%level, .true., F%nvars, F%shape, F%ctx)

    !!!! set initial conditions

    call copy(y0F, q0)

    if(present(nsteps)) then
       nblocks = nsteps
    else
       nblocks = int(ceiling(tend/dt))
    end if

    !!!! time step loop
    do i = 1,nblocks

       step = i-1
       t0   = step * dt

       pf%state%step = step
       pf%state%t0   = t0

       ! y0newF = y0F + 0.5d0 * dt * f0
       call eval_f1(y0F, t0, F%level, F%ctx, f0)
       call copy(y0newF, y0F)
       call axpy(y0newF, 0.5d0*dt, f0)

       ! y0newF = y0F + 0.5d0 * dt * f1
       call eval_f1(y0newF, t0, F%level, F%ctx, f1)
       call copy(y0newF, y0F)
       call axpy(y0newF, 0.5d0*dt, f1)

       ! y0newF = y0F + dt * f2
       call eval_f1(y0newF, t0, F%level, F%ctx, f2)
       call copy(y0newF, y0F)
       call axpy(y0newF, dt, f2)

       ! y0F = y0F + dt/6 * (f0 + 2.0d0 * (f1 + f2) + f3)
       call eval_f1(y0newF, t0, F%level, F%ctx, f3)
       call axpy(y0F, dt/6, f0)
       call axpy(y0F, dt/3, f1)
       call axpy(y0F, dt/3, f2)
       call axpy(y0F, dt/6, f3)

       call copy(F%qend, y0F)
       call call_hooks(pf, F%level, PF_POST_STEP)

    end do

    call destroy(y0F)
    call destroy(y0newF)
    call destroy(f0)
    call destroy(f1)
    call destroy(f2)
    call destroy(f3)

    call end_timer(pf, TTOTAL)

    if (present(qend)) then
       F => pf%levels(1)
       call copy(qend, F%qend)
    end if

  end subroutine rk4_run

  subroutine ark4_run(pf, q0, dt, tend, nsteps, qend)
    use pf_mod_timer
    use pf_mod_hooks
    use feval, only: eval_f1, eval_f2, comp_f2

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_encap_t),  intent(in)    :: q0
    double precision,  intent(in)    :: dt, tend
    integer,           intent(in), optional    :: nsteps
    type(pf_encap_t),  intent(inout), optional :: qend

    integer    :: i, j, k, step, nblocks
    real(pfdp) :: t0

    type(pf_level_t), pointer :: F
    type(pf_encap_t) :: qRK(0:5), fRK(0:1,0:5), rhs

    real(pfdp) :: be(0:5), ae(0:5,0:5), ai(0:5,0:5)

    be = 0.0d0
    be(0) = 82889.0d0/524892.0d0
    be(2) = 15625.0d0/83664.0d0
    be(3) = 69875.0d0/102672.0d0
    be(4) =-2260.0d0/8211.0d0
    be(5) = 1.0d0/4.0d0

    ae = 0.0d0
    ae(1,0) = 1.0d0/2.0d0
    ae(2,0) = 13861.0d0/62500.0d0
    ae(2,1) = 6889.0d0/62500.0d0
    ae(3,0) =-116923316275.0d0/2393684061468.0d0
    ae(3,1) =-2731218467317.0d0/15368042101831.0d0
    ae(3,2) = 9408046702089.0d0/11113171139209.0d0
    ae(4,0) =-451086348788.0d0/2902428689909.0d0
    ae(4,1) =-2682348792572.0d0/7519795681897.0d0
    ae(4,2) = 12662868775082.0d0/11960479115383.0d0
    ae(4,3) = 3355817975965.0d0/11060851509271.0d0
    ae(5,0) = 647845179188.0d0/3216320057751.0d0
    ae(5,1) = 73281519250.0d0/8382639484533.0d0
    ae(5,2) = 552539513391.0d0/3454668386233.0d0
    ae(5,3) = 3354512671639.0d0/8306763924573.0d0
    ae(5,4) = 4040.0d0/17871.0d0

    ai = 0.0d0
    ai(1,0) = 1.0d0/4.0d0
    ai(1,1) = 1.0d0/4.0d0
    ai(2,0) = 8611.0d0/62500.0d0
    ai(2,1) =-1743.0d0/31250.0d0
    ai(2,2) = 1.0d0/4.0d0
    ai(3,0) = 5012029.0d0/34652500.0d0
    ai(3,1) =-654441.0d0/2922500.0d0
    ai(3,2) = 174375.0d0/388108.0d0
    ai(3,3) = 1.0d0/4.0d0
    ai(4,0) = 15267082809.0d0/155376265600.0d0
    ai(4,1) =-71443401.0d0/120774400.0d0
    ai(4,2) = 730878875.0d0/902184768.0d0
    ai(4,3) = 2285395.0d0/8070912.0d0
    ai(4,4) = 1.0d0/4.0d0
    ai(5,0) = 82889.0d0/524892.0d0
    ai(5,2) = 15625.0d0/83664.0d0
    ai(5,3) = 69875.0d0/102672.0d0
    ai(5,4) =-2260.0d0/8211.0d0
    ai(5,5) = 1.0d0/4.0d0

    F => pf%levels(1)
    pf%state%dt = dt

    call start_timer(pf, TTOTAL)

    do k = 0, 5
       do i = 0, 1
          call create(fRK(i,k), F%level, .true., F%nvars, F%shape, F%ctx)
       end do
       call create(qRK(k), F%level, .false., F%nvars, F%shape, F%ctx)
    end do

    call create(rhs, F%level, .false., F%nvars, F%shape, F%ctx)

    !!!! set initial conditions

    call copy(qRK(0), q0)

    if(present(nsteps)) then
       nblocks = nsteps
    else
       nblocks = int(ceiling(tend/dt))
    end if

    !!!! time step loop
    do i = 1, nblocks

       step = i-1
       t0   = step * dt

       pf%state%step = step
       pf%state%t0   = t0

       call eval_f1(qRK(0), t0, F%level, F%ctx, fRK(0,0))
       call eval_f2(qRK(0), t0, F%level, F%ctx, fRK(1,0))

       do k = 1, 5
          call copy(rhs, qRK(0))
          do j = 0, k-1
             call axpy(rhs, dt*ae(k,j), fRK(0,j))
             call axpy(rhs, dt*ai(k,j), fRK(1,j))
          end do

          call comp_f2(qRK(k), t0, ai(k,k)*dt, rhs, F%level, F%ctx, fRK(1,k))
          call eval_f1(qRK(k), t0, F%level, F%ctx, fRK(0,k))
       end do

       do j = 0, 5
          call axpy(qRK(0), be(j)*dt, fRK(0,j))
          call axpy(qRK(0), be(j)*dt, fRK(1,j))
       end do

       call copy(F%qend, qRK(0))
       call call_hooks(pf, F%level, PF_POST_STEP)

    end do

    !!!! done

    do k = 0, 5
       do i = 0, 1
          call destroy(fRK(i,k))
       end do
       call destroy(qRK(k))
    end do

    call destroy(rhs)

    call end_timer(pf, TTOTAL)

    if (present(qend)) then
       F => pf%levels(1)
       call copy(qend, F%qend)
    end if

  end subroutine ark4_run

  subroutine sdc_integrate(qSDC,fSDC, dt, F, fintSDC)
  ! Compute SDC integral
    type(pf_level_t),  intent(in)    :: F
    type(pf_encap_t),  intent(in)    :: qSDC(F%nnodes,1) ! Solution
    type(pf_encap_t),  intent(in)    :: fSDC(F%nnodes,1) ! Function values
    type(pf_encap_t),  intent(inout) :: fintSDC(F%nnodes-1)    ! Integrals of f
    real(pfdp),        intent(in)    :: dt

    integer :: n, m, p

    print *,' Warning, this routine sdc_integrate is not implemented for RK'
  end subroutine sdc_integrate

end module pf_mod_rk
