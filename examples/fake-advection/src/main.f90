!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

program fpfasst
  use pfasst
  use feval
  use hooks
  use transfer

  implicit none

  type(pf_pfasst_t),  pointer :: pf, pf0
  type(pf_comm_t),    target  :: comm
  type(pf_sweeper_t), target  :: sweeper
  type(pf_encap_t),   target  :: encap

  integer            :: nprocs, nlevs, nvars(3), nnodes(3), l, p
  double precision   :: dt

  type(ndarray), pointer :: q0

  !
  ! initialize pfasst
  !

  nprocs = 8
  nvars  = [ 32, 64, 128 ]
  nnodes = [ 2, 3, 5 ]
  dt     = 0.1_pfdp
  nlevs  = 3

  call ndarray_encap_create(encap)
  call pf_fake_create(comm, nprocs)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)

  do p = 1, nprocs
     allocate(pf)
     call pf_pfasst_create(pf, comm, nlevs)

     pf%rank   = p - 1
     pf%niters = 12
     pf%qtype  = SDC_GAUSS_LOBATTO + SDC_PROPER_NODES

     pf%echo_timings = .false.

     pf%window      = PF_WINDOW_RING
     pf%abs_res_tol = 1.d-8

     if (nlevs > 1) then
        pf%levels(1)%nsweeps = 2
     end if

     do l = 1, nlevs
        pf%levels(l)%nvars  = nvars(3-nlevs+l)
        pf%levels(l)%nnodes = nnodes(3-nlevs+l)
        pf%levels(l)%interpolate => interpolate
        pf%levels(l)%restrict    => restrict
        pf%levels(l)%encap       => encap
        pf%levels(l)%sweeper     => sweeper

        allocate(pf%levels(l)%shape(1))
        pf%levels(l)%shape(1)    = nvars(l)
     end do

     call pf_fake_setup(comm, pf)
     call pf_pfasst_setup(pf)

     call pf_add_hook(pf, nlevs, PF_POST_ITERATION, echo_error)
     ! call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)

     if (pf%rank == 0) then
        print *, 'nvars: ', pf%levels(:)%nvars
        print *, 'nnodes:', pf%levels(:)%nnodes
     end if

  end do

  !
  ! create workspaces on first processor and attach them to all fake processors
  !
  call c_f_pointer(comm%pfs(1), pf0)
  do l = 1, nlevs
     call feval_create_workspace(pf0%levels(l)%ctx, pf0%levels(l)%nvars)
  end do

  do p = 2, nprocs
     call c_f_pointer(comm%pfs(p), pf)
     do l = 1, nlevs
        pf%levels(l)%ctx = pf0%levels(l)%ctx
     end do
  end do


  !
  ! run
  !
  allocate(q0)
  call ndarray_create_simple(q0, [ pf%levels(nlevs)%nvars ])
  call initial(q0)

  call pf_fake_run(pf, c_loc(q0), dt, nsteps=2*nprocs)


  !
  ! cleanup
  !
  call ndarray_destroy(c_loc(q0))

  do l = 1, nlevs
     call feval_destroy_workspace(pf0%levels(l)%ctx)
  end do

  call pf_imex_destroy(sweeper)
  call pf_fake_destroy(comm)
  call fftw_cleanup()

end program fpfasst
