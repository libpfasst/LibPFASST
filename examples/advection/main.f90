!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

program fpfasst
  use pfasst
  use feval
  use hooks
  use encap_array1d
  use transfer
  use pf_mod_mpi, only: MPI_COMM_WORLD

  implicit none

  type(pf_pfasst_t)  :: pf
  type(pf_comm_t)    :: comm
  type(pf_sweeper_t) :: sweeper
  type(pf_encap_t)   :: encap
  integer            :: ierror, nlevs, nvars(3), nnodes(3), l
  double precision   :: dt

  type(array1d), pointer :: q0
  type(c_ptr) :: q0p

  !!!! initialize mpi
  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !!!! initialize pfasst
  nvars  = [ 32,64,128 ]
  nnodes = [ 3, 5, 9 ]
  dt     = 0.01_pfdp
  nlevs  = 3

  call array1d_create(encap)

  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)
  call pf_pfasst_create(pf, comm, encap, sweeper, nlevs, nvars(1:nlevs), nnodes(1:nlevs))

  pf%niters  = 12
  pf%qtype   = 1

  pf%echo_timings = .false.
  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, pf%nlevels
     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
  end do

  call pf_mpi_setup(comm, pf)
  call setup(pf)

  call add_hook(pf, nlevs, PF_POST_ITERATION, echo_error)

  ! do l = 1, nlevs
  !    call add_hook(pf, l, PF_POST_PREDICTOR, echo_error)
  ! end do


  !!!! initialize advection/diffusion
  call feval_init(size(nvars), nvars)
  call create(q0p, nlevs, .false., nvars(nlevs), [0], c_null_ptr)

  call c_f_pointer(q0p, q0)
  call initial(q0)

  !!!! run
  call pfasst_run(pf, q0p, dt, 0.0_pfdp, 2*comm%nproc)


  !!!! done
  call destroy(q0p)
  call destroy(pf)
  call pf_mpi_destroy(comm)
  call feval_finalize()
  call mpi_finalize(ierror)

end program fpfasst
