program main
  use pfasst
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use pf_mod_ndarray

  use feval
  use hooks
  use probin
  use solutions
  use transfer

  implicit none

  type(pf_pfasst_t) :: pf
  type(pf_comm_t)   :: comm

  type(ndarray),      target :: q1
  type(pf_sweeper_t), target :: sweeper
  type(pf_encap_t),   target :: encap

  integer        :: ierror, l
  character(256) :: probin_fname


  !
  ! read options
  !

  if (command_argument_count() == 1) then
     call get_command_argument(1, value=probin_fname)
  else
     probin_fname = "probin.nml"
  end if
  call probin_init(probin_fname)


  !
  ! initialize mpi
  !

  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst
  !

  call ndarray_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  select case(problem)
  case (PROB_WAVE)
     call pf_explicit_create(sweeper, f1eval2)
  case (PROB_SHEAR)
     call pf_imex_create(sweeper, f1eval2, f2eval2, f2comp2)
  case default
     call pf_imex_create(sweeper, f1eval1, f2eval1, f2comp1)
  end select
  call pf_pfasst_create(pf, comm, nlevs)

  pf%niters = niters
  pf%qtype  = SDC_GAUSS_LOBATTO ! + SDC_PROPER_NODES
  pf%window = wtype

  pf%abs_res_tol = abs_tol
  pf%rel_res_tol = rel_tol
  

  pf%echo_timings = .true.

  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     allocate(pf%levels(l)%shape(dim))

     if (problem == PROB_WAVE) then
        pf%levels(l)%shape  = [ nvars(l), 2 ]
     else
        pf%levels(l)%shape  = nvars(l)
     end if

     pf%levels(l)%nvars  = product(pf%levels(l)%shape)
     pf%levels(l)%nnodes = nnodes(maxlevs-nlevs+l)

     call feval_create_workspace(pf%levels(l)%ctx, pf%levels(l)%shape(1))

     pf%levels(l)%encap       => encap
     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%sweeper     => sweeper
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  if (pf%rank == 0) then
     print *, 'nlevs:  ', nlevs
     print *, 'nvars:  ', pf%levels(:)%nvars
     print *, 'nnodes: ', pf%levels(:)%nnodes
     print *, 'output: ', output
  end if


  !
  ! run
  !

  call ndarray_create_simple(q1, pf%levels(nlevs)%shape)
  call initial(q1)

  if (problem == PROB_AD) then
     call pf_add_hook(pf, nlevs, PF_POST_ITERATION, echo_error_hook)
  end if

  if (len_trim(output) > 0) then
     call dump_mkdir(output, len_trim(output))
     call pf_add_hook(pf, nlevs, PF_POST_SWEEP, dump_hook)
  end if

  if (nsteps < comm%nproc) then
     nsteps = comm%nproc
  end if

  call pf_pfasst_run(pf, c_loc(q1), dt, 0.0_pfdp, nsteps=nsteps)


  !
  ! cleanup
  !

  deallocate(q1%flatarray)

  do l = 1, nlevs
     call feval_destroy_workspace(pf%levels(l)%ctx)
  end do

  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(ierror)
  call fftw_cleanup()

end program main
