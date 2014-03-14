program main
  use pfasst
  use pf_mod_mpi, only: MPI_THREAD_FUNNELED, MPI_COMM_WORLD
  use pf_mod_ndarray

  use feval
  use hooks
  use probin
  use solutions
  use transfer

  implicit none

  type(pf_pfasst_t)        :: pf
  type(pf_comm_t)          :: comm
  type(ndarray),    target :: q1
  type(pf_encap_t), target :: encap

  integer        :: ierror, iprovided, l
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

  call mpi_init_thread(mpi_thread_funneled, iprovided, ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst
  !

  call ndarray_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
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
     else if (problem == PROB_SHEAR) then
        pf%levels(l)%shape  = [ nvars(l), nvars(l) ]
     else
        pf%levels(l)%shape  = nvars(l)
     end if

     pf%levels(l)%nvars  = product(pf%levels(l)%shape)
     pf%levels(l)%nnodes = nnodes(l)

     if (dim == 1 .or. problem == PROB_WAVE) then
        call create_work1(pf%levels(l)%levelctx, pf%levels(l)%shape(1))
     else if (dim == 2) then
        call create_work2(pf%levels(l)%levelctx, pf%levels(l)%shape(1))
     end if

     pf%levels(l)%encap       => encap
     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict

     select case(problem)
     case (PROB_HEAT)
        call pf_implicit_create(pf%levels(l)%sweeper, f2eval1, f2comp1)
     case (PROB_WAVE)
        call pf_explicit_create(pf%levels(l)%sweeper, f1eval1wave)
     case (PROB_SHEAR)
        call pf_imex_create(pf%levels(l)%sweeper, f1eval2, f2eval2, f2comp2)
     case default
        call pf_imex_create(pf%levels(l)%sweeper, f1eval1, f2eval1, f2comp1)
     end select

  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  pf%outdir = output

  ! if (pf%rank == 0) then
  !    print *, 'nlevs:  ', nlevs
  !    print *, 'nvars:  ', pf%levels(:)%nvars
  !    print *, 'nnodes: ', pf%levels(:)%nnodes
  !    print *, 'output: ', output
  ! end if


  !
  ! run
  !

  call ndarray_create_simple(q1, pf%levels(nlevs)%shape)
  call initial(q1)

  if (problem == PROB_AD) then
     call pf_add_hook(pf, nlevs, PF_POST_ITERATION, echo_error_hook)
  end if

  if (len_trim(output) > 0) then
     call ndarray_mkdir(output, len_trim(output))
     ! call pf_add_hook(pf, nlevs, PF_POST_SWEEP, ndarray_dump_hook)
     call pf_add_hook(pf, -1, PF_POST_SWEEP, ndarray_dump_hook)
  end if

  if (dim == 1) then
     call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_residual_hook)
  end if

  if (nsteps < comm%nproc) then
     nsteps = comm%nproc
  end if

  call pf_print_options(pf)

  call mpi_barrier(mpi_comm_world, ierror)
  call pf_pfasst_run(pf, c_loc(q1), dt, 0.0_pfdp, nsteps=nsteps)


  !
  ! cleanup
  !

  deallocate(q1%flatarray)

  do l = 1, nlevs
     if (dim == 1) then
        call destroy_work1(pf%levels(l)%levelctx)
     else
        ! call destroy_work2(pf%levels(l)%levelctx)
     end if
  end do

  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(ierror)
  call fftw_cleanup()

end program main