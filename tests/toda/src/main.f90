program pfasst_rttddft
  use pf_mod_mpi
  integer :: rank, err

  call mpi_init(err)
  if (err .ne. 0) stop "ERROR: Can't initialize MPI."
  call mpi_comm_rank(MPI_COMM_WORLD, rank, err)

  call rttddft()
  call mpi_finalize(err)

contains
  subroutine rttddft()
      use pfasst     !< library module defining highest-level pfasst object
      use pf_mod_mpi !< library module for mpi-related business

      use probin     !< should be library module for reading/parsing problem parameters

      use factory    !< prog-specified module containing solution type information
      use sweeper    !< prog-specified module containing sweeper information
      use hooks      !< prog-specified module containing program hooks

      implicit none

      type(pf_pfasst_t)  :: pf          !< pfasst data structure
      type(pf_comm_t)    :: comm        !< mpi communicator

      type(zndarray) :: dmat_t0, dmat_tfinal

      character(256) :: probin_fname       !<  file name for input
      integer    :: err, l
      real(pfdp) :: start, finish

      probin_fname = "probin.nml"
      if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
      call probin_init(probin_fname)

      call pf_mpi_create(comm, MPI_COMM_WORLD)
      call pf_pfasst_create(pf, comm, fname=probin_fname)

!---- Create all the things -------------------------------------------------------

      do l = 1, pf%nlevels
          allocate(pf%levels(l)%shape(1))

          allocate(magpicard_context::pf%levels(l)%ulevel)
          allocate(zndarray_factory::pf%levels(l)%ulevel%factory)
          allocate(magpicard_sweeper_t::pf%levels(l)%ulevel%sweeper)

          call initialize_magpicard_sweeper(pf%levels(l)%ulevel%sweeper, l, pf%qtype, &
               pf%debug, pf%levels(l)%shape, pf%levels(l)%nvars)

          if (pf%qtype == 5) then
            pf%levels(l)%nnodes = nnodes(l)+2
          else
            pf%levels(l)%nnodes = nnodes(l)
          endif
          pf%levels(l)%nsweeps = nsweeps(l)
          pf%levels(l)%nsweeps_pred = nsweeps_pred(l)
      end do

      print *,'Initializing mpi and pfasst...'
      call pf_mpi_setup(comm, pf, err)
      call pf_pfasst_setup(pf)

      call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
      if (save_solutions) call pf_add_hook(pf, -1, PF_POST_ITERATION, save_solution)

      call zndarray_build(dmat_t0, pf%levels(pf%nlevels)%shape(1))
      call zndarray_build(dmat_tfinal, pf%levels(pf%nlevels)%shape(1))
      call initial(dmat_t0)

      call mpi_barrier(MPI_COMM_WORLD, err)

      start = MPI_Wtime()
      print*, 'Running pfasst...'
      call pf_pfasst_run(pf, dmat_t0, dt, 0.0_pfdp, nsteps, dmat_tfinal)

      finish = MPI_Wtime()
      print '("processor ", i2, " returned from pfasst CPU time: ", f16.10, " seconds")', &
           pf%rank, finish-start

      call mpi_barrier(MPI_COMM_WORLD, err)

      if(pf%rank == comm%nproc-1) then
        call dmat_tfinal%write_to_disk('final_solution') !necessary for pfasst.py
        if (pf%debug) call dmat_tfinal%eprint() !only for debug purpose
      endif

      call zndarray_destroy(dmat_t0)
      call zndarray_destroy(dmat_tfinal)

      print *,'destroying pf'
      call pf_pfasst_destroy(pf)

      print *,'destroying MPI'
      call pf_mpi_destroy(comm)

    end subroutine rttddft

  end program pfasst_rttddft
