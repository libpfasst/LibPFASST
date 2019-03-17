program pfasst_imk
  use pf_mod_mpi
  integer :: rank, err

  call mpi_init(err)
  if (err .ne. 0) stop "ERROR: Can't initialize MPI."
  call mpi_comm_rank(MPI_COMM_WORLD, rank, err)

  call simulate()
  if(rank == 0) print *,'shutting down mpi'

  call mpi_finalize(err)
  if(rank == 0) print *,'all done - bye bye'
  

contains
  subroutine simulate()
      use pfasst     !< library module defining highest-level pfasst object
      use pf_mod_mpi !< library module for mpi-related business

      use probin     !< should be library module for reading/parsing problem parameters

      use mod_zmkpair    !< prog-specified module containing solution type information
      use sweeper    !< prog-specified module containing sweeper information
      use hooks      !< prog-specified module containing program hooks
      use utils

      implicit none

      type(pf_pfasst_t)  :: pf          !< pfasst data structure
      type(pf_comm_t)    :: comm        !< mpi communicator

      type(zmkpair) :: q0, qend

      character(256) :: probin_fname       !<  file name for input
      integer    :: err, l
      real(pfdp) :: start, finish

      probin_fname = "probin.nml"
      if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
      call probin_init(probin_fname)

      call pf_mpi_create(comm, MPI_COMM_WORLD)
      call pf_pfasst_create(pf, comm, fname=probin_fname)


      !---- Create the levels -------------------------------------------------------
      do l = 1, pf%nlevels
          allocate(imk_context::pf%levels(l)%ulevel)
          allocate(zmkpair_factory::pf%levels(l)%ulevel%factory)

          allocate(imk_sweeper_t::pf%levels(l)%ulevel%sweeper)

          call initialize_imk_sweeper(pf%levels(l)%ulevel%sweeper, &
               l, pf%debug, use_sdc, rk, mkrk, pf%qtype, nterms(l))

          call pf_level_set_size(pf,l,[nparticles])
          pf%levels(l)%mpibuflen = nparticles * nparticles * 2
      end do

      if(pf%rank == 0) print *,'Initializing mpi and pfasst...'
      call pf_pfasst_setup(pf)
!      call pf_add_hook(pf, -1, PF_PRE_ITERATION, echo_residual)
!      call pf_add_hook(pf, -1, PF_PRE_SWEEP, echo_residual)      
!      call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
      call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_residual)
      !      if (save_solutions) call pf_add_hook(pf, -1, PF_POST_ITERATION, save_solution)
      if (save_solutions) call pf_add_hook(pf, -1, PF_POST_CONVERGENCE, save_solution)      


      !>  output the run options 
      call pf_print_options(pf,un_opt=6)
      
      call zmkpair_build(q0, pf%levels(pf%nlevels)%shape(1))
      call zmkpair_build(qend, pf%levels(pf%nlevels)%shape(1))
      call initial(q0)

      call mpi_barrier(MPI_COMM_WORLD, err)

      start = MPI_Wtime()
      if(pf%rank == 0) print*, 'Running pfasst...'
      call pf_pfasst_run(pf, q0, dt, 0.0_pfdp, nsteps, qend)

      finish = MPI_Wtime()
      print '("processor ", i2, " returned from pfasst CPU time: ", f16.10, " seconds")', &
           pf%rank, finish-start

      call mpi_barrier(MPI_COMM_WORLD, err)

      if(pf%rank == comm%nproc-1) then
        call qend%write_to_disk('final_solution') !necessary for pfasst.py
        if (pf%debug) call qend%eprint() !only for debug purpose
      call qend%eprint() !only for debug purpose        
      endif

      call zmkpair_destroy(q0)
      call zmkpair_destroy(qend)

      !> deallocate local sweeper stuff
      do l = 1, pf%nlevels
         call destroy_imk_sweeper(pf%levels(l)%ulevel%sweeper)
      end do
      
      if(pf%rank == 0) print *,'destroying pf'
      call pf_pfasst_destroy(pf)

    end subroutine simulate

  end program pfasst_imk
