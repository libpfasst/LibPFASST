!
! Block MPI advection/diffusion test.
!

program main
  use pfasst
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use spec_ad

  implicit none

  integer, parameter :: nlevs = 3

  type(pf_pfasst_t)           :: pf
  type(pf_comm_t)             :: comm
  type(pf_sweeper_t), target  :: sweeper
  type(pf_encap_t),   target  :: encap
  type(ndarray),      pointer :: q0
  integer                     :: err, nvars(nlevs), nnodes(nlevs), l
  double precision            :: dt


  !
  ! initialize mpi
  !

  call mpi_init(err)
  if (err .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst using three levels
  !

  nvars  = [ 32, 64, 128 ]   ! number of dofs on the time/space levels
  nnodes = [ 2, 3, 5 ]       ! number of sdc nodes on time/space levels
  dt     = 0.1_pfdp

  call ndarray_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)
  call pf_pfasst_create(pf, comm, nlevs)

  pf%niters = 12
  pf%qtype  = SDC_GAUSS_LOBATTO

  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     pf%levels(l)%nvars  = nvars(l)
     pf%levels(l)%nnodes = nnodes(l)

     call feval_create_workspace(pf%levels(l)%ctx, pf%levels(l)%nvars)

     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%encap       => encap
     pf%levels(l)%sweeper     => sweeper

     allocate(pf%levels(l)%shape(1))
     pf%levels(l)%shape(1)    = nvars(l)
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  if (pf%rank == 0) then
     print *, 'nvars: ', pf%levels(:)%nvars
     print *, 'nnodes:', pf%levels(:)%nnodes
  end if


  !
  ! compute initial condition, add hooks, run
  !

  allocate(q0)
  call ndarray_create_simple(q0, [ nvars(3) ])
  call initial(q0)

  call pf_add_hook(pf, nlevs, PF_POST_ITERATION, echo_error)
  call pf_pfasst_run(pf, c_loc(q0), dt, tend=0.d0, nsteps=4*comm%nproc)


  !
  ! cleanup
  !

  call ndarray_destroy(c_loc(q0))

  do l = 1, nlevs
     call feval_destroy_workspace(pf%levels(l)%ctx)
  end do

  call pf_imex_destroy(sweeper)
  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(err)
  call fftw_cleanup()

contains 

  subroutine echo_error(pf, level, state, ctx)
    use spec_ad, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)

    qend => array1(level%qend)

    call exact(state%t0+state%dt, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(qend-yexact))
  end subroutine echo_error

end program main
