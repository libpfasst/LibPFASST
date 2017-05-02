! ################################################################
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
! ################################################################

#include <petsc/finclude/petscdef.h>

program main

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ! `main': Main program
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  ! ================================================================
  ! Modules
  ! ================================================================

  ! Fortran
  use iso_c_binding,only : c_ptr

  ! External
  use mpi,only : MPI_COMM_WORLD,MPI_wTime,MPI_MAX,MPI_Comm_free,MPI_DOUBLE_PRECISION
  use petsc,only : PETSC_COMM_WORLD,PETSC_NULL_CHARACTER
  use pfasst,only : pfdp,pf_comm_t,pf_mpi_create,pf_mpi_setup,pf_mpi_destroy,pf_pfasst_t,pf_pfasst_create,pf_pfasst_setup,pf_pfasst_destroy,pf_encap_t,pf_explicitQ_create,pf_implicitQ_create,pf_imexQ_create,pf_pfasst_run

  ! Project
  use mod_pf,only : num_level,num_iter,qType,window,abs_res_tol,rel_res_tol,pipeline_g,pfasst_pred,calc_residual,echo_timing,taui0,out_dir,num_step_time,size_step_time,orig_dom_time,all_num_node_time,all_num_sweep_pred,all_num_sweep_corr,create_pf,destroy_pf
  use mod_mpi,only : max_len_char_,num_dim_space,create_mpi,destroy_mpi,mpi_time_,mpi_space_,mpi_world_,mpi_exit_gracefully
  use mod_ctx,only : num_var_dep,create_ctx,destroy_ctx,all_f_ctx_
  use mod_prob,only : prob,create_prob,destroy_prob,prob_set_ini
  use mod_hook,only : create_hook,destroy_hook
  use mod_deriv,only : create_deriv,destroy_deriv
  use mod_encap,only : create_encap,destroy_encap,encap_create_pf,encap_destroy_pf
  use mod_petsc,only : create_petsc,destroy_petsc
  use mod_feval,only : sweeper,use_LU,create_feval,destroy_feval,feval_f1eval,feval_f2eval,feval_f2comp
  use mod_interp,only : create_interp,destroy_interp
  use mod_transfer,only : create_transfer,destroy_transfer,transfer_interpolate,transfer_restrict

  ! ================================================================
  ! Safety
  ! ================================================================

  implicit none

  ! ================================================================ 
  ! Parameters
  ! ================================================================ 

  integer :: io_unit_ = 9 ! Unit for reading
  integer :: bin_name_len_ = 4 ! Length of string denoting executable (i.e. `len_trim(main) = 4')
  character(len = max_len_char_) :: name_file_input_petsc_stub_ = "arg/petsc.txt" ! Name stub of input file for PETSc solvers
  character(len = max_len_char_) :: name_file_input_pfasst_stub_ = "arg/pfasst.txt" ! File name stub for non-PETSc arguments

  ! ================================================================
  ! Arguments
  ! ================================================================

  ! ================================================================
  ! Locals
  ! ================================================================

  type(pf_comm_t) :: comm_pf
  type(pf_pfasst_t) :: pfasst_pf
  type(pf_encap_t),target :: encap_pf
  type(c_ptr) :: c_encap_tic,c_encap_toc
  real(pfdp) :: tot_runtime,runtime,max_runtime
  integer :: mpi_stat,petsc_stat,i_level,alloc_stat,dealloc_stat
  character(len = max_len_char_) :: bin_dir_name,name_file_input_pfasst,name_file_input_petsc

  ! ================================================================
  ! Work
  ! ================================================================

  ! Global MPI initialization
  call MPI_Init(mpi_stat)

  ! Get executable's location
  call get_command_argument(0,value = bin_dir_name)
  ! Get input file name
  if(command_argument_count() > 0) then
    ! Input file given on command line
    call get_command_argument(1,value = name_file_input_pfasst)
  else
    ! Input file *not* given on command line
    name_file_input_pfasst = bin_dir_name(1:len_trim(bin_dir_name)-bin_name_len_)//trim(name_file_input_pfasst_stub_)
  end if

  ! Input file name for PETSc parameters
  name_file_input_petsc = bin_dir_name(1:len_trim(bin_dir_name)-bin_name_len_)//trim(name_file_input_petsc_stub_)

  ! Define parallel programming world
  call create_mpi(MPI_COMM_WORLD,io_unit_,name_file_input_pfasst)
  ! MPI world communicator for PETSc
  call MPI_Comm_dup(mpi_space_%comm,PETSC_COMM_WORLD,mpi_stat)
  ! Initialize PETSc
  call PetscInitialize(PETSC_NULL_CHARACTER,petsc_stat)

  ! Read PFASST and time domain input arguments
  call create_pf(io_unit_,name_file_input_pfasst)

  ! Create PFASST's MPI time communicator
  call pf_mpi_create(comm_pf,mpi_time_%comm)
  ! Create PFASST's main container
  call pf_pfasst_create(pfasst_pf,comm_pf,num_level,noCmd = .true.)

  ! Complete value-initialization of `pfasst_pf' members
  pfasst_pf%nIters = num_iter ! Number of PFASST iterations
  pfasst_pf%qType = qType ! Legendre (`1`) or Radau (`2`) Gauss-Lobatto nodes in time (use `256+qType` if `all_num_node_time` [see below] requires proper inter-level interpolation and transfer_restriction)
  pfasst_pf%window = window ! Window mode for convergence check
  pfasst_pf%abs_res_tol = abs_res_tol ! Absolute tolerance for iterations
  pfasst_pf%rel_res_tol = rel_res_tol ! Relative tolerance for iterations
  pfasst_pf%pipeline_g = pipeline_g ! Pipelined prediction
  pfasst_pf%pfasst_pred = pfasst_pred ! Predictor activation
  pfasst_pf%calc_residuals = calc_residual ! Calculate residual when entering a hook
  pfasst_pf%echo_timings = echo_timing ! Performance output 
  pfasst_pf%taui0 = taui0 ! Tau correction
  pfasst_pf%outDir = out_dir ! Output directory

  ! Prepare derivative coefficients
  call create_deriv(io_unit_,name_file_input_pfasst)
  ! Set up interpolation weights
  call create_interp(io_unit_,name_file_input_pfasst)
  ! Define degrees of freedom (independent and dependent)
  call create_ctx(io_unit_,name_file_input_pfasst,pfasst_pf)

  ! Define temporal deegrees of freedom for PFASST levels
  do i_level = 1,num_level
    pfasst_pf%levels(i_level)%nNodes = all_num_node_time(i_level) ! Number of `qType' nodes
    pfasst_pf%levels(i_level)%nSweeps_pred = all_num_sweep_pred(i_level) ! Number of predictor sweeps
    pfasst_pf%levels(i_level)%nSweeps = all_num_sweep_corr(i_level) ! Number of main sweeps
  end do ! `i_level'

  ! Define spatial degrees of freedom for PFASST levels
  do i_level = 1,num_level
    ! Create `shape' array for distribution of DoFs
    allocate(pfasst_pf%levels(i_level)%shape(1),stat = alloc_stat)
    if(alloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Allocation of `pfasst_pf%levels(...)%shape' failed.")
    end if ! `alloc_stat .ne. 0'
    
    ! Define `shape' values
    pfasst_pf%levels(i_level)%shape(1) = num_var_dep*all_f_ctx_(i_level)%f_ctx%count(1)
    ! Total number of variables in solution
    pfasst_pf%levels(i_level)%nVars = product(pfasst_pf%levels(i_level)%shape)
  end do ! `i_level'

  ! Create PFASST encapsulation interface to access such routines as `copy', `aXpY', ...
  call encap_create_pf(encap_pf)
  ! Let each PFASST encapsulation point to the encapsulation example created just now
  do i_level = 1,num_level
    ! Assign instance of PFASST encapsulation
    pfasst_pf%levels(i_level)%encap => encap_pf
  end do ! `i_level'

  ! Create level transfer environment
  call create_transfer()

  ! Routines for inter-PFASST level data transfer
  do i_level = 1,num_level
    ! Interpolation routine
    pfasst_pf%levels(i_level)%interpolate => transfer_interpolate
    ! Restriction routine
    pfasst_pf%levels(i_level)%restrict => transfer_restrict
  end do ! `i_level'

  ! Create problem and solver objects
  call create_prob(io_unit_,name_file_input_pfasst,name_file_input_petsc,pfasst_pf)
  ! Prepare output
  call create_hook(io_unit_,name_file_input_pfasst,pfasst_pf)
  ! Create PETSc solvers
  call create_petsc(io_unit_,name_file_input_pfasst,name_file_input_petsc,pfasst_pf)
  ! Create time steppers
  call create_feval(io_unit_,name_file_input_pfasst)
  ! Create sweeper
  select case(num_dim_space)
  case(1)
    call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
  case(2)
    select case(prob)
    case(1)
      do i_level = 1,num_level
        select case(sweeper)
        case(1)
          call pf_explicitQ_create(pfasst_pf%levels(i_level)%sweeper,feval_f1eval)
        case(2)
          call pf_implicitQ_create(pfasst_pf%levels(i_level)%sweeper,feval_f2eval,feval_f2comp,use_LU)
        case(3)
          call pf_imexQ_create(pfasst_pf%levels(i_level)%sweeper,feval_f1eval,feval_f2eval,feval_f2comp,use_LU)
        case default
          call mpi_exit_gracefully(__FILE__,__LINE__,"Sweeper not available.")
        end select ! `sweeper'
      end do ! `i_level'
    case(2)
      call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
    case default
      call mpi_exit_gracefully(__FILE__,__LINE__,"Problem not available.")
    end select ! `prob'
  case default
    call mpi_exit_gracefully(__FILE__,__LINE__,"Number of spatial dimensions not available.")
  end select ! `num_dim_space'

  ! Reserve space for send and receive requests (one for each PFASST level)
  call pf_mpi_setup(comm_pf,pfasst_pf)
  ! Allocate memory for nodal matrices
  call pf_pfasst_setup(pfasst_pf)

  ! Create initial solution
  call create_encap(c_encap_tic,num_level,1,pfasst_pf%levels(num_level)%nVars,pfasst_pf%levels(num_level)%shape,pfasst_pf%levels(num_level)%ctx)
  ! Assign initial data
  call prob_set_ini(c_encap_tic)
  ! Create final solution
  call create_encap(c_encap_toc,num_level,1,pfasst_pf%levels(num_level)%nVars,pfasst_pf%levels(num_level)%shape,pfasst_pf%levels(num_level)%ctx)

  ! Initialize total runtime
  tot_runtime = 0d0

  ! All processes shall start the integration simultaneously
  call MPI_Barrier(mpi_world_%comm,mpi_stat)
    
  ! Start ticking
  runtime = MPI_wTime()

  ! Run PFASST solver
  call pf_pfasst_run(pfasst_pf,c_encap_tic,size_step_time,orig_dom_time,num_step_time,c_encap_toc)

  ! Stop ticking
  tot_runtime = tot_runtime+MPI_wTime()-runtime

  ! Get runtime
  call mpi_reduce(tot_runtime,max_runtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,mpi_world_%size-1,mpi_world_%comm,mpi_stat)

  ! Print runtime
  if(mpi_world_%rank .eq. mpi_world_%size-1) then
    write(6,"(a,a,a,1i0,a,1es24.16e3,a)") "`",trim(__FILE__),":",__LINE__,"': runtime = ",max_runtime," (s)"
    call flush(6)
  end if ! `mpi_world_%rank .eq. mpi_world_%size-1'

  ! Destroy final solution
  call destroy_encap(c_encap_toc)
  ! Destroy initial solution
  call destroy_encap(c_encap_tic)

  ! Destroy time steppers
  call destroy_feval()
  ! Destroy PETSc solvers
  call destroy_petsc()
  ! Clean up output environment
  call destroy_hook()
  ! Free memory occupied by problem its solver objects
  call destroy_prob()

  ! Create level transfer environment
  call destroy_transfer()

  ! Remove PFASST encapsulation interface
  call encap_destroy_pf(encap_pf)

  ! Remove container for spatial degrees of freedom for PFASST levels
  do i_level = 1,num_level
    ! Destroy `shape' array for distribution of DoFs
    deallocate(pfasst_pf%levels(i_level)%shape,stat = dealloc_stat)
    if(dealloc_stat .ne. 0) then
      call mpi_exit_gracefully(__FILE__,__LINE__,"Deallocation of `pfasst_pf%levels(...)%shape' failed.")
    end if ! `dealloc_stat .ne. 0'
  end do ! `i_level'

  ! Destroy spatial PETSc meshes
  call destroy_ctx(pfasst_pf)
  ! Remove containers for interpolation weights
  call destroy_interp()
  ! Erase arrays for derivative coefficients
  call destroy_deriv()

  ! Destroy PFASST's main container
  call pf_pfasst_destroy(pfasst_pf)
  ! Destroy PFASST's MPI time communicator
  call pf_mpi_destroy(comm_pf)

  ! Free memory occupied for definition of time domain specifics
  call destroy_pf()

  ! Finalize PETSc
  call PetscFinalize(petsc_stat)
  ! Free PETSc communicator
  call MPI_Comm_free(PETSC_COMM_WORLD,mpi_stat)
  ! Destroy space occupied by MPI world
  call destroy_mpi()
  ! Global MPI finalization
  call MPI_Finalize(mpi_stat)

  ! ================================================================
  ! Exit
  ! ================================================================

  return

end program main
