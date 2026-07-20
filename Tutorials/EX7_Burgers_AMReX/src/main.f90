!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the PDEs in periodic domains using a
!!  spectral representation of the solution
!! 
!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
program main
  use pfasst         !<  This module has include statements for the main pfasst routines
  use pf_mod_mpi

  integer ::  ierror, gRank, gSize
  interface
     subroutine add_parameters () bind(c)
     end subroutine add_parameters
  end interface

   !> initialize MPI
   print *,'main.f90:initialize mpi'
   call mpi_init(ierror)
   if (ierror /= 0) call pf_stop(__FILE__,__LINE__,'Can not initialize MPI, ierrer',ierror)
   call mpi_comm_rank(MPI_COMM_WORLD, gRank, ierror)
   call mpi_comm_size(MPI_COMM_WORLD, gSize, ierror)
      

   !> call the routine to do PFASST
   call run_pfasst()
  
   !> close mpi
   if (gRank .eq. 0) print *,'main.f90:finalize mpi'
   call mpi_finalize(ierror)
   if (ierror /= 0) call pf_stop(__FILE__,__LINE__,'Can not finalize MPI, ierrer',ierror)
   if (gRank .eq. 0) print *,'main.f90:all done'

contains
  !>  This subroutine implements pfasst to PDEs in spectral space
  subroutine run_pfasst()  
    use pf_my_sweeper       !< Local module for defining sweeper and function evaluations
    use pf_my_level         !< Local module defining the levels and restriction/interpolation
    use pf_mod_AMReX_mfab   !< Libpfasst encapsulation for  AMReX multifabs
    use pf_space_time_comm  !< Local module for space and time mpi communicator
    use hooks               !< Local module for diagnostics and i/o
    use probin              !< Local module reading/parsing problem parameters
    use pf_mod_rutils       !< module with a few helpful subroutines
    use amrex_base_module
    
    implicit none

    !>  Local variables
    type(pf_pfasst_t)     :: pf             !<  the main pfasst structure
    type(pf_comm_t)       :: comm           !<  the communicator (here it is mpi)
    type(pf_amrex_mfab_t) :: y_0            !<  the initial condition
    type(pf_amrex_mfab_t) :: y_end          !<  the solution at the final time
    character(256)        :: pf_fname       !<  file name for input parameters
    class(my_sweeper_t), pointer :: sweeper
    
    !> AMReX level specifics
    type(amrex_box)                             :: domain
    type(amrex_boxarray), allocatable,  target  :: ba(:)
    type(amrex_distromap), allocatable, target  :: dm(:)
    type(amrex_geometry), allocatable,  target  :: geom(:)      ! need one per level, since comtains domain & dx 
    
    !> local parameters
    integer                :: l                    ! loop variable over levels
    integer, allocatable   :: grid_size_amrex(:)   ! extended to account for amrex need of nghost (nx,ny,nz,ncomp,nghost)
    integer                :: nghost
    integer                :: box_per_proc         ! number of boxes handled per proc - needed for mpibuflen
    integer                :: mpibuflen            ! length of MPI buffer
    integer                :: ref_ratio
    type(amrex_string)     :: am_string
    
    !> mpi parameters
    integer :: space_comm, time_comm, space_color, time_color, err
    integer :: tRank, tSize, sRank, sSize, io
    logical :: group_space
    real    :: startT, endT
    
    !> read problem parameters
    call probin_init(pf_fname)
    
    !>  set up communicator
    group_space = .FALSE.
    call create_space_time_communicators(nspace, ntime, space_comm, time_comm, space_color, time_color, group_space)

    !call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_mpi_create(comm, time_comm)
    if (gRank .eq. 0) print *,'main.f90:intializing amrex'    
    !call amrex_init(comm=MPI_COMM_WORLD, arg_parmparse=.false., proc_parmparse=add_parameters)
    call amrex_init(comm=space_comm, arg_parmparse=.false., proc_parmparse=add_parameters)
    if (gRank .eq. 0) print *,'main.f90:amrex is initialized'
    
    !>  create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)
    call mpi_comm_size(space_comm, sSize, err)
    call mpi_comm_rank(space_comm, sRank,  err)
    call mpi_comm_size(time_comm, tSize, err)
    call mpi_comm_rank(time_comm, tRank,  err)
    call mpi_barrier(MPI_COMM_WORLD, err)
    
    !> Define grid_size & amrex_grid_size
    allocate(grid_size_amrex(Ndim+2))       ! [nx(,ny,nz),ncomp,nghost]
    nghost = 2
    grid_size(Ndim+1) = 1                   ! = ncomp
    grid_size(1:Ndim) = 0                   ! just default initalization
    grid_size_amrex = [grid_size, nghost]   ! extended to account for amrex need of nghost (nx,ny,nz,ncomp,nghost)

    !> Define geometry parameters
    call amrex_geometry_set_coord_sys(0)
    call amrex_geometry_set_prob_domain((/0.0d0,0.0d0,0.0d0/), (/Lx,Ly,Lz/))
    call amrex_geometry_set_periodic ((/ .true.,.true.,.true./))
    
    !> Define AMReX level specifics
    allocate(ba(pf%nlevels))
    allocate(dm(pf%nlevels))
    allocate(geom(pf%nlevels))
    ! loop over levels - could also be merged with lower level-loop
    do l = 1, pf%nlevels
      !> Define domain
      if (Ndim == 1) then
         domain = amrex_box((/0,0,0/), (/nx(l)-1,0,0/))
         if (l .gt. 1) then
            if (modulo(nx(l),nx(l-1)) .ne. 0) call pf_stop(__FILE__,__LINE__,'no integer-type refine-ratio')
         end if
      else if (Ndim == 2) then
         domain = amrex_box((/0,0,0/), (/nx(l)-1,ny(l)-1,0/))
         if (l .gt. 1) then
            if (modulo(nx(l),nx(l-1)) .ne. 0) call pf_stop(__FILE__,__LINE__,'no integer-type refine-ratio')
            ref_ratio = nx(l) / nx(l-1)
            if (ny(l)/ny(l-1) .ne. ref_ratio) call pf_stop(__FILE__,__LINE__,'different refine-ratio in x- and y-dimension')
         end if 
      else if (Ndim == 3) then
         domain = amrex_box((/0,0,0/), (/nx(l)-1,ny(l)-1,nz(l)-1/))
         if (l .gt. 1) then
            if (modulo(nx(l),nx(l-1)) .ne. 0) call pf_stop(__FILE__,__LINE__,'no integer-type refine-ratio')
            ref_ratio = nx(l) / nx(l-1)
            if (ny(l)/ny(l-1) .ne. ref_ratio) call pf_stop(__FILE__,__LINE__,'different refine-ratio in x- and y-dimension')
            if (nz(l)/nz(l-1) .ne. ref_ratio) call pf_stop(__FILE__,__LINE__,'different refine-ratio in x- and z-dimension')
         end if 
      end if

      !> Initialize the boxarray "ba" from the single box "domain"
      call amrex_boxarray_build(ba(l), domain)
      ! Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction 
      ! max_grid_size read from input-file - AMReX defaults: 2D -> 128, 3D -> 32
      call ba(l)%maxSize(max_grid_size)
      
      !> Build a DistributionMapping for the boxarray
      call amrex_distromap_build(dm(l), ba(l))

      !> This defines a amrex_geometry object.
      call amrex_geometry_build(geom(l), domain)

    end do
    
    !> loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)
      
       !>  Assign the factory for making the solution encapsulation
       allocate(pf_amrex_mfab_factory_t::pf%levels(l)%ulevel%factory)
       select type (fac => pf%levels(l)%ulevel%factory)
        type is (pf_amrex_mfab_factory_t)
          !>  Assign boxarray and distribution mapping to factory for AMReX mfabs
          fac%ba   => ba(l)
          fac%dm   => dm(l)
          fac%geom => geom(l)
       end select

       !>  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !>  Allocate the shape array for level 
       grid_size(1) = nx(l)                    ! modified to match definition of (nx,ny,nz,ncomp)
       if (Ndim > 1) grid_size(2) = ny(l)
       if (Ndim > 2) grid_size(3) = nz(l)
       grid_size_amrex(:) = [grid_size, nghost]     ! extended to account for amrex need of nghost (nx,ny,nz,ncomp,nghost)
       if (mod(ba(l)%nboxes(),sSize) .ne. 0) then
         if (gRank == 0) & 
               print *, "WARNING: best-performance if number of boxes is multiple of nspace!"
         if (gRank == 0) & 
               print *, "Number of boxes: ", ba(l)%nboxes(), " number of space-procs: ", sSize 
         box_per_proc = ba(l)%nboxes() / sSize + 1
       else
         box_per_proc = ba(l)%nboxes() / sSize
       end if 
       ! mpibuflen = length of data vector handled by processor 
       mpibuflen = box_per_proc * (max_grid_size ** Ndim) * grid_size(Ndim+1)  
       if (gRank .eq. 0) print *, 'Level: ', l
       if (gRank .eq. 0) print *, 'BoxArray: '
       if (gRank == 0) call amrex_print(ba(l))
       if (gRank .eq. 0) print *, 'Distribution Map: '
       if (gRank == 0) call amrex_print(dm(l))
       if (Ndim .eq. 1) then
        if (gRank .eq. 0) print *,'grid_size_amrex: ', grid_size_amrex, ' [nx,ncomp,nghost]'
       else if (Ndim .eq. 2) then
        if (gRank .eq. 0) print *,'grid_size_amrex: ', grid_size_amrex, ' [nx,ny,ncomp,nghost]'
       else if (Ndim .eq. 3) then 
        if (gRank .eq. 0) print *,'grid_size_amrex: ', grid_size_amrex, ' [nx,ny,nz,ncomp,nghost]'
       end if 
       if (gRank .eq. 0) print *,'mpibuflen: ', mpibuflen, ' = number of DoFs handled per processor without ghost cells - required for sending/receiving'
       call pf_level_set_size(pf,l,grid_size_amrex,mpibuflen)
    end do
    
    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)
    
    !> add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
    call pf_add_hook(pf, -1, PF_POST_CONVERGENCE, echo_error)
    call pf_add_hook(pf, -1, PF_POST_BLOCK, write_plotfile)
    
    
    !>  output the run options 
    if (gRank .eq. 0) call pf_print_options(pf,un_opt=6)

    !>  output local parameters
    if (gRank .eq. 0) call print_loc_options(pf,un_opt=6)
    
    !> allocate initial and final solutions
    call amrex_mfab_build(y_0, grid_size_amrex, ba(pf%nlevels), dm(pf%nlevels), geom(pf%nlevels))
    call amrex_mfab_build(y_end, grid_size_amrex, ba(pf%nlevels), dm(pf%nlevels), geom(pf%nlevels))
       
    !> compute initial condition on finest level (owned by sweeper)
    sweeper => as_my_sweeper(pf%levels(pf%nlevels)%ulevel%sweeper)        
    call set_ic(sweeper,y_0)    
    
    !> do the PFASST stepping
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call cpu_time(startT)
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)
    call cpu_time(endT)
    print 1000, "RANK ", gRank, "Elasped time: ", endT - startT, " [s]"
    1000 format (A5,I2,A20,F10.5,A4)
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    
    !> print final
    call amrex_string_build(am_string, "conc")
    !call amrex_write_plotfile("dat/"//trim(pf%outdir)//"/plt_file_initial_sol", &
    !                          1, [y_0%mfab], [am_string], [y_0%geom], real(0.0,8), [1], [1])
    !call amrex_write_plotfile("dat/"//trim(pf%outdir)//"/plt_file_final_sol", &
    !                          1, [y_end%mfab], [am_string], [y_end%geom], real(Tfin,8), [1], [1])

    !>  wait for everyone to be done
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    
    !>  deallocate initial condition, final solution and AMReX stuff
    call amrex_mfab_destroy(y_0)
    call amrex_mfab_destroy(y_end)
    call pf_pfasst_destroy(pf)
    call amrex_distromap_destroy(dm)
    call amrex_geometry_destroy(geom)
    
    !>  deallocate pfasst structure
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    !if (gRank .eq. 0) call pf_pfasst_destroy(pf)
    !call pf_pfasst_destroy(pf)
    deallocate(grid_size_amrex)
    
    !>  close down amrex
    call amrex_finalize()
    if (gRank .eq. 0) print *,'main.f90:amrex finalized'    
    
  end subroutine run_pfasst

  !> Routine to set initial condition.
  subroutine set_ic(this,y_0)
    use pf_my_sweeper  !<  Local module for defining sweeper and function evaluations
    use pf_mod_rutils    

    class(my_sweeper_t), intent(inout) :: this
    type(pf_amrex_mfab_t), intent(inout) :: y_0

    !call y_0%setval(1.0_pfdp)
    call exact(0.0_pfdp,y_0)

  end subroutine set_ic
end program main
  



!! need to look into this
subroutine add_parameters () bind(c)
  use amrex_base_module
  implicit none

  type(amrex_parmparse) :: pp

!!$  ! prefix "amrex"
!!$  call amrex_parmparse_build(pp,"amrex")
!!$  call pp%add("fpe_trap_invalid", 1)   ! turn on NaN trapping, which is off by default.
!!$  call amrex_parmparse_destroy(pp)
!!$
!!$  ! anonymous prefix
!!$  call amrex_parmparse_build(pp)
!!$  call pp%add("an_int_scalar", 2)      ! integer scalar: an_int_scalar
!!$  call pp%add("a_bool_scalar", .true.) ! logical scalar: a_bool_scalar
!!$  call pp%addarr("a_real_array", [1._amrex_real, 2._amrex_real, 3._amrex_real]) ! real array: a_real_array
!!$  call amrex_parmparse_destroy(pp)
!!$
!!$  ! prefix "a_prefix"
!!$  call amrex_parmparse_build(pp, "a_prefix")
!!$  call pp%addarr("an_int_array", [2, 3, 4])      ! integer array: a_prefix.an_int_array
!!$  call pp%add("a_real_scalar", 3.14_amrex_real)  ! real scalar  : a_prefix.a_real_scalar
!!$  call pp%add("a_string", "vonNeumann")          ! character    : a_prefix.a_string
!!$  call amrex_parmparse_destroy(pp)
end subroutine add_parameters
  


