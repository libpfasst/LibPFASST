!!  High level routines for PFASST data type
!
! This file is part of LIBPFASST.
!
!>  Module containing the routines to create, setup, and destroy the main data structure in PFASST
module pf_mod_pfasst
  use pf_mod_dtype
  use pf_mod_comm_mpi
  use pf_mod_utils
  use pf_mod_results
  
  implicit none
contains

  !> Create a PFASST object
  subroutine pf_pfasst_create(pf, comm, nlevels, fname, nocmd)
    use pf_mod_hooks, only: PF_MAX_HOOK


    type(pf_pfasst_t), intent(inout)           :: pf        !! Main pfasst object
    type(pf_comm_t),   intent(inout), target   :: comm      !! Communicator
    integer,           intent(in   ), optional :: nlevels   !! number of pfasst levels
    character(len=*),  intent(in   ), optional :: fname     !! Input file for pfasst parameters
    logical,           intent(in   ), optional :: nocmd     !! Determines if command line variables are to be read

    logical :: read_cmd              !! Local version of nocmd
    integer :: ierr                  !! Record system call error
    integer :: l                     !!  Loop variable for levels
    integer :: system                !!  For opening directory
    character(len=5) :: dirname     ! used to append output directory    
    if (present(nlevels)) pf%nlevels = nlevels


    !> gather some input from a file and command line
    read_cmd = .true.
    if (present(nocmd)) then
         if (nocmd) read_cmd = .false.
    end if
    if (present(fname)) then      !!  fname  present,  read inputs from a file (and maybe command line)
       call pf_read_opts(pf, read_cmd, fname)
    else                           !!  fname not present, only call read_opts if we want command line read
       if (read_cmd) call pf_read_opts(pf, read_cmd)
    end if

    !>  set communicator
    pf%comm => comm

    !>  Set up the mpi communicator
    call pf_mpi_setup(pf%comm, pf,ierr) 
    if (ierr /=0 )  call pf_stop(__FILE__,__LINE__,"ERROR: mpi_setup failed")

    if (pf%rank < 0) then
       call pf_stop(__FILE__,__LINE__,&
            "Invalid PF rank: did you call setup correctly?")
    end if

    !>  allocate level pointers
    allocate(pf%levels(pf%nlevels),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error",pf%nlevels)
    !>  loop over levels to set parameters
    do l = 1, pf%nlevels
       pf%levels(l)%index = l
       pf%levels(l)%nsweeps = pf%nsweeps(l)
       pf%levels(l)%nsweeps_pred = pf%nsweeps_pred(l)
       pf%levels(l)%nnodes = pf%nnodes(l)
       pf%levels(l)%Finterp = pf%Finterp
    end do
    
    !>  allocate hooks
    allocate(pf%hooks(pf%nlevels, PF_MAX_HOOK, PF_MAX_HOOKS),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error hooks")
    allocate(pf%nhooks(pf%nlevels, PF_MAX_HOOK),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error nhooks")
    pf%nhooks = 0

    !>  allocate status
    allocate(pf%state,stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error state")
    pf%state%pstatus = 0
    pf%state%status  = 0

    ! Create the output directory if it is not there
    ierr= system('mkdir -p dat')
    if (ierr .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory dat")       

    !  Stick the number of processors on the end of the output directory 
    write (dirname, "(A1,I0.4)") 'P',pf%comm%nproc
    pf%outdir       = trim(pf%outdir)//trim(dirname)
    ierr= system('mkdir -p dat/' // trim(pf%outdir))
    if (ierr .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make base directory")    

  end subroutine pf_pfasst_create

  !> Helper routine to set the size and mpi buffer length for regular grids
  subroutine pf_level_set_size(pf,level_index,shape_in,buflen_in)
    type(pf_pfasst_t), intent(inout) :: pf   !!  Main pfasst structure
    integer, intent(in)  ::  level_index
    integer, intent(in)  ::  shape_in(:)
    integer, intent(in),optional  ::  buflen_in

    integer ::  buflen_local,ierr

    ! Allocate and set shape array for the level
    allocate(pf%levels(level_index)%lev_shape(SIZE(shape_in)),stat=ierr)
    pf%levels(level_index)%lev_shape = shape_in

    !  Set the size of mpi buffer
    buflen_local= product(shape_in)
    if (present(buflen_in))     buflen_local= buflen_in
    
    pf%levels(level_index)%mpibuflen = buflen_local

    ! 
  end subroutine pf_level_set_size
  

  !> Setup both the PFASST object and the comm object
  subroutine pf_pfasst_setup(pf)
    type(pf_pfasst_t), intent(inout), target :: pf   !!  Main pfasst structure

    class(pf_level_t), pointer :: f_lev, c_lev  !!  Pointers to level structures for brevity
    integer                   :: l                      !!  Level loop index
    integer                   :: ierr                   !!  error flag


    !>  loop over levels to set parameters
    do l = 1, pf%nlevels
       call pf_level_setup(pf, l)
    end do
    !>  set default finest level
    pf%state%finest_level=pf%nlevels
    !>  Loop over levels setting interpolation and restriction matrices (in time)
    do l = pf%nlevels, 2, -1
       f_lev => pf%levels(l); c_lev => pf%levels(l-1)
       allocate(f_lev%tmat(f_lev%nnodes,c_lev%nnodes),stat=ierr)
       if (ierr /= 0)  call pf_stop(__FILE__,__LINE__,"allocate fail",f_lev%nnodes)

       allocate(f_lev%rmat(c_lev%nnodes,f_lev%nnodes),stat=ierr)
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate fail",f_lev%nnodes)

       
       ! with the RK stepper, no need to interpolate and restrict in time
       ! we only copy the first node and last node betweem levels
       if (pf%use_rk_stepper .eqv. .true.) then
          f_lev%tmat = 0.0_pfdp
          f_lev%rmat = 0.0_pfdp

          f_lev%tmat(1,1) = 1.0_pfdp
          f_lev%tmat(f_lev%nnodes,c_lev%nnodes) = 1.0_pfdp

          f_lev%rmat(1,1) = 1.0_pfdp
          f_lev%rmat(c_lev%nnodes,f_lev%nnodes) = 1.0_pfdp
       else         ! else compute the interpolation matrix
          call pf_time_interpolation_matrix(f_lev%nodes, f_lev%nnodes, c_lev%nodes, c_lev%nnodes, f_lev%tmat)
          call pf_time_interpolation_matrix(c_lev%nodes, c_lev%nnodes, f_lev%nodes, f_lev%nnodes, f_lev%rmat)
       endif
    end do
    

  end subroutine pf_pfasst_setup

  !> Setup (allocate) PFASST level
  !! If the level is already setup, calling this again will allocate
  !! (or deallocate) tauQ appropriately.
  subroutine pf_level_setup(pf, level_index)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(inout),target :: pf   !!  Main pfasst structure
    integer,           intent(in)    :: level_index  !!  level to set up

    class(pf_level_t),  pointer :: lev  !!  Level to set up

    integer :: mpibuflen, nnodes, npieces, nnodes0
    integer :: i,ierr

    lev => pf%levels(level_index)   !!  Assign level pointer
    
    !> do some sanity checks
    mpibuflen  = lev%mpibuflen
    if (mpibuflen <= 0) call pf_stop(__FILE__,__LINE__,'bad value for mpibulen=',mpibuflen)

    nnodes = lev%nnodes
    if (nnodes <= 0) call pf_stop(__FILE__,__LINE__,'bad value for nnodes=',nnodes)    

    lev%residual = -1.0_pfdp

    !> (re)allocate tauQ 
    if ((lev%index < pf%nlevels) .and. (.not. allocated(lev%tauQ))) then
       call lev%ulevel%factory%create_array(lev%tauQ, nnodes-1, lev%index,  lev%lev_shape)
    end if

    !> skip the rest if we're already allocated
    if (lev%allocated) return
    lev%allocated = .true.

    !> allocate flat buffers for send, and recv
    allocate(lev%send(mpibuflen),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate fail")
    allocate(lev%recv(mpibuflen),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate fail")
    !> allocate nodes, flags, and integration matrices
    allocate(lev%nodes(nnodes),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate fail")
    allocate(lev%nflags(nnodes),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate fail")
    lev%nflags=0
    !>  initialize sweeper
    if (pf%use_sdc_sweeper) then
       !>  Allocate and compute all the matrices
       allocate(lev%sdcmats,stat=ierr)
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error sdcmats")
       call pf_init_sdcmats(pf,lev%sdcmats, nnodes,lev%nflags)
       lev%nodes = lev%sdcmats%qnodes
       
       lev%ulevel%sweeper%use_LUq=pf%use_LUq
       call lev%ulevel%sweeper%initialize(pf,level_index)
    end if
    if (pf%use_rk_stepper)  call lev%ulevel%stepper%initialize(pf,level_index)
    !>  Allocate space for solutions 
    call lev%ulevel%factory%create_array(lev%Q, nnodes, lev%index,  lev%lev_shape)

    !> allocate solution and function arrays for sdc sweepers
    if (pf%use_sdc_sweeper) then
       npieces = lev%ulevel%sweeper%npieces
       call lev%ulevel%factory%create_array(lev%I, nnodes-1, lev%index,  lev%lev_shape)

       !  Space for function values
       call lev%ulevel%factory%create_array(lev%Fflt, nnodes*npieces, lev%index,  lev%lev_shape)
       do i = 1, nnodes*npieces
          call lev%Fflt(i)%setval(0.0_pfdp, 0)
       end do
       lev%F(1:nnodes,1:npieces) => lev%Fflt

       !  Need space for old function values in im sweepers
       call lev%ulevel%factory%create_array(lev%pFflt, nnodes*npieces, lev%index, lev%lev_shape)
       lev%pF(1:nnodes,1:npieces) => lev%pFflt
       if (lev%index < pf%nlevels) then
          call lev%ulevel%factory%create_array(lev%pQ, nnodes, lev%index,  lev%lev_shape)
       end if
    end if

    !>  Allocate space for residual and things need for sdc and rk
    call lev%ulevel%factory%create_array(lev%R, nnodes-1, lev%index,  lev%lev_shape)
    
    call lev%ulevel%factory%create_single(lev%qend, lev%index,   lev%lev_shape)
    call lev%ulevel%factory%create_single(lev%q0, lev%index,   lev%lev_shape)
    call lev%ulevel%factory%create_single(lev%delta_q0, lev%index,   lev%lev_shape)

  end subroutine pf_level_setup


  !> Deallocate PFASST object
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf  !!  Main pfasst structure

    integer :: l

    !>   deallocate results data
    call destroy_results(pf%results)

    !>  destroy all levels
    do l = 1, pf%nlevels
       call pf_level_destroy(pf,l)
    end do
    
    !>  deallocate pfasst pointer arrays

    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)
    deallocate(pf%state)
    call pf_mpi_destroy(pf%comm)

  end subroutine pf_pfasst_destroy


  !> Deallocate PFASST level
  subroutine pf_level_destroy(pf,level_index)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(inout),target :: pf  !!  Main pfasst structure    
    integer, intent(in)              :: level_index

    class(pf_level_t), pointer :: lev    !!  points to current level    
    lev => pf%levels(level_index)   !!  Assign level pointer

    if (.not. lev%allocated) return

    !> deallocate flat buffers for communcition
    deallocate(lev%send)
    deallocate(lev%recv)

    !> deallocate nodes, flags, and integration matrices
    deallocate(lev%nodes)
    deallocate(lev%nflags)
    
    !> deallocate solution and function storage
    if ((lev%index < pf%nlevels) .and. allocated(lev%tauQ)) then
       call lev%ulevel%factory%destroy_array(lev%tauQ)
    end if

    if (pf%use_sdc_sweeper) then
       call pf_destroy_sdcmats(lev%sdcmats)
       deallocate(lev%sdcmats)
       call lev%ulevel%factory%destroy_array(lev%Fflt)
       call lev%ulevel%factory%destroy_array(lev%I)
       call lev%ulevel%factory%destroy_array(lev%pFflt)
       if (lev%index < pf%nlevels) then
          call lev%ulevel%factory%destroy_array(lev%pQ)
       end if
    end if
    
    if (lev%interp_workspace_allocated   .eqv. .true.) then      
       call lev%ulevel%factory%destroy_array(lev%c_delta)
       call lev%ulevel%factory%destroy_array(lev%cf_delta)
       lev%interp_workspace_allocated =.false.
    endif
    if (lev%restrict_workspace_allocated   .eqv. .true.) then
        call lev%ulevel%factory%destroy_array(lev%f_encap_array_c)
        lev%restrict_workspace_allocated =.false.
    endif
    call lev%ulevel%factory%destroy_array(lev%Q)
    call lev%ulevel%factory%destroy_array(lev%R)
    call lev%ulevel%factory%destroy_single(lev%qend)
    call lev%ulevel%factory%destroy_single(lev%q0)
    call lev%ulevel%factory%destroy_single(lev%delta_q0)

    !> destroy the sweeper
    if (pf%use_sdc_sweeper)  call lev%ulevel%sweeper%destroy(pf,level_index)


    !> deallocate misc. arrays
    if (allocated(lev%lev_shape)) then
       deallocate(lev%lev_shape)
    end if

    if (allocated(lev%tmat)) then
       deallocate(lev%tmat)
    end if

    if (allocated(lev%rmat)) then
       deallocate(lev%rmat)
    end if

  end subroutine pf_level_destroy

  !>  Subroutine to read pfasst options from file and command line
  subroutine pf_read_opts(pf, read_cmd, fname)
    type(pf_pfasst_t), intent(inout)           :: pf
    logical,           intent(in   )           :: read_cmd
    character(len=*),  intent(in   ), optional :: fname
    
    ! local versions of pfasst parameters
    integer :: niters, MINiters,nlevels, qtype
    integer :: nsweeps(PF_MAXLEVS)
    integer :: nsweeps_pred(PF_MAXLEVS) 
    integer :: nnodes(PF_MAXLEVS)

    real(pfdp) :: abs_res_tol, rel_res_tol
    logical    :: PFASST_pred, RK_pred, pipeline_pred
    integer    ::  nsweeps_burn, q0_style, taui0
    logical    ::  Vcycle,use_pysdc_V,Finterp, use_LUq, use_Sform
    logical    :: debug, use_rk_stepper, use_sdc_sweeper, sweep_at_conv
    logical    :: save_residuals,save_delta_q0, save_errors
    integer    :: save_timings, save_solutions
    logical    :: use_no_left_q,use_composite_nodes,use_proper_nodes
    
    ! stuff for reading the command line
    integer, parameter :: un = 9
    integer            :: i, ios,stat
    character(len=128)  :: arg
    character(len=256) :: istring   ! stores command line argument
    character(len=1024) :: message  ! use for i/o error messages
    character(len=256) :: outdir    ! base name for output directory


    
    !> define the namelist for reading
    namelist /pf_params/ niters,MINiters, nlevels, qtype, nsweeps, nsweeps_pred, nnodes, abs_res_tol, rel_res_tol
    namelist /pf_params/ PFASST_pred, RK_pred, pipeline_pred, nsweeps_burn, q0_style, taui0
    namelist /pf_params/ Vcycle,Finterp,  debug, save_timings,save_residuals,save_delta_q0, save_errors, save_solutions
    namelist /pf_params/ use_sdc_sweeper,sweep_at_conv,use_pysdc_V,use_LUq, use_Sform
    namelist /pf_params/ use_no_left_q,use_composite_nodes,use_proper_nodes, use_rk_stepper, outdir

    !> set local variables to pf_pfasst defaults
    nlevels      = pf%nlevels
    niters       = pf%niters
    MINiters     = pf%MINiters
    qtype        = pf%qtype
    nsweeps      = pf%nsweeps
    nsweeps_pred = pf%nsweeps_pred
    nnodes       = pf%nnodes

    abs_res_tol  = pf%abs_res_tol
    rel_res_tol  = pf%rel_res_tol
    pfasst_pred  = pf%pfasst_pred
    pipeline_pred= pf%pipeline_pred
    nsweeps_burn = pf%nsweeps_burn
    q0_style     = pf%q0_style
    Vcycle       = pf%Vcycle
    Finterp      = pf%Finterp
    use_LUq      = pf%use_LUq
    use_Sform    = pf%use_Sform
    taui0        = pf%taui0
    outdir       = 'outdir'
    debug        = pf%debug
    save_residuals = pf%save_residuals
    save_delta_q0 = pf%save_delta_q0
    save_errors = pf%save_errors
    save_solutions = pf%save_solutions
    save_timings = pf%save_timings


    rk_pred      = pf%rk_pred
    use_rk_stepper= pf%use_rk_stepper
    use_sdc_sweeper= pf%use_sdc_sweeper
    use_pysdc_V= pf%use_pysdc_V
    sweep_at_conv= pf%sweep_at_conv
    
    use_no_left_q      = pf%use_no_left_q
    use_composite_nodes= pf%use_composite_nodes
    use_proper_nodes   = pf%use_proper_nodes

    !> open the file "fname" and read the pfasst namelist
    if (present(fname))  then
       open(unit=un, file=fname, status='old', action='read')
       read(unit=un, nml=pf_params)
       close(unit=un)
    end if

    !> overwrite parameters defined on  command line
    if (read_cmd) then
       i = 0
       do
          call get_command_argument(i, arg,status=stat)
          if (len_trim(arg) == 0) exit
          if (i > 1) then
             istring="&pf_params " // trim(arg) // " /"
             read(istring, nml=pf_params, iostat=ios, iomsg=message) ! internal read of namelist
          end if
          i = i+1
       end do
    end if

    !> re-assign the pfasst internals
    pf%nlevels      = nlevels
    pf%niters       = niters
    pf%MINiters     = MINiters
    pf%qtype        = qtype
    pf%nsweeps      = nsweeps
    pf%nsweeps_pred = nsweeps_pred
    pf%nnodes       = nnodes
    pf%abs_res_tol  = abs_res_tol
    pf%rel_res_tol  = rel_res_tol

    pf%pfasst_pred  = pfasst_pred
    pf%pipeline_pred= pipeline_pred
    pf%nsweeps_burn = nsweeps_burn
    pf%q0_style     = q0_style
    pf%Vcycle       = Vcycle
    pf%Finterp      = Finterp
    pf%use_LUq      = use_LUq
    pf%use_Sform    = use_Sform
    pf%taui0        = taui0

    pf%outdir       = outdir
    pf%debug        = debug
    pf%save_residuals = save_residuals
    pf%save_delta_q0 = save_delta_q0
    pf%save_timings = save_timings
    pf%save_errors = save_errors
    pf%save_solutions = save_solutions

    pf%use_rk_stepper=use_rk_stepper
    pf%use_sdc_sweeper=use_sdc_sweeper
    pf%rk_pred      = rk_pred
    pf%use_pysdc_V= use_pysdc_V
    pf%sweep_at_conv= sweep_at_conv

    pf%use_no_left_q       = use_no_left_q
    pf%use_composite_nodes = use_composite_nodes
    pf%use_proper_nodes    = use_proper_nodes

    !>  Sanity check
    if (pf%nlevels < 1) then
       call pf_stop(__FILE__,__LINE__,'Bad specification for nlevels',pf%nlevels)
    endif

  end subroutine pf_read_opts

  !>  Subroutine to write out run parameters
  subroutine pf_print_options(pf, un_opt, show_mats_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: un_opt
    logical,           intent(in   ), optional :: show_mats_opt

    integer :: un = 6
    logical :: show_mats = .FALSE.
    integer :: l, i,istat,system
    character(8)   :: date
    character(10)  :: time
    character(len = 128) :: fname  !!  output file name for residuals
    character(len = 128) :: datpath  !!  path to output files

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt
    write(un,*) '================================================'
    write(un,*) '----------- LibPFASST Parameters ---------------'

    call date_and_time(date=date, time=time)
    write(un,*) 'date:        ', date
    write(un,*) 'time:        ', time
    if(pf%use_sdc_sweeper) then
       write(un,*) 'method:        ',' PFASST'
    else
       write(un,*) 'method:        ',' parareal'
    end if
    write(un,*) 'double precision:   ', pfdp   ,'  bytes'
    write(un,*) 'quad precision:   ', pfqp   ,'  bytes'    
    write(un,*) 'Output directory: ', trim(pf%outdir)    
    write(un,*) 'Nprocs:      ', pf%comm%nproc, '! number of pfasst "time" processors'
    write(un,*) 'Nlevels:     ', pf%nlevels, '! number of levels'
    if(pf%use_sdc_sweeper) then
       write(un,*) 'Niters:      ', pf%niters, '! maximum number of sdc/pfasst iterations'
       write(un,*) 'MINiters:    ', pf%MINiters, '! Minimum number of sdc/pfasst iterations'
    else
       write(un,*) 'Niters:      ', pf%niters, '! maximum number of parareal iterations'
       write(un,*) 'MINiters:    ', pf%MINiters, '! Minimum number of parareal iterations'
    end if
    if (pf%use_sdc_sweeper) then
       if (pf%comm%nproc == 1) then
          write(un,*) '            ', '             ', ' ! since 1 time proc is being used, this is a serial sdc run'
       else
          write(un,*) '            ', '             ', ' ! since >1 time procs are being used, this is a parallel pfasst run'
       end if
       write(un,*) 'Nnodes:      ', pf%levels(1:pf%nlevels)%nnodes, '! number of sdc nodes per level'
       write(un,*) 'Nsweeps:     ', pf%levels(1:pf%nlevels)%nsweeps, '! number of sdc sweeps performed per visit to each level'
       write(un,*) 'Nsweeps_pred:     ', pf%levels(1:pf%nlevels)%nsweeps_pred, '! number of sdc sweeps in predictor'
       select case(pf%qtype)
       case (SDC_GAUSS_LEGENDRE)
          write(un,*) 'qtype:',pf%qtype, '! Gauss Legendre nodes are used'
       case (SDC_GAUSS_LOBATTO)
          write(un,*) 'qtype:',pf%qtype,'! Gauss Lobatto nodes are used'
       case (SDC_GAUSS_RADAU)
          write(un,*) 'qtype:',pf%qtype,'! Gauss Radau nodes are used'
       case (SDC_CLENSHAW_CURTIS)
          write(un,*) 'qtype:',pf%qtype,'! Clenshaw Curtis nodes are used'
       case (SDC_UNIFORM)
          write(un,*) 'qtype:', pf%qtype,'! Uniform  nodes are used'
       case (SDC_CHEBYSHEV)
          write(un,*) 'qtype:', pf%qtype,'! Chebyshev  nodes are used'
       case DEFAULT
          call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',pf%qtype)
       end select
       if (pf%use_Luq) then
          write(un,*) 'Implicit matrix is LU  '
       else
          write(un,*) 'Implicit matrix is backward Euler  '
       end if
       if (pf%use_Sform) then
          write(un,*) 'The Smat form of stepping is being done'
       else
          write(un,*) 'The Qmat form of stepping is being done'       
       end if
       if (present(show_mats_opt)) show_mats=show_mats_opt
       if (show_mats) then
          do l = 1, pf%nlevels
             print *, "Level", l
             print *, "-----------------"
             print *, "  nodes"
             print *, pf%levels(l)%nodes
             print *, "  Q"
             do i = 1, pf%levels(l)%nnodes-1
                print *, pf%levels(l)%sdcmats%qmat(i,:)
             end do
          end do
       end if
       if (pf%use_proper_nodes)  write(un,*) 'Using proper node nesting'
       if (pf%use_composite_nodes)  write(un,*) 'Using composite node nesting'
       if (pf%use_no_left_q)  write(un,*) ' Skipping left end point in quadruture rule '        
       write(un,*) 'taui0:     ',   pf%taui0, '! cutoff for tau correction'
       write(un,*) 'abs_res_tol:', pf%abs_res_tol, '! absolute residual tolerance: '
       write(un,*) 'rel_res_tol:', pf%rel_res_tol, '! relative residual tolerance: '
       write(un,*) 'mpibuflen:   ', pf%levels(1:pf%nlevels)%mpibuflen, '! size of data send between time steps'
       if (pf%Vcycle) then
          write(un,*) 'V-cycling is on'
       else
          write(un,*) 'V-cycling is off, fine level is pipelining'
       end if
       
       if (pf%rk_pred) then
          write(un,*) 'Runge-Kutta used for predictor'
       else
          if (pf%pipeline_pred) then
             write(un,*) 'Predictor pipelining is ON    '
          else
             write(un,*) 'Predictor pipelining is OFF    '
          end if
          if (pf%PFASST_pred) then
             write(un,*) 'PFASST Predictor style  '
          else
             write(un,*) 'Serial Predictor style  '
          end if
       endif
       
       if (pf%debug) write(un,*) 'Debug mode is on '
    end if

  end subroutine pf_print_options

  !>  Subroutine to dump stats to disk
  subroutine pf_dump_stats(pf)
    type(pf_pfasst_t), intent(inout)           :: pf   
    
    integer :: l, i,istat,system,un
    character(8)   :: date
    character(10)  :: time
    character(len = 128) :: fname  !!  output file name for residuals
    character(len = 128) :: datpath  !!  path to output files

    
    istat= system('mkdir -p dat/' // trim(pf%outdir))       
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in pf_print_options")
    datpath= 'dat/' // trim(pf%outdir)     
    !  Save the statistics before returning
    if (pf%save_timings > 0)  call dump_timingsl(pf%results,pf)
    
    if (pf%save_json) call pf_dump_json(pf)
    
    call dump_results(pf%results)


  end subroutine pf_dump_stats
  
  !>  Subroutine to write out run parameters
  subroutine pf_dump_json(pf)
    type(pf_pfasst_t), intent(inout)           :: pf   
    
    integer :: l, i,istat,system,un
    character(8)   :: date
    character(10)  :: time
    character(len = 128) :: fname  !!  output file name for residuals
    character(len = 128) :: datpath  !!  path to output files
    
    if (pf%rank /= 0) return
    
    istat= system('mkdir -p dat/' // trim(pf%outdir))       
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in pf_print_options")
    datpath= 'dat/' // trim(pf%outdir) 
    fname=trim(datpath) // '/pfasst_params.json'
    un=321
    open(unit=un, file=trim(fname), form='formatted')
    write(un,*) '{'
    if(pf%use_sdc_sweeper) then 
       write(un,*) '      "method" :  "PFASST",'
    else
       write(un,*) '      "method" :  "parareal",'
       write(un,123)  '"nsteps_rk" :',adjustr(convert_int_array(pf%nsteps_rk(1:pf%nlevels),pf%nlevels)), ','
       write(un,123)  '"rk_order" :',adjustr(convert_int_array(pf%rk_order(1:pf%nlevels),pf%nlevels)), ','
       write(un,123)  '"rk_nstages" :',adjustr(convert_int_array(pf%rk_nstages(1:pf%nlevels),pf%nlevels)), ','
    end if
    write(un,123)  '"use_rk_stepper" :',     convert_logical(pf%use_rk_stepper), ','
    write(un,123)  '"use_sdc_sweeper" :',     convert_logical(pf%use_sdc_sweeper), ','
    write(un,122)  '"nproc" :',       pf%comm%nproc, ','
    write(un,122)  '"nsteps" :',      pf%state%nsteps, ','
    write(un,122)  '"nlevels" :',     pf%nlevels, ','
    write(un,122)  '"niters" :',      pf%niters, ','
    write(un,122)  '"miniters" :',    pf%miniters, ','
    write(un,123)  '"nnodes" :',      adjustr(convert_int_array(pf%nnodes(1:pf%nlevels),pf%nlevels)), ','
    write(un,122)  '"q0_style" :',    pf%q0_style, ','
    write(un,123)  '"nsweeps" :',     adjustr(convert_int_array(pf%nsweeps(1:pf%nlevels),pf%nlevels)), ','
    if(pf%use_sdc_sweeper) then 
       write(un,122)  '"qtype" :',       pf%qtype, ','
       write(un,123)  '"nsweeps_pred" :',adjustr(convert_int_array(pf%nsweeps_pred(1:pf%nlevels),pf%nlevels)), ','
       write(un,122)  '"nsweeps_burn" :',pf%nsweeps_burn, ','
       write(un,122)  '"taui0" :',       pf%taui0, ','
    end if
    
    write(un,124) '"abs_res_tol" :',pf%abs_res_tol, ','
    write(un,124) '"rel_res_tol" :',pf%abs_res_tol, ','
    if(pf%use_sdc_sweeper) then        
       write(un,123)  '"use_proper_nodes" :',   convert_logical(pf%use_proper_nodes), ','
       write(un,123)  '"use_composite_nodes" :',convert_logical(pf%use_composite_nodes), ','
       write(un,123)  '"use_no_left_q" :',      convert_logical(pf%use_no_left_q), ','
       write(un,123)  '"PFASST_pred" :',        convert_logical(pf%PFASST_pred), ','
       write(un,123)  '"pipeline_pred" :',      convert_logical(pf%pipeline_pred), ','
       write(un,123)  '"sweep_at_conv" :',      convert_logical(pf%sweep_at_conv), ','
       write(un,123)  '"use_LUq" :',            convert_logical(pf%use_LUq), ','
       write(un,123)  '"use_Sform" :',          convert_logical(pf%use_Sform), ','
       write(un,123)  '"Finterp" :',            convert_logical(pf%Finterp), ','
    end if
    write(un,123)  '"Vcycle" :',             convert_logical(pf%Vcycle), ','
    write(un,123)  '"RK_pred" :',            convert_logical(pf%RK_pred), ','
    write(un,123)  '"save_residuals" :',     convert_logical(pf%save_residuals), ','
    write(un,123)  '"save_delta_q0" :',      convert_logical(pf%save_delta_q0), ','
    write(un,122)  '"save_timings" :',       pf%save_timings, ','
    write(un,122)  '"save_solutions" :',     pf%save_solutions, ','
    write(un,123)  '"save_errors" :',        convert_logical(pf%save_errors), ','    
    write(un,123)  '"debug" :',                 convert_logical(pf%debug),','    
    write(un,"(A24,A64)")  '"outdir" : ',       adjustr('"'//trim(pf%outdir)//'"')
    write(un,*) '}'
122 FORMAT (A24,I15,A1)
123 FORMAT (A24,A15,A1)
124 FORMAT (A24,e15.6,A1)
    close(unit=un)       
  end subroutine pf_dump_json
  !> Subroutine to make the matrices for interpolation  between noodes
  subroutine pf_time_interpolation_matrix(f_nodes, f_nnodes, c_nodes, c_nnodes, tmat)
    integer,    intent(in)  :: f_nnodes  !!  number of nodes on fine level
    integer,    intent(in)  :: c_nnodes  !!  number of nodes on coarse  level
    real(pfdp), intent(in)  :: f_nodes(0:f_nnodes-1)  !!  quadrature nodes on fine  level
    real(pfdp), intent(in)  :: c_nodes(0:c_nnodes-1)  !!  quadrature nodes on coarse  level
    real(pfdp), intent(out) :: tmat(0:f_nnodes-1,0:c_nnodes-1)  !!  Interpolation matrix to compute
    
    integer    :: i, j, k
    real(pfqp) :: xi, num, den
    
    do i = 0, f_nnodes-1
       xi = real(f_nodes(i), pfqp)
       
       do j = 0, c_nnodes-1
          den = 1.0_pfqp
          num = 1.0_pfqp
          
          do k = 0, c_nnodes-1
             if (k == j) cycle
             den = den * real(c_nodes(j) - c_nodes(k),pfqp)
             num = num * real(xi        - c_nodes(k),pfqp)
          end do
          tmat(i, j) = real(num/den,pfdp)
       end do
    end do
  end subroutine pf_time_interpolation_matrix

end module pf_mod_pfasst
