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
    integer :: ierr
    integer :: l                     !!  Loop variable for levels
    if (present(nlevels)) pf%nlevels = nlevels

    pf%outdir = "dat/"

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
       pf%levels(l)%nsteps_rk = pf%nsteps_rk(l)
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
    if (mpibuflen <= 0) call pf_stop(__FILE__,__LINE__,'allocate fail',mpibuflen)

    nnodes = lev%nnodes
    if (nnodes <= 0) call pf_stop(__FILE__,__LINE__,'allocate fail',nnodes)    

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
    !>  Allocate and compute all the matrices
    allocate(lev%sdcmats,stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error sdcmats")
    call pf_init_sdcmats(pf,lev%sdcmats, nnodes,lev%nflags)
    lev%nodes = lev%sdcmats%qnodes

    !>  initialize sweeper
    lev%ulevel%sweeper%use_LUq=pf%use_LUq
    call lev%ulevel%sweeper%initialize(pf,level_index)

    if (pf%use_rk_stepper)  call lev%ulevel%stepper%initialize(pf,level_index)

    !> allocate solution and function arrays
    npieces = lev%ulevel%sweeper%npieces

    call lev%ulevel%factory%create_array(lev%Q, nnodes, lev%index,  lev%lev_shape)
    call lev%ulevel%factory%create_array(lev%I, nnodes-1, lev%index,  lev%lev_shape)

    
    call lev%ulevel%factory%create_array(lev%Fflt, nnodes*npieces, lev%index,  lev%lev_shape)

    do i = 1, nnodes*npieces
       call lev%Fflt(i)%setval(0.0_pfdp, 0)
    end do

    lev%F(1:nnodes,1:npieces) => lev%Fflt

    call lev%ulevel%factory%create_array(lev%R, nnodes-1, lev%index,  lev%lev_shape)

    !  Need space for old function values in im sweepers
    call lev%ulevel%factory%create_array(lev%pFflt, nnodes*npieces, lev%index, lev%lev_shape)
    lev%pF(1:nnodes,1:npieces) => lev%pFflt
    if (lev%index < pf%nlevels) then
       call lev%ulevel%factory%create_array(lev%pQ, nnodes, lev%index,  lev%lev_shape)
    end if

    call lev%ulevel%factory%create_single(lev%qend, lev%index,   lev%lev_shape)
    call lev%ulevel%factory%create_single(lev%q0, lev%index,   lev%lev_shape)
    call lev%ulevel%factory%create_single(lev%q0_delta, lev%index,   lev%lev_shape)

  end subroutine pf_level_setup


  !> Deallocate PFASST object
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf  !!  Main pfasst structure

    integer :: l

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

    integer                          :: npieces  !!  local copy of number of function pieces
    class(pf_level_t), pointer :: lev    !!  points to current level    
    lev => pf%levels(level_index)   !!  Assign level pointer

    if (.not. lev%allocated) return

    !> deallocate flat buffers for communcition
    deallocate(lev%send)
    deallocate(lev%recv)

    !> deallocate nodes, flags, and integration matrices
    deallocate(lev%nodes)
    deallocate(lev%nflags)

    call pf_destroy_sdcmats(lev%sdcmats)
    deallocate(lev%sdcmats)

    !> deallocate solution and function storage
    npieces = lev%ulevel%sweeper%npieces

    if ((lev%index < pf%nlevels) .and. allocated(lev%tauQ)) then
       call lev%ulevel%factory%destroy_array(lev%tauQ)
    end if

    call lev%ulevel%factory%destroy_array(lev%Q)
    call lev%ulevel%factory%destroy_array(lev%Fflt)
    call lev%ulevel%factory%destroy_array(lev%I)
    call lev%ulevel%factory%destroy_array(lev%R)
    call lev%ulevel%factory%destroy_array(lev%pFflt)
    if (lev%index < pf%nlevels) then
       call lev%ulevel%factory%destroy_array(lev%pQ)
    end if
    if (lev%interp_workspace_allocated   .eqv. .true.) then      
       call lev%ulevel%factory%destroy_array(lev%c_delta)
       call lev%ulevel%factory%destroy_array(lev%cf_delta)
       lev%interp_workspace_allocated =.false.
    endif
 

    !> destroy the sweeper 
    call lev%ulevel%sweeper%destroy(pf,level_index)

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
    integer :: niters, nlevels, qtype
    integer :: nsweeps(PF_MAXLEVS)
    integer :: nsweeps_pred(PF_MAXLEVS) 
    integer :: nnodes(PF_MAXLEVS)
    integer :: nsteps_rk(PF_MAXLEVS)

    real(pfdp) :: abs_res_tol, rel_res_tol
    logical    :: PFASST_pred, RK_pred, pipeline_pred
    integer    ::  nsweeps_burn, q0_style, taui0
    logical    ::  Vcycle,Finterp, use_LUq, use_Sform
    logical    :: debug, use_rk_stepper
    logical    :: save_residuals, save_errors
    integer    :: save_timings
    logical    :: use_no_left_q,use_composite_nodes,use_proper_nodes
    
    ! stuff for reading the command line
    integer, parameter :: un = 9
    integer            :: i, ios,stat
    character(len=128)  :: arg
    character(len=256) :: istring  ! stores command line argument
    character(len=1024) :: message  ! use for i/o error messages
    character(len=256) :: outdir

    
    !> define the namelist for reading
    namelist /pf_params/ niters, nlevels, qtype, nsweeps, nsweeps_pred, nnodes, nsteps_rk, abs_res_tol, rel_res_tol
    namelist /pf_params/ PFASST_pred, RK_pred, pipeline_pred, nsweeps_burn, q0_style, taui0
    namelist /pf_params/ Vcycle,Finterp, use_LUq, use_Sform, debug, save_timings,save_residuals, save_errors, use_rk_stepper
    namelist /pf_params/ use_no_left_q,use_composite_nodes,use_proper_nodes, outdir

    !> set local variables to pf_pfasst defaults
    nlevels      = pf%nlevels
    niters       = pf%niters
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
    outdir       = pf%outdir
    debug        = pf%debug
    save_residuals = pf%save_residuals
    save_errors = pf%save_errors
    save_timings = pf%save_timings


    nsteps_rk    = pf%nsteps_rk
    rk_pred      = pf%rk_pred
    use_rk_stepper= pf%use_rk_stepper
    
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
    pf%save_timings = save_timings
    pf%save_errors = save_errors

    pf%use_rk_stepper=use_rk_stepper
    pf%nsteps_rk    = nsteps_rk    
    pf%rk_pred      = rk_pred

    pf%use_no_left_q       = use_no_left_q
    pf%use_composite_nodes = use_composite_nodes
    pf%use_proper_nodes    = use_proper_nodes

    !>  Sanity check
    if (pf%nlevels < 1) then
       call pf_stop(__FILE__,__LINE__,'Bad specification for nlevels',pf%nlevels)
    endif

  end subroutine pf_read_opts

  !>  Subroutine to write out run parameters
  subroutine pf_print_options(pf, un_opt, show_mats_opt,json_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: un_opt
    logical,           intent(in   ), optional :: show_mats_opt
    logical,           intent(in   ), optional :: json_opt

    integer :: un = 6
    logical :: show_mats = .FALSE.
    logical :: dump_json = .TRUE.
    integer :: l, i,istat
    character(8)   :: date
    character(10)  :: time
    character(len = 128) :: fname  !!  output file name for residuals
    character(len = 128) :: datpath  !!  path to output files

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt
    write(un,*) '=================================================='
    write(un,*) 'PFASST Configuration'
    write(un,*) '--------------------'

    call date_and_time(date=date, time=time)
    write(un,*) 'date:        ', date
    write(un,*) 'time:        ', time

    write(un,*) 'double precision:   ', pfdp   ,'  bytes'
    write(un,*) 'quad precision:   ', pfqp   ,'  bytes'    

    write(un,*) 'nlevels:     ', pf%nlevels, '! number of pfasst levels'
    write(un,*) 'nprocs:      ', pf%comm%nproc, '! number of pfasst "time" processors'
    
    if (pf%comm%nproc == 1) then
       write(un,*) '            ', '             ', ' ! since 1 time proc is being used, this is a serial sdc run'
    else
       write(un,*) '            ', '             ', ' ! since >1 time procs are being used, this is a parallel pfasst run'
    end if
    write(un,*) 'niters:      ', pf%niters, '! maximum number of sdc/pfasst iterations'
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

    if (pf%use_proper_nodes)  write(un,*) 'Using proper node nesting'
    if (pf%use_composite_nodes)  write(un,*) 'Using composite node nesting'
    if (pf%use_no_left_q)  write(un,*) ' Skipping left end point in quadruture rule '        
    write(un,*) 'nnodes:      ', pf%levels(1:pf%nlevels)%nnodes, '! number of sdc nodes per level'
    write(un,*) 'mpibuflen:   ', pf%levels(1:pf%nlevels)%mpibuflen, '! size of data send between time steps'
    write(un,*) 'nsweeps:     ', pf%levels(1:pf%nlevels)%nsweeps, '! number of sdc sweeps performed per visit to each level'
    write(un,*) 'nsweeps_pred:     ', pf%levels(1:pf%nlevels)%nsweeps_pred, '! number of sdc sweeps in predictor'
    write(un,*) 'taui0:     ',   pf%taui0, '! cutoff for tau correction'
    write(un,*) 'abs_res_tol:', pf%abs_res_tol, '! absolute residual tolerance: '
    write(un,*) 'rel_res_tol:', pf%rel_res_tol, '! relative residual tolerance: '
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

    write(un,*) 'Output directory ', pf%outdir    
    write(un,*) ''

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
    if (present(json_opt)) dump_json=json_opt
    if (dump_json) then
       ! Create a json file of all the pfasst parameters
       istat= system('mkdir -p dat')
       if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in pf_print_options")       
       istat= system('mkdir -p dat/' // trim(pf%outdir))       
       if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in pf_print_options")
       datpath= 'dat/' // trim(pf%outdir) 
       fname=trim(datpath) // 'pfasst_params.json'

       open(unit=321, file=trim(fname), form='formatted')
       write(321,*) '{'
       write(321,"(A24,I15,A1)")  '"nproc" :',       pf%comm%nproc, ','
       write(321,"(A24,I15,A1)")  '"nlevels" :',     pf%nlevels, ','
       write(321,"(A24,I15,A1)")  '"niters" :',      pf%niters, ','
       write(321,"(A24,I15,A1)")  '"qtype" :',       pf%qtype, ','
       write(321,"(A24,I15,A1)")  '"q0_style" :',    pf%q0_style, ','
       write(321,"(A24,I15,A1)")  '"taui0" :',       pf%taui0, ','
       write(321,"(A24,I15,A1)")  '"nsweeps_burn" :',pf%nsweeps_burn, ','
       write(321,"(A24,A15,A1)")  '"nnodes" :',      adjustr(convert_int_array(pf%nnodes(1:pf%nlevels),pf%nlevels)), ','
       write(321,"(A24,A15,A1)")  '"nsweeps" :',     adjustr(convert_int_array(pf%nsweeps(1:pf%nlevels),pf%nlevels)), ','
       write(321,"(A24,A15,A1)")  '"nsweeps_pred" :',adjustr(convert_int_array(pf%nsweeps_pred(1:pf%nlevels),pf%nlevels)), ','
       write(321,"(A24,A15,A1)")  '"nsteps_rk" :',   adjustr(convert_int_array(pf%nsteps_rk(1:pf%nlevels),pf%nlevels)), ','
       write(321,"(A24,e15.6,A1)") '"abs_res_tol" :',pf%abs_res_tol, ','
       write(321,"(A24,e15.6,A1)") '"rel_res_tol" :',pf%abs_res_tol, ','
       
       write(321,"(A24,A15,A1)")  '"use_proper_nodes" :',   convert_logical(pf%use_proper_nodes), ','
       write(321,"(A24,A15,A1)")  '"use_composite_nodes" :',convert_logical(pf%use_composite_nodes), ','
       write(321,"(A24,A15,A1)")  '"use_no_left_q" :',      convert_logical(pf%use_no_left_q), ','
       write(321,"(A24,A15,A1)")  '"PFASST_pred" :',        convert_logical(pf%PFASST_pred), ','
       write(321,"(A24,A15,A1)")  '"pipeline_pred" :',      convert_logical(pf%pipeline_pred), ','
       write(321,"(A24,A15,A1)")  '"Vcycle" :',             convert_logical(pf%Vcycle), ','
       write(321,"(A24,A15,A1)")  '"sweep_at_conv" :',      convert_logical(pf%sweep_at_conv), ','
       write(321,"(A24,A15,A1)")  '"Finterp" :',            convert_logical(pf%Finterp), ','
       write(321,"(A24,A15,A1)")  '"use_LUq" :',            convert_logical(pf%use_LUq), ','
       write(321,"(A24,A15,A1)")  '"use_Sform" :',          convert_logical(pf%use_Sform), ','
       write(321,"(A24,A15,A1)")  '"use_rk_stepper" :',     convert_logical(pf%use_rk_stepper), ','
       write(321,"(A24,A15,A1)")  '"RK_pred" :',            convert_logical(pf%RK_pred), ','
       write(321,"(A24,A15,A1)")  '"save_residuals" :',     convert_logical(pf%save_residuals), ','
       write(321,"(A24,I15,A1)")  '"save_timings" :',       pf%save_timings, ','
       write(321,"(A24,A15,A1)")  '"save_errors" :',        convert_logical(pf%save_errors), ','    
       write(321,"(A24,A15)")  '"debug" :',                 convert_logical(pf%debug)
       write(321,*) '}'    
    end if
  end subroutine pf_print_options

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


  !>  Subroutine to write out run parameters
  subroutine pf_initialize_results(pf)
  
  type(pf_pfasst_t), intent(inout)           :: pf

  integer :: level_index
  ALLOCATE(pf%results(pf%nlevels))
  do level_index = 1,pf%nlevels
     call  initialize_results(pf%results(level_index),pf%state%nsteps, pf%niters, pf%comm%nproc, pf%nsweeps(level_index),pf%rank,level_index,pf%outdir,pf%save_residuals)
  end do
  end subroutine pf_initialize_results


  !>  Subroutine to write out run parameters
  subroutine pf_dump_results(pf)
    
    type(pf_pfasst_t), intent(inout)           :: pf
    
    integer :: level_index
    
    if (pf%save_residuals) then
       do level_index = 1,pf%nlevels
          call  dump_resids(pf%results(level_index))
       end do
    end if
    
    if (pf%save_errors) then
       do level_index = 1,pf%nlevels
          call  dump_errors(pf%results(level_index))
       end do
    end if
    
    if (pf%save_timings > 0) then
       call  dump_timings(pf%results(pf%nlevels),pf)
    end if
  

end subroutine pf_dump_results

!>  Subroutine to destroy the results
  subroutine pf_destroy_results(pf)
    
    type(pf_pfasst_t), intent(inout)           :: pf
    
    integer :: level_index
    
       do level_index = 1,pf%nlevels
          call  destroy_results(pf%results(level_index))
       end do
  

end subroutine pf_destroy_results
  
  
    
  
end module pf_mod_pfasst
