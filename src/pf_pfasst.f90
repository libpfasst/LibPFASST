!! Main data structure routines 
!
! This file is part of LIBPFASST.
!
!>  Module containing the routines to create, setup, and destroy the main data structure in PFASST
!!  See pf_dtype.f90 for the type definition
module pf_mod_pfasst
  use pf_mod_dtype
  use pf_mod_comm_mpi
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
    integer :: ierror
    integer :: l                     !!  Loop variable for levels
    if (present(nlevels)) pf%nlevels = nlevels

    pf%outdir = ""

    !> gather some input from a file and command line
    read_cmd = .true.
    if (present(nocmd)) then
         if (nocmd) read_cmd = .false.
    end if
    if (present(fname)) then
       call pf_read_opts(pf, read_cmd, fname)
    else
       if (read_cmd) call pf_read_opts(pf, read_cmd)
    end if

    !>  set communicator
    pf%comm => comm

    !>  Set up the mpi communicator
    call pf_mpi_setup(pf%comm, pf,ierror) 
    if (ierror /=0 )        stop "ERROR: mpi_setup failed"
    

    if (pf%rank < 0) then
       stop 'Invalid PF rank: did you call setup correctly?'
    end if

    !>  allocate level pointers
    allocate(pf%levels(pf%nlevels))

    !>  loop over levels to set parameters
    do l = 1, pf%nlevels
       pf%levels(l)%index = l
       pf%levels(l)%nsweeps = pf%nsweeps(l)
       pf%levels(l)%nsweeps_pred = pf%nsweeps_pred(l)
       pf%levels(l)%nnodes = pf%nnodes(l)              
    end do
    
    !>  allocate hooks
    allocate(pf%hooks(pf%nlevels, PF_MAX_HOOK, PF_MAX_HOOKS))
    allocate(pf%nhooks(pf%nlevels, PF_MAX_HOOK))
    pf%nhooks = 0

    !>  allocate status
    allocate(pf%state)
    pf%state%pstatus = 0
    pf%state%status  = 0


  end subroutine pf_pfasst_create


  !> Setup both the PFASST object and the comm object
  subroutine pf_pfasst_setup(pf)
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout), target :: pf   !!  Main pfasst structure

    class(pf_level_t), pointer :: lev_fine, lev_coarse  !!  Pointers to level structures for brevity
    integer                   :: l                      !!  Level loop index

   !>  loop over levels to set parameters
    do l = 1, pf%nlevels
       call pf_level_setup(pf, pf%levels(l))
    end do

    !>  Loop over levels setting interpolation and restriction matrices (in time)
    do l = pf%nlevels, 2, -1
       lev_fine => pf%levels(l); lev_coarse => pf%levels(l-1)
       allocate(lev_fine%tmat(lev_fine%nnodes,lev_coarse%nnodes))
       allocate(lev_fine%rmat(lev_coarse%nnodes,lev_fine%nnodes))
       ! with the RK stepper, no need to interpolate and restrict in time
       ! we only copy the first node and last node betweem levels
       if (pf%use_rk_stepper .eqv. .true.) then
          lev_fine%tmat = 0.0_pfdp
          lev_fine%rmat = 0.0_pfdp

          lev_fine%tmat(1,1) = 1.0_pfdp
          lev_fine%tmat(lev_fine%nnodes,lev_coarse%nnodes) = 1.0_pfdp

          lev_fine%rmat(1,1) = 1.0_pfdp
          lev_fine%rmat(lev_coarse%nnodes,lev_fine%nnodes) = 1.0_pfdp
       else         ! else compute the interpolation matrix
          call pf_time_interpolation_matrix(lev_fine%nodes, lev_fine%nnodes, lev_coarse%nodes, lev_coarse%nnodes, lev_fine%tmat)
          call pf_time_interpolation_matrix(lev_coarse%nodes, lev_coarse%nnodes, lev_fine%nodes, lev_fine%nnodes, lev_fine%rmat)
       endif
    end do

  end subroutine pf_pfasst_setup

  !
  !> Setup (allocate) PFASST level
  !! If the level is already setup, calling this again will allocate
  !! (or deallocate) tauQ appropriately.
  subroutine pf_level_setup(pf, lev)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(in   )         :: pf   !!  Main pfasst structure
    class(pf_level_t), intent(inout), target :: lev  !!  Level to set up

    integer :: mpibuflen, nnodes, npieces, nnodes0
    integer :: i

    !> do some sanity checks
    if (lev%mpibuflen <= 0) stop "ERROR: Invalid mpibuflen/dofs (pf_pfasst.f90)."
    if (lev%nnodes <= 0) stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (lev%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."

    mpibuflen  = lev%mpibuflen
    nnodes = lev%nnodes

    lev%residual = -1.0_pfdp


    !> (re)allocate tauQ (may to need create/destroy tauQ dynamically  when doing AMR)
    if ((lev%index < pf%nlevels) .and. (.not. allocated(lev%tauQ))) then
       call lev%ulevel%factory%create_array(lev%tauQ, nnodes-1, lev%index,  lev%shape)
    else if ((lev%index >= pf%nlevels) .and. (allocated(lev%tauQ))) then
       deallocate(lev%tauQ)
    end if

    !> skip the rest if we're already allocated
    if (lev%allocated) return
    lev%allocated = .true.

    !> allocate flat buffers for send, and recv
    allocate(lev%send(mpibuflen))
    allocate(lev%recv(mpibuflen))


    !> allocate nodes, flags, and integration matrices
    allocate(lev%nodes(nnodes))
    allocate(lev%nflags(nnodes))
    
    !> make quadrature matrices
    if (btest(pf%qtype, 8)) then
       nnodes0=pf%levels(1)%nnodes
    else
       nnodes0=pf%levels(pf%nlevels)%nnodes
    end if
    !>  Allocate and compute all the matrices
    call pf_init_sdcmats(lev%sdcmats,pf%qtype, nnodes,nnodes0,lev%nflags)
    lev%nodes = lev%sdcmats%qnodes
    
    !>  initialize sweeper
    lev%ulevel%sweeper%use_LUq=pf%use_LUq
    call lev%ulevel%sweeper%initialize(lev)

    
    if (pf%use_rk_stepper)  call lev%ulevel%stepper%initialize(lev)

    !> allocate solution and function arrays
    npieces = lev%ulevel%sweeper%npieces

    call lev%ulevel%factory%create_array(lev%Q, nnodes, lev%index,  lev%shape)
    call lev%ulevel%factory%create_array(lev%Fflt, nnodes*npieces, lev%index,  lev%shape)
    do i = 1, nnodes*npieces
       call lev%Fflt(i)%setval(0.0_pfdp, 0)
    end do
    lev%F(1:nnodes,1:npieces) => lev%Fflt
    call lev%ulevel%factory%create_array(lev%I, nnodes-1, lev%index,  lev%shape)
    call lev%ulevel%factory%create_array(lev%R, nnodes-1, lev%index,  lev%shape)

    !  Need space for old function values in imexR sweepers
    call lev%ulevel%factory%create_array(lev%pFflt, nnodes*npieces, lev%index, lev%shape)
    lev%pF(1:nnodes,1:npieces) => lev%pFflt
    if (lev%index < pf%nlevels) then
       call lev%ulevel%factory%create_array(lev%pQ, nnodes, lev%index,  lev%shape)
    end if
    call lev%ulevel%factory%create_single(lev%qend, lev%index,   lev%shape)
    call lev%ulevel%factory%create_single(lev%q0, lev%index,   lev%shape)

    
  end subroutine pf_level_setup


  !> Deallocate PFASST object
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf  !!  Main pfasst structure

    integer :: l

    !>  destroy all levels
    do l = 1, pf%nlevels
       call pf_level_destroy(pf%levels(l),pf%nlevels)
    end do
    !>  deallocate pfasst pointer arrays
    call  pf%results%destroy(pf%results)
    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)
    deallocate(pf%state)
    call pf_mpi_destroy(pf%comm)

  end subroutine pf_pfasst_destroy


  !> Deallocate PFASST level
  subroutine pf_level_destroy(lev,nlevels)
    use pf_mod_quadrature
    class(pf_level_t), intent(inout) :: lev      !!  level to destroy
    integer                          :: nlevels  !!  number of pfasst levels


    integer                          :: npieces  !!  local copy of number of function pieces

    if (.not. lev%allocated) return

    !> deallocate flat buffers for communcition
    deallocate(lev%send)
    deallocate(lev%recv)

    !> deallocate nodes, flags, and integration matrices
    deallocate(lev%nodes)
    deallocate(lev%nflags)

    call pf_destroy_sdcmats(lev%sdcmats)
    !> deallocate solution and function storage
    npieces = lev%ulevel%sweeper%npieces

    if ((lev%index < nlevels) .and. allocated(lev%tauQ)) then
       call lev%ulevel%factory%destroy_array(lev%tauQ, lev%nnodes-1, lev%index,   lev%shape)
    end if

    call lev%ulevel%factory%destroy_array(lev%Q, lev%nnodes, lev%index,   lev%shape)
    call lev%ulevel%factory%destroy_array(lev%Fflt, lev%nnodes*npieces, lev%index,   lev%shape)
    call lev%ulevel%factory%destroy_array(lev%I, lev%nnodes-1, lev%index,  lev%shape)
    call lev%ulevel%factory%destroy_array(lev%R, lev%nnodes-1, lev%index,  lev%shape)
    call lev%ulevel%factory%destroy_array(lev%pFflt, lev%nnodes*npieces, lev%index, lev%shape)
    if (lev%index < nlevels) then
       call lev%ulevel%factory%destroy_array(lev%pQ, lev%nnodes, lev%index,   lev%shape)
    end if
    call lev%ulevel%factory%destroy_single(lev%qend, lev%index,  lev%shape)
    call lev%ulevel%factory%destroy_single(lev%q0, lev%index,   lev%shape)

    !> destroy the sweeper 
    call lev%ulevel%sweeper%destroy(lev)

    !> deallocate misc. arrays
    if (allocated(lev%shape)) then
       deallocate(lev%shape)
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
    integer :: nnodes_rk(PF_MAXLEVS)

    real(pfdp) :: abs_res_tol, rel_res_tol
    logical    :: PFASST_pred, RK_pred, pipeline_pred
    integer    ::  nsweeps_burn, q0_style, taui0
    logical    ::  Vcycle,Finterp, use_LUq
    logical    :: echo_timings, debug, save_results, use_rk_stepper
    
    ! stuff for reading the command line
    integer, parameter :: un = 9
    integer            :: i, ios
    character(len=32)  :: arg
    character(len=255) :: istring  ! stores command line argument
    character(len=255) :: message  ! use for i/o error messages
    character(len=512) :: outdir

    
    !> define the namelist for reading
    namelist /pf_params/ niters, nlevels, qtype, nsweeps, nsweeps_pred, nnodes, nnodes_rk, abs_res_tol, rel_res_tol
    namelist /pf_params/ PFASST_pred, RK_pred, pipeline_pred, nsweeps_burn, q0_style, taui0
    namelist /pf_params/ Vcycle,Finterp, use_LUq, echo_timings, debug, save_results, use_rk_stepper


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
    taui0        = pf%taui0
    outdir       = pf%outdir
    debug        = pf%debug
    save_results = pf%save_results
    echo_timings = pf%echo_timings

    nnodes_rk    = pf%nnodes_rk
    rk_pred      = pf%rk_pred
    use_rk_stepper= pf%use_rk_stepper

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
          call get_command_argument(i, arg)
          if (len_trim(arg) == 0) exit
          if (i > 0) then
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
    pf%taui0        = taui0

    pf%echo_timings = echo_timings
    pf%outdir       = outdir
    pf%debug        = debug
    pf%save_results = save_results
    pf%echo_timings = echo_timings

    pf%use_rk_stepper=use_rk_stepper
    pf%nnodes_rk    = nnodes_rk    
    pf%rk_pred      = rk_pred

    !>  Sanity check
    if (pf%nlevels < 1) then
       write(*,*) 'Bad specification for nlevels=', pf%nlevels
       stop
    endif
  end subroutine pf_read_opts

  !>  Subroutine to write out run parameters
  subroutine pf_print_options(pf, un_opt, show_mats_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: un_opt
    logical,           intent(in   ), optional :: show_mats_opt

    integer :: un = 6
    logical :: show_mats = .FALSE.
    integer :: l, i
    character(8)   :: date
    character(10)  :: time

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt
    
    write(un,*) 'PFASST Configuration'
    write(un,*) '===================='

    call date_and_time(date=date, time=time)
    write(un,*) 'date:        ', date
    write(un,*) 'time:        ', time

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
       write(un,*) 'qtype:',pf%qtype,'! Gauss Radua nodes are used'
    case (SDC_CLENSHAW_CURTIS)
       write(un,*) 'qtype:',pf%qtype,'! Clenshaw Curtis nodes are used'
    case (SDC_UNIFORM)
       write(un,*) 'qtype:', pf%qtype,'! Uniform  nodes are used'
    case default
       print *,'qtype = ',pf%qtype
       stop "ERROR: Invalid qtype"
    end select

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
 


  end subroutine pf_print_options

  !> Subroutine to make the matrices for interpolation  between noodes
  subroutine pf_time_interpolation_matrix(f_nodes, f_nnodes, c_nodes, c_nnodes, tmat)
    integer,    intent(in)  :: f_nnodes  !!  number of nodes on fine level
    integer,    intent(in)  :: c_nnodes  !!  number of nodes on coarse  level
    real(pfdp), intent(in)  :: f_nodes(0:f_nnodes-1)  !!  quadrature nodes on fine  level
    real(pfdp), intent(in)  :: c_nodes(0:c_nnodes-1)  !!  quadrature nodes on coarse  level
    real(pfdp), intent(out) :: tmat(0:f_nnodes-1,0:c_nnodes-1)  !!  Interpolation matrix to compute
    
    integer    :: i, j, k
    real(pfdp) :: xi, num, den
    
    do i = 0, f_nnodes-1
       xi = f_nodes(i)
       
       do j = 0, c_nnodes-1
          den = 1.0_pfdp
          num = 1.0_pfdp
          
          do k = 0, c_nnodes-1
             if (k == j) cycle
             den = den * (c_nodes(j) - c_nodes(k))
             num = num * (xi        - c_nodes(k))
          end do
          
          tmat(i, j) = num/den
       end do
    end do
  end subroutine pf_time_interpolation_matrix


end module pf_mod_pfasst
