!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!
!>  Module containing the routines to create, setup, and destroy the main data structure in PFASST
!!  See pf_dtype.f90 for the type definition
module pf_mod_pfasst
  use pf_mod_dtype
  implicit none
contains


  !> Create a PFASST object
  subroutine pf_pfasst_create(pf, comm, nlevels, fname, nocmd)
    use pf_mod_hooks, only: PF_MAX_HOOK


    type(pf_pfasst_t), intent(inout)           :: pf        !< Main pfasst object
    type(pf_comm_t),   intent(inout), target   :: comm      !< Communicator
    integer,           intent(in   ), optional :: nlevels   !< number of pfasst levels
    character(len=*),  intent(in   ), optional :: fname     !< Input file for pfasst parameters
    logical,           intent(in   ), optional :: nocmd     !< Determines if command line variables are to be read

    logical :: read_cmd              !< Local version of nocmd

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

    !>  allocate level pointers
    allocate(pf%levels(pf%nlevels))
    
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
    type(pf_pfasst_t), intent(inout), target :: pf   !<  Main pfasst structure

    class(pf_level_t), pointer :: lev_fine, lev_coarse  !<  Pointers to level structures for brevity
    integer                   :: l                      !<  Level loop index

    if (pf%rank < 0) then
       stop 'Invalid PF rank: did you call setup correctly?'
    end if

    !>  loop over levels to set parameters
    do l = 1, pf%nlevels
       pf%levels(l)%index = l
       call pf_level_setup(pf, pf%levels(l))
    end do

    !l  Loop over levels setting interpolation and restriction matrices (in time)
    do l = pf%nlevels, 2, -1
       lev_fine => pf%levels(l); lev_coarse => pf%levels(l-1)
       allocate(lev_fine%tmat(lev_fine%nnodes,lev_coarse%nnodes))
       allocate(lev_fine%rmat(lev_coarse%nnodes,lev_fine%nnodes))
       call pf_time_interpolation_matrix(lev_fine%nodes, lev_fine%nnodes, lev_coarse%nodes, lev_coarse%nnodes, lev_fine%tmat)
       call pf_time_interpolation_matrix(lev_coarse%nodes, lev_coarse%nnodes, lev_fine%nodes, lev_fine%nnodes, lev_fine%rmat)
    end do

  end subroutine pf_pfasst_setup

  !
  !> Setup (allocate) PFASST level
  !! If the level is already setup, calling this again will allocate
  !! (or deallocate) tauQ appropriately.
  subroutine pf_level_setup(pf, lev)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(in   )         :: pf   !<  Main pfasst structure
    class(pf_level_t), intent(inout), target :: lev  !<  Level to set up

    integer :: nvars, nnodes, npieces
    integer :: i

    !> do some sanity checks
    if (lev%nvars <= 0) stop "ERROR: Invalid nvars/dofs (pf_pfasst.f90)."
    if (lev%nnodes <= 0) stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (lev%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."

    nvars  = lev%nvars
    nnodes = lev%nnodes

    lev%residual = -1.0_pfdp


    !> (re)allocate tauQ (may to need create/destroy tauQ dynamically  when doing AMR)
    if ((lev%index < pf%nlevels) .and. (.not. allocated(lev%tauQ))) then
       call lev%ulevel%factory%create_array(lev%tauQ, nnodes-1, lev%index, SDC_KIND_INTEGRAL, nvars, lev%shape)
    else if ((lev%index >= pf%nlevels) .and. (allocated(lev%tauQ))) then
       deallocate(lev%tauQ)
    end if

    !> skip the rest if we're already allocated
    if (lev%allocated) return
    lev%allocated = .true.

    !> allocate flat buffers for send, and recv
    allocate(lev%send(nvars))
    allocate(lev%recv(nvars))


    !> allocate nodes, flags, and integration matrices
    allocate(lev%nodes(nnodes))
    allocate(lev%nflags(nnodes))
    allocate(lev%s0mat(nnodes-1,nnodes))
    allocate(lev%qmat(nnodes-1,nnodes))
    allocate(lev%qmatFE(nnodes-1,nnodes))
    allocate(lev%qmatBE(nnodes-1,nnodes))
    allocate(lev%LUmat(nnodes-1,nnodes))

    !> make quadrature matrices
    if (btest(pf%qtype, 8)) then
       call pf_quadrature(pf%qtype, nnodes, pf%levels(1)%nnodes, &
            lev%nodes, lev%nflags, lev%s0mat, lev%qmat,lev%qmatFE,lev%qmatBE)
    else
       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            lev%nodes, lev%nflags, lev%s0mat, lev%qmat,lev%qmatFE,lev%qmatBE)
    end if

    !>  initialize sweeper
    call lev%ulevel%sweeper%initialize(lev)



    !> allocate solution and function arrays
    npieces = lev%ulevel%sweeper%npieces

    call lev%ulevel%factory%create_array(lev%Q, nnodes, lev%index, SDC_KIND_SOL_FEVAL, nvars, lev%shape)
    call lev%ulevel%factory%create_array(lev%Fflt, nnodes*npieces, lev%index, SDC_KIND_FEVAL, nvars, lev%shape)
    do i = 1, nnodes*npieces
       call lev%Fflt(i)%setval(0.0_pfdp)
    end do
    lev%F(1:nnodes,1:npieces) => lev%Fflt
    call lev%ulevel%factory%create_array(lev%I, nnodes-1, lev%index, SDC_KIND_INTEGRAL, nvars, lev%shape)
    call lev%ulevel%factory%create_array(lev%R, nnodes-1, lev%index, SDC_KIND_INTEGRAL, nvars, lev%shape)

    if (lev%index < pf%nlevels) then
      if (lev%Finterp) then
          call lev%ulevel%factory%create_array(lev%pFflt, nnodes*npieces, lev%index, SDC_KIND_FEVAL, nvars, lev%shape)
          lev%pF(1:nnodes,1:npieces) => lev%pFflt
       end if
       call lev%ulevel%factory%create_array(lev%pQ, nnodes, lev%index, SDC_KIND_SOL_NO_FEVAL, nvars, lev%shape)
    end if
    call lev%ulevel%factory%create_single(lev%qend, lev%index, SDC_KIND_FEVAL, nvars, lev%shape)
    call lev%ulevel%factory%create_single(lev%q0, lev%index, SDC_KIND_FEVAL, nvars, lev%shape)


  end subroutine pf_level_setup


  !> Deallocate PFASST object
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf  !<  Main pfasst structure

    integer :: l

    !>  destroy all levels
    do l = 1, pf%nlevels
       call pf_level_destroy(pf%levels(l),pf%nlevels)
    end do
    !>  deallocate pfasst pointer arrays
    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)
    deallocate(pf%state)
  end subroutine pf_pfasst_destroy



  !> Deallocate PFASST level
  subroutine pf_level_destroy(lev,nlevels)
    class(pf_level_t), intent(inout) :: lev      !<  level to destroy
    integer                          :: nlevels  !<  number of pfasst levels

    
    integer                          :: npieces  !<  local copy of number of function pieces

    if (.not. lev%allocated) return

    !> deallocate flat buffers for communcition
    deallocate(lev%send)
    deallocate(lev%recv)

    !> deallocate nodes, flags, and integration matrices
    deallocate(lev%nodes)
    deallocate(lev%nflags)
    deallocate(lev%qmat)
    deallocate(lev%qmatFE)
    deallocate(lev%qmatBE)
    deallocate(lev%s0mat)
    deallocate(lev%LUmat)

    !> deallocate solution and function storage
    npieces = lev%ulevel%sweeper%npieces

    if ((lev%index < nlevels) .and. allocated(lev%tauQ)) then
       call lev%ulevel%factory%destroy_array(lev%tauQ, lev%nnodes-1, lev%index, SDC_KIND_INTEGRAL, lev%nvars, lev%shape)
    end if

    call lev%ulevel%factory%destroy_array(lev%Q, lev%nnodes, lev%index, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%destroy_array(lev%Fflt, lev%nnodes*npieces, lev%index, SDC_KIND_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%destroy_array(lev%I, lev%nnodes-1, lev%index, SDC_KIND_INTEGRAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%destroy_array(lev%R, lev%nnodes-1, lev%index, SDC_KIND_INTEGRAL, lev%nvars, lev%shape)
    if (lev%index < nlevels) then
       if (lev%Finterp) then
          call lev%ulevel%factory%destroy_array(lev%pFflt, lev%nnodes*npieces, lev%index, SDC_KIND_FEVAL, lev%nvars, lev%shape)
       end if
       call lev%ulevel%factory%destroy_array(lev%pQ, lev%nnodes, lev%index, SDC_KIND_SOL_NO_FEVAL, lev%nvars, lev%shape)
    end if
    call lev%ulevel%factory%destroy_single(lev%qend, lev%index, SDC_KIND_FEVAL, lev%nvars, lev%shape)
    call lev%ulevel%factory%destroy_single(lev%q0, lev%index, SDC_KIND_FEVAL, lev%nvars, lev%shape)

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
    integer          :: niters, nlevels, qtype, taui0
    double precision :: abs_res_tol, rel_res_tol
    logical          :: pipeline_g , pfasst_pred, echo_timings, debug, Vcycle

    ! stuff for reading the command line
    integer, parameter :: un = 9
    integer            :: i, ios
    character(len=32)  :: arg
    character(len=255) :: istring  ! stores command line argument
    character(len=255) :: message  ! use for i/o error messages
    character(len=512) :: outdir

    !> define the namelist for reading
    namelist /pf_params/ niters, nlevels, qtype, abs_res_tol, rel_res_tol, debug
    namelist /pf_params/ pipeline_g, pfasst_pred, echo_timings, taui0, outdir, Vcycle

    !> set local variables to pf_pfasst defaults
    nlevels      = pf%nlevels
    niters       = pf%niters
    qtype        = pf%qtype
    abs_res_tol  = pf%abs_res_tol
    rel_res_tol  = pf%rel_res_tol
    pipeline_g   = pf%pipeline_g
    pfasst_pred  = pf%pfasst_pred
    echo_timings = pf%echo_timings
    taui0        = pf%taui0
    outdir       = pf%outdir
    debug        = pf%debug
    Vcycle       = pf%Vcycle

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
    pf%niters       = niters
    pf%nlevels      = nlevels
    pf%qtype        = qtype
    pf%abs_res_tol  = abs_res_tol
    pf%rel_res_tol  = rel_res_tol
    pf%pipeline_g   = pipeline_g
    pf%pfasst_pred  = pfasst_pred
    pf%echo_timings = echo_timings
    pf%taui0        = taui0
    pf%outdir       = outdir
    pf%debug        = debug
    pf%Vcycle       = Vcycle

    !>  Sanity check
    if (pf%nlevels < 1) then
       write(*,*) 'Bad specification for nlevels=', pf%nlevels
       stop
    endif
  end subroutine pf_read_opts

  !>  Subroutine to write out run parameters
  subroutine pf_print_options(pf, unitno, show_mats)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: unitno
    logical,           intent(in   ), optional :: show_mats

    integer :: un = 6
    integer :: l, i
    character(8)   :: date
    character(10)  :: time

    print *,'Print options', pf%rank
    if (pf%rank /= 0) return
    if (present(unitno)) un = unitno

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
    write(un,*) 'nnodes:      ', pf%levels(1:pf%nlevels)%nnodes, '! number of sdc nodes per level'
    write(un,*) 'nvars:       ', pf%levels(1:pf%nlevels)%nvars, '! number of degrees of freedom per level'
    write(un,*) 'nsweeps:     ', pf%levels(1:pf%nlevels)%nsweeps, '! number of sdc sweeps performed per visit to each level'
    write(un,*) 'nsweeps_pred:     ', pf%levels(1:pf%nlevels)%nsweeps_pred, '! number of sdc sweeps in predictor'
    write(un,*) 'taui0:     ',   pf%taui0, '! cutoff for tau correction'

    if (pf%Pipeline_G) then
       write(un,*) 'Predictor Pipelining is ON    '
    else
       write(un,*) 'Predictor Pipelining is OFF    '
    end if
    if (pf%PFASST_pred) then
       write(un,*) 'PFASST Predictor style  '
    else
       write(un,*) 'Serial Predictor style  '
    end if
    if (pf%debug) write(un,*) 'Debug mode is on '

    write(un,*) ''

    if (present(show_mats)) then
       if (show_mats) then
          do l = 1, pf%nlevels
             print *, "Level", l
             print *, "-----------------"
             print *, "  nodes"
             print *, pf%levels(l)%nodes
             print *, "  Q"
             do i = 1, pf%levels(l)%nnodes-1
                 print *, pf%levels(l)%qmat(i,:)
             end do
          end do
       end if
    end if


  end subroutine pf_print_options

  !> Subroutine to make the matrices for interpolation  between noodes
  subroutine pf_time_interpolation_matrix(f_nodes, f_nnodes, c_nodes, c_nnodes, tmat)
    integer,    intent(in)  :: f_nnodes  !>  number of nodes on fine level
    integer,    intent(in)  :: c_nnodes  !>  number of nodes on coarse  level
    real(pfdp), intent(in)  :: f_nodes(0:f_nnodes-1)  !>  quadrature nodes on fine  level
    real(pfdp), intent(in)  :: c_nodes(0:c_nnodes-1)  !>  quadrature nodes on coarse  level
    real(pfdp), intent(out) :: tmat(0:f_nnodes-1,0:c_nnodes-1)  !>  Interpolation matrix to compute
    
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

  ! Subroutine to spread initial condition for predictors.
  subroutine spreadq0(lev, t0)
    class(pf_level_t), intent(inout) :: lev  !<  Level on which to spread
    real(pfdp),       intent(in)    :: t0    !<  time at beginning of interval

    integer :: m, p

    call lev%Q(1)%copy(lev%q0)

    call lev%ulevel%sweeper%evaluate(lev, t0, 1)

    do m = 2, lev%nnodes
       call lev%Q(m)%copy(lev%Q(1))
       do p = 1, lev%ulevel%sweeper%npieces
         call lev%F(m,p)%copy(lev%F(1,p))
       end do
    end do
  end subroutine spreadq0

end module pf_mod_pfasst
