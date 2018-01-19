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
!>  This module contains the routines to create, setup, and destroy the main data structure in PFASST
!!  See pf_dtype.f90 for the type definition
module pf_mod_pfasst
  use pf_mod_dtype
  implicit none
contains


  !> Create a PFASST object
  subroutine pf_pfasst_create(pf, comm, nlevels, fname, nocmd)
    use pf_mod_hooks, only: PF_MAX_HOOK

    use pf_mod_options
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

end module pf_mod_pfasst
