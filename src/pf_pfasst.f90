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

module pf_mod_pfasst
  use pf_mod_dtype
  implicit none
contains

  !
  ! Create a PFASST object
  !
  subroutine pf_pfasst_create(pf, comm, nlevels, fname, nocmd)
    use pf_mod_hooks, only: PF_MAX_HOOK

    use pf_mod_options
    type(pf_pfasst_t), intent(inout)           :: pf
    type(pf_comm_t),   intent(inout), target   :: comm
    integer,           intent(in   ), optional :: nlevels
    character(len=*),  intent(in   ), optional :: fname
    logical,           intent(in   ), optional :: nocmd

    logical :: read_cmd

    if (present(nlevels)) pf%nlevels = nlevels

    pf%outdir = ""

    ! gather some input from a file and command line
    read_cmd = .true.
    if (present(nocmd)) then
         if (nocmd) read_cmd = .false.
    end if
    if (present(fname)) then
       call pf_read_opts(pf, read_cmd, fname)
    else
       if (read_cmd) call pf_read_opts(pf, read_cmd)
    end if

    pf%comm => comm

    allocate(pf%levels(pf%nlevels))
    allocate(pf%hooks(pf%nlevels, PF_MAX_HOOK, PF_MAX_HOOKS))
    allocate(pf%nhooks(pf%nlevels, PF_MAX_HOOK))
    pf%nhooks = 0

    allocate(pf%state)
    pf%state%pstatus = 0
    pf%state%status  = 0
  end subroutine pf_pfasst_create

  !
  ! Setup PFASST object
  !
  subroutine pf_pfasst_setup(pf)
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout), target :: pf

    class(pf_level_t), pointer :: F, G
    integer                   :: l

    if (pf%rank < 0) then
       stop 'Invalid PF rank: did you call setup correctly?'
    end if

    do l = 1, pf%nlevels
       pf%levels(l)%level = l
       call pf_level_setup(pf, pf%levels(l))
    end do

    do l = pf%nlevels, 2, -1
       F => pf%levels(l); G => pf%levels(l-1)
       allocate(F%tmat(F%nnodes,G%nnodes))
       allocate(F%rmat(G%nnodes,F%nnodes))
       call pf_time_interpolation_matrix(F%nodes, F%nnodes, G%nodes, G%nnodes, F%tmat)
       call pf_time_interpolation_matrix(G%nodes, G%nnodes, F%nodes, F%nnodes, F%rmat)
    end do

  end subroutine pf_pfasst_setup

  !
  ! Setup (allocate) PFASST level
  !
  ! If the level is already setup, calling this again will allocate
  ! (or deallocate) tau appropriately.
  !
  subroutine pf_level_setup(pf, F)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(in   )         :: pf
    class(pf_level_t), intent(inout), target :: F

    integer :: nvars, nnodes, npieces

    !
    ! do some sanity checks
    !

    if (F%nvars <= 0) stop "ERROR: Invalid nvars/dofs (pf_pfasst.f90)."
    if (F%nnodes <= 0) stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (F%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."

    nvars  = F%nvars
    nnodes = F%nnodes

    F%residual = -1.0_pfdp

    !
    ! (re)allocate tau (may to need create/destroy tau dynamically
    !                   when doing AMR)
    !
    if ((F%level < pf%nlevels) .and. (.not. allocated(F%tau))) then
       call F%ulevel%factory%create1(F%tau, nnodes-1, F%level, SDC_KIND_INTEGRAL, nvars, F%shape)
    else if ((F%level >= pf%nlevels) .and. (allocated(F%tau))) then
       deallocate(F%tau)
    end if

    !
    ! (re)allocate tauQ (may to need create/destroy tau dynamically
    !                   when doing AMR)
    !
    if ((F%level < pf%nlevels) .and. (.not. allocated(F%tauQ))) then
       call F%ulevel%factory%create1(F%tauQ, nnodes-1, F%level, SDC_KIND_INTEGRAL, nvars, F%shape)
    else if ((F%level >= pf%nlevels) .and. (allocated(F%tauQ))) then
       deallocate(F%tauQ)
    end if

    !
    ! skip the rest if we're already allocated
    !
    if (F%allocated) return
    F%allocated = .true.

    !
    ! allocate flat buffers (q0, send, and recv)
    !
    allocate(F%q0(nvars))
    allocate(F%send(nvars))
    allocate(F%recv(nvars))

    !
    ! nodes, flags, and integration matrices
    !
    allocate(F%nodes(nnodes))
    allocate(F%nflags(nnodes))
    allocate(F%s0mat(nnodes-1,nnodes))
    allocate(F%qmat(nnodes-1,nnodes))

    if (btest(pf%qtype, 8)) then
       call pf_quadrature(pf%qtype, nnodes, pf%levels(1)%nnodes, &
            F%nodes, F%nflags, F%s0mat, F%qmat)
    else
       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            F%nodes, F%nflags, F%s0mat, F%qmat)
    end if

    call F%ulevel%sweeper%initialize(F)


    !
    ! encaps
    !
    npieces = F%ulevel%sweeper%npieces
    call F%ulevel%factory%create1(F%Q, nnodes, F%level, SDC_KIND_SOL_FEVAL, nvars, F%shape)
    call F%ulevel%factory%create1(F%Fflt, nnodes*npieces, F%level, SDC_KIND_FEVAL, nvars, F%shape)
    F%F(1:nnodes,1:npieces) => F%Fflt
    call F%ulevel%factory%create1(F%S, nnodes-1, F%level, SDC_KIND_INTEGRAL, nvars, F%shape)
    call F%ulevel%factory%create1(F%I, nnodes-1, F%level, SDC_KIND_INTEGRAL, nvars, F%shape)
    call F%ulevel%factory%create1(F%R, nnodes-1, F%level, SDC_KIND_INTEGRAL, nvars, F%shape)
    if (F%level < pf%nlevels) then
       if (F%Finterp) then
          call F%ulevel%factory%create1(F%pFflt, nnodes*npieces, F%level, SDC_KIND_FEVAL, nvars, F%shape)
          F%pF(1:nnodes,1:npieces) => F%pFflt
       end if
       call F%ulevel%factory%create1(F%pQ, nnodes, F%level, SDC_KIND_SOL_NO_FEVAL, nvars, F%shape)
    end if
    call F%ulevel%factory%create0(F%qend, F%level, SDC_KIND_FEVAL, nvars, F%shape)

  end subroutine pf_level_setup

  !
  ! Deallocate PFASST object
  !
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: l
    do l = 1, pf%nlevels
       call pf_level_destroy(pf%levels(l))
    end do
    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)
    deallocate(pf%state)
  end subroutine pf_pfasst_destroy


  !
  ! Deallocate PFASST level
  !
  subroutine pf_level_destroy(F)
    class(pf_level_t), intent(inout) :: F

    if (.not. F%allocated) return

    ! flat buffers
    deallocate(F%q0)
    deallocate(F%send)
    deallocate(F%recv)

    ! nodes, flags, and integration matrices
    deallocate(F%nodes)
    deallocate(F%nflags)
    deallocate(F%qmat)
    deallocate(F%s0mat)

    ! encaps
    deallocate(F%Q)
    deallocate(F%F)
    deallocate(F%S)
    deallocate(F%I)
    deallocate(F%R)
    if (allocated(F%pQ)) then
       if (F%Finterp) then
          deallocate(F%pF)
       end if
       deallocate(F%pQ)
    end if
    deallocate(F%qend)
    if (allocated(F%tau)) then
       deallocate(F%tau)
    end if
    if (allocated(F%tauQ)) then
       deallocate(F%tauQ)
    end if

    ! other
    if (allocated(F%shape)) then
       deallocate(F%shape)
    end if

    if (allocated(F%tmat)) then
       deallocate(F%tmat)
    end if

    if (allocated(F%rmat)) then
       deallocate(F%rmat)
    end if
  end subroutine pf_level_destroy

end module pf_mod_pfasst
