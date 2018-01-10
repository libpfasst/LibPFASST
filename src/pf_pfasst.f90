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
  !< Setup both the PFASST object and the comm object
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
       pf%levels(l)%index = l
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
  ! (or deallocate) tauQ appropriately.
  !
  subroutine pf_level_setup(pf, F)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(in   )         :: pf
    class(pf_level_t), intent(inout), target :: F

    integer :: nvars, nnodes, npieces
    integer :: i
    !
    ! do some sanity checks
    !

    if (F%nvars <= 0) stop "ERROR: Invalid nvars/dofs (pf_pfasst.f90)."
    if (F%nnodes <= 0) stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (F%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."

    print*, 'nnodes = ', F%nnodes
    nvars  = F%nvars
    nnodes = F%nnodes

    F%residual = -1.0_pfdp

    !
    ! (re)allocate tauQ (may to need create/destroy tauQ dynamically
    !                   when doing AMR)
    !
    if ((F%index < pf%nlevels) .and. (.not. allocated(F%tauQ))) then
       call F%ulevel%factory%create_array(F%tauQ, nnodes-1, F%index, SDC_KIND_INTEGRAL, nvars, F%shape)
    else if ((F%index >= pf%nlevels) .and. (allocated(F%tauQ))) then
       deallocate(F%tauQ)
    end if

    !
    ! skip the rest if we're already allocated
    !
    if (F%allocated) return
    F%allocated = .true.

    !
    ! allocate flat buffers (send, and recv)
    !
    allocate(F%send(nvars))
    allocate(F%recv(nvars))

    !
    ! nodes, flags, and integration matrices
    !
    allocate(F%nodes(nnodes))
    allocate(F%nflags(nnodes))
    allocate(F%s0mat(nnodes-1,nnodes))
    allocate(F%qmat(nnodes-1,nnodes))
    allocate(F%qmatFE(nnodes-1,nnodes))
    allocate(F%qmatBE(nnodes-1,nnodes))
    allocate(F%LUmat(nnodes-1,nnodes))

    if (btest(pf%qtype, 8)) then
       call pf_quadrature(pf%qtype, nnodes, pf%levels(1)%nnodes, &
            F%nodes, F%nflags, F%s0mat, F%qmat,F%qmatFE,F%qmatBE)
    else
       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            F%nodes, F%nflags, F%s0mat, F%qmat,F%qmatFE,F%qmatBE)
    end if

    call F%ulevel%sweeper%initialize(F)


    !
    ! encaps
    !
    npieces = F%ulevel%sweeper%npieces

    call F%ulevel%factory%create_array(F%Q, nnodes, F%index, SDC_KIND_SOL_FEVAL, nvars, F%shape)
    call F%ulevel%factory%create_array(F%Fflt, nnodes*npieces, F%index, SDC_KIND_FEVAL, nvars, F%shape)
    do i = 1, nnodes*npieces
       call F%Fflt(i)%setval(0.0_pfdp)
    end do
    F%F(1:nnodes,1:npieces) => F%Fflt
    call F%ulevel%factory%create_array(F%I, nnodes-1, F%index, SDC_KIND_INTEGRAL, nvars, F%shape)
    call F%ulevel%factory%create_array(F%R, nnodes-1, F%index, SDC_KIND_INTEGRAL, nvars, F%shape)
    print *,'Finter in factory',F%Finterp, F%index
    if (F%index < pf%nlevels) then
        print *,'Finter in factory',F%Finterp, F%index
      if (F%Finterp) then
          call F%ulevel%factory%create_array(F%pFflt, nnodes*npieces, F%index, SDC_KIND_FEVAL, nvars, F%shape)
          F%pF(1:nnodes,1:npieces) => F%pFflt
       end if
       call F%ulevel%factory%create_array(F%pQ, nnodes, F%index, SDC_KIND_SOL_NO_FEVAL, nvars, F%shape)
    end if
    call F%ulevel%factory%create_single(F%qend, F%index, SDC_KIND_FEVAL, nvars, F%shape)
    call F%ulevel%factory%create_single(F%q0, F%index, SDC_KIND_FEVAL, nvars, F%shape)


  end subroutine pf_level_setup

  !
  ! Deallocate PFASST object
  !
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: l
    do l = 1, pf%nlevels
       call pf_level_destroy(pf%levels(l),pf%nlevels)
    end do
    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)
    deallocate(pf%state)
  end subroutine pf_pfasst_destroy

  !
  ! Deallocate PFASST level
  !
  subroutine pf_level_destroy(F,nlevels)
    class(pf_level_t), intent(inout) :: F
    integer                          :: nlevels, npieces

    if (.not. F%allocated) return

    ! flat buffers
    deallocate(F%send)
    deallocate(F%recv)

    ! nodes, flags, and integration matrices
    deallocate(F%nodes)
    deallocate(F%nflags)
    deallocate(F%qmat)
    deallocate(F%qmatFE)
    deallocate(F%qmatBE)
    deallocate(F%s0mat)
    deallocate(F%LUmat)

    ! encaps
    npieces = F%ulevel%sweeper%npieces

    if ((F%index < nlevels) .and. allocated(F%tauQ)) then
       call F%ulevel%factory%destroy_array(F%tauQ, F%nnodes-1, F%index, SDC_KIND_INTEGRAL, F%nvars, F%shape)
    end if

    call F%ulevel%factory%destroy_array(F%Q, F%nnodes, F%index, SDC_KIND_SOL_FEVAL, F%nvars, F%shape)
    call F%ulevel%factory%destroy_array(F%Fflt, F%nnodes*npieces, F%index, SDC_KIND_FEVAL, F%nvars, F%shape)
    call F%ulevel%factory%destroy_array(F%I, F%nnodes-1, F%index, SDC_KIND_INTEGRAL, F%nvars, F%shape)
    call F%ulevel%factory%destroy_array(F%R, F%nnodes-1, F%index, SDC_KIND_INTEGRAL, F%nvars, F%shape)
    if (F%index < nlevels) then
       if (F%Finterp) then
          call F%ulevel%factory%destroy_array(F%pFflt, F%nnodes*npieces, F%index, SDC_KIND_FEVAL, F%nvars, F%shape)
       end if
       call F%ulevel%factory%destroy_array(F%pQ, F%nnodes, F%index, SDC_KIND_SOL_NO_FEVAL, F%nvars, F%shape)
    end if
    call F%ulevel%factory%destroy_single(F%qend, F%index, SDC_KIND_FEVAL, F%nvars, F%shape)
    call F%ulevel%factory%destroy_single(F%q0, F%index, SDC_KIND_FEVAL, F%nvars, F%shape)

    ! destroy the sweeper 
    call F%ulevel%sweeper%destroy(F)

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
