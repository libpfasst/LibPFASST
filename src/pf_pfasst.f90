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

  interface create
     module procedure pf_pfasst_create
     module procedure pf_level_create
  end interface create

  interface setup
     module procedure pf_pfasst_setup
     module procedure pf_level_setup
  end interface setup

  interface destroy
     module procedure pf_pfasst_destroy
     module procedure pf_level_destroy
  end interface destroy

contains

  ! Create a PFASST object
  !
  ! Passing the 'nvars' and 'nnodes' arrays is optional, but these
  ! should be set appropriately before calling setup.
  subroutine pf_pfasst_create(pf, pf_comm, nlevels, nvars, nnodes, maxlevels)
    type(pf_pfasst_t), intent(inout)         :: pf
    type(pf_comm_t),   intent(inout), target :: pf_comm
    integer,           intent(in)            :: nlevels
    integer,           intent(in), optional  :: nvars(nlevels), nnodes(nlevels), maxlevels

    integer :: l, nlevs

    if (present(maxlevels)) then
       nlevs = maxlevels
    else
       nlevs = nlevels
    end if

    pf%comm => pf_comm

    pf%nlevels = nlevels
    allocate(pf%levels(nlevs))
    allocate(pf%hooks(nlevs,PF_MAX_HOOKS))
    allocate(pf%nhooks(nlevs))

    pf%nhooks = 0

    do l = 1, nlevels
       call create(pf%levels(l), l)
    end do

    if (present(nvars)) then
       do l = 1, nlevels
          pf%levels(l)%nvars = nvars(l)
       end do
    end if

    if (present(nnodes)) then
       do l = 1, nlevels
          pf%levels(l)%nnodes = nnodes(l)
       end do
    end if
  end subroutine pf_pfasst_create


  ! Create a PFASST level object
  subroutine pf_level_create(level, nlevel)
    type(pf_level_t), intent(inout) :: level
    integer,          intent(in)    :: nlevel

    level%level = nlevel
  end subroutine pf_level_create


  ! Setup (allocate) PFASST object
  subroutine pf_pfasst_setup(pf, logunit, logname)
    use pf_mod_utils
    use pf_mod_version

    type(pf_pfasst_t), intent(inout)        :: pf
    integer,           intent(in), optional :: logunit
    character(len=*),  intent(in), optional :: logname

    type(pf_level_t), pointer :: F, G
    integer :: l

    character(len=128) :: fname

    if (pf%rank < 0) then
       stop 'invalid PF rank: did you call setup correctly?'
    end if

    do l = 1, pf%nlevels
       call setup(pf%levels(l), pf)
    end do

    do l = pf%nlevels, 2, -1
       F => pf%levels(l)
       G => pf%levels(l-1)

       allocate(F%tmat(F%nnodes,G%nnodes))
       call pf_time_interpolation_matrix(F%nodes, F%nnodes, G%nodes, G%nnodes, F%tmat)

       allocate(F%rmat(G%nnodes,F%nnodes))
       call pf_time_interpolation_matrix(G%nodes, G%nnodes, F%nodes, F%nnodes, F%rmat)
    end do

    if (present(logunit)) then
       pf%log = logunit
       write (fname, "(a,i0)"), logname // '.', pf%rank
       open(unit=pf%log, file=fname)

       write (pf%log, *) 'PFASST RUN STARTED AT: ', ctime(time8())
       write (pf%log, *) 'PFASST GIT VERSION:    ', PF_GIT_VERSION
    end if
  end subroutine pf_pfasst_setup


  ! Setup (allocate) PFASST object
  !
  ! If the level is already setup, calling this again will allocate
  ! (or deallocate) tau appropriately.
  subroutine pf_level_setup(level, pf)
    use pf_mod_quadrature
    use pf_mod_sweep
    use pf_mod_sweep, only: npieces
    type(pf_level_t),  intent(inout) :: level
    type(pf_pfasst_t), intent(in)    :: pf

    integer :: m, p, nvars, nnodes

    nvars   = level%nvars
    nnodes  = level%nnodes

    if (.not. level%allocated) then

       allocate(level%q0(nvars))
       allocate(level%send(nvars))
       allocate(level%recv(nvars))
       allocate(level%qSDC(nnodes))
       allocate(level%fSDC(nnodes,npieces))
       allocate(level%nodes(nnodes))
       allocate(level%nmask(nnodes))
       allocate(level%s0mat(nnodes-1,nnodes))
       allocate(level%qmat(nnodes-1,nnodes))

       if (level%Finterp) then
          if (level%level < pf%nlevels) then
             allocate(level%pfSDC(nnodes,npieces))
             allocate(level%pSDC(1))
          end if
       else
          if (level%level < pf%nlevels) then
             allocate(level%pSDC(nnodes))
          end if
       end if



       call create(level%qend, level%level, .false., nvars, level%shape, level%ctx)
       call create(level%qex, level%level, .false., nvars, level%shape, level%ctx)

       do m = 1, nnodes
          call create(level%qSDC(m), level%level, .false., nvars, level%shape, level%ctx)
          do p = 1, npieces
             call create(level%fSDC(m,p), level%level, .true., nvars, level%shape, level%ctx)
          end do
       end do

       ! Create space to store previous iteration info
       if (level%level < pf%nlevels) then
          if (level%Finterp) then  !  Doing store of f and qSDC(1) only
             do m = 1, nnodes
                do p = 1, npieces
                   call create(level%pfSDC(m,p), level%level, .true., nvars, level%shape, level%ctx)
                end do
             end do
             call create(level%pSDC(1), level%level, .false., nvars, level%shape, level%ctx)
          else   !  Storing all qSDC
             do m = 1, nnodes
                call create(level%pSDC(m), level%level, .false., nvars, level%shape, level%ctx)
             end do
          end if
       end if


       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            level%nodes, level%nmask, level%s0mat, level%qmat)

       call sdcinit(level)

    end if

    if ((level%level < pf%nlevels) .and. (.not. associated(level%tau))) then
       allocate(level%tau(nnodes-1))
       do m = 1, nnodes-1
          call create(level%tau(m), level%level, .false., nvars, level%shape, level%ctx)
       end do
    else if ((level%level >= pf%nlevels) .and. (associated(level%tau))) then
       do m = 1, nnodes-1
          call destroy(level%tau(m))
       end do
       deallocate(level%tau)
       nullify(level%tau)
    end if

    level%allocated = .true.
  end subroutine pf_level_setup

  ! Deallocate PFASST object
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: l

    do l = 1, pf%nlevels
       call destroy(pf%levels(l))
    end do

    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)

    if (pf%log > 0) then
       close(pf%log)
    end if
  end subroutine pf_pfasst_destroy

  ! Deallocate PFASST level object
  subroutine pf_level_destroy(level)
    type(pf_level_t), intent(inout) :: level

    integer :: m, p

    if (level%allocated) then
       deallocate(level%q0)
       deallocate(level%send)
       deallocate(level%recv)
       deallocate(level%nodes)
       deallocate(level%nmask)
       deallocate(level%qmat)
       deallocate(level%s0mat)

       if (associated(level%shape)) then
          deallocate(level%shape)
       end if

       if (associated(level%smat)) then
          deallocate(level%smat)
       end if

       call destroy(level%qend)
       call destroy(level%qex)

       do m = 1, level%nnodes
          call destroy(level%qSDC(m))
          do p = 1, size(level%fSDC(m,:))
             call destroy(level%fSDC(m,p))
          end do
       end do

       if (level%Finterp) then
          do m = 1, level%nnodes
             do p = 1, size(level%fSDC(m,:))
                if (associated(level%pfSDC)) then
                   call destroy(level%pfSDC(m,p))
                end if
             end do
          end do
          if (associated(level%pSDC)) then
             call destroy(level%pSDC(1))
          end if
       else
          do m = 1, level%nnodes
             if (associated(level%pSDC)) then
                call destroy(level%pSDC(m))
             end if
          end do
       end if


       deallocate(level%qSDC)
       deallocate(level%fSDC)
       if (level%Finterp) then
          if (associated(level%pfSDC)) deallocate(level%pfSDC)
       else
          if (associated(level%pSDC)) deallocate(level%pSDC)
       end if

       if (associated(level%tau)) then
          do m = 1, level%nnodes-1
             call destroy(level%tau(m))
          end do
          deallocate(level%tau)
          nullify(level%tau)
       end if

       if (associated(level%tmat)) then
          deallocate(level%tmat)
          nullify(level%tmat)
       end if

       if (associated(level%rmat)) then
          deallocate(level%rmat)
          nullify(level%rmat)
       end if

    end if

    level%allocated = .false.

  end subroutine pf_level_destroy

end module pf_mod_pfasst
