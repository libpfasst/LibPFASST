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

  ! Create a PFASST object
  !
  ! Passing the 'nvars' and 'nnodes' arrays is optional, but these
  ! should be set appropriately before calling setup.
  subroutine pf_pfasst_create(pf, comm, encap, sweeper, nlevels, nvars, nnodes, maxlevels)
    type(pf_pfasst_t),   intent(inout)         :: pf
    type(pf_comm_t),     intent(inout), target :: comm
    type(pf_encap_t),    intent(in)            :: encap
    type(pf_sweeper_t),  intent(in)            :: sweeper
    integer,             intent(in)            :: nlevels
    integer,             intent(in), optional  :: nvars(nlevels), nnodes(nlevels), maxlevels

    integer :: l, nlevs

    if (present(maxlevels)) then
       nlevs = maxlevels
    else
       nlevs = nlevels
    end if

    pf%comm => comm

    pf%nlevels = nlevels
    allocate(pf%levels(nlevs))
    allocate(pf%hooks(nlevs,PF_MAX_HOOKS))
    allocate(pf%nhooks(nlevs))

    pf%nhooks = 0

    do l = 1, nlevels
       call pf_level_create(pf%levels(l), l)
       pf%levels(l)%sweeper = sweeper
       pf%levels(l)%encap = encap
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
    level%tau => null()
    level%shape => null()
    level%smat => null()
    level%pfSDC => null()
    level%pSDC => null()
    level%tmat => null()
    level%rmat => null()
  end subroutine pf_level_create


  ! Setup (allocate) PFASST object
  subroutine pf_pfasst_setup(pf)
    use pf_mod_utils
    use pf_mod_version

    type(pf_pfasst_t), intent(inout) :: pf

    type(pf_level_t), pointer :: F, G
    integer :: l

    if (pf%rank < 0) then
       stop 'invalid PF rank: did you call setup correctly?'
    end if

    do l = 1, pf%nlevels
       call setup(pf, pf%levels(l))
    end do

    do l = pf%nlevels, 2, -1
       F => pf%levels(l)
       G => pf%levels(l-1)

       allocate(F%tmat(F%nnodes,G%nnodes))
       call pf_time_interpolation_matrix(F%nodes, F%nnodes, G%nodes, G%nnodes, F%tmat)

       allocate(F%rmat(G%nnodes,F%nnodes))
       call pf_time_interpolation_matrix(G%nodes, G%nnodes, F%nodes, F%nnodes, F%rmat)
    end do

  end subroutine pf_pfasst_setup


  ! Setup (allocate) PFASST object
  !
  ! If the level is already setup, calling this again will allocate
  ! (or deallocate) tau appropriately.
  subroutine pf_level_setup(pf, F)
    use pf_mod_quadrature

    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: F

    integer :: m, p, nvars, nnodes, npieces

    nvars   = F%nvars
    nnodes  = F%nnodes
    npieces = F%sweeper%npieces

    if (.not. F%allocated) then

       allocate(F%q0(nvars))
       allocate(F%send(nvars))
       allocate(F%recv(nvars))
       allocate(F%qSDC(nnodes))
       allocate(F%fSDC(nnodes,npieces))
       allocate(F%nodes(nnodes))
       allocate(F%nflags(nnodes))
       allocate(F%s0mat(nnodes-1,nnodes))
       allocate(F%qmat(nnodes-1,nnodes))

       if (F%Finterp) then
          if (F%level < pf%nlevels) then
             allocate(F%pfSDC(nnodes,npieces))
             allocate(F%pSDC(1))
          end if
       else
          if (F%level < pf%nlevels) then
             allocate(F%pSDC(nnodes))
          end if
       end if

       call F%encap%create(F%qend, F%level, .false., nvars, F%shape, F%ctx)
       ! call F%encap%create(F%qex, F%level, .false., nvars, F%shape, F%ctx)

       do m = 1, nnodes
          call F%encap%create(F%qSDC(m), F%level, .false., nvars, F%shape, F%ctx)
          do p = 1, npieces
             call F%encap%create(F%fSDC(m,p), F%level, .true., nvars, F%shape, F%ctx)
          end do
       end do

       ! create space to store previous iteration info
       if (F%level < pf%nlevels) then

          if (F%Finterp) then  !  Doing store of f and qSDC(1) only
             do m = 1, nnodes
                do p = 1, npieces
                   call F%encap%create(F%pfSDC(m,p), F%level, .true., nvars, F%shape, F%ctx)
                end do
             end do
             call F%encap%create(F%pSDC(1), F%level, .false., nvars, F%shape, F%ctx)
          else   !  Storing all qSDC
             do m = 1, nnodes
                call F%encap%create(F%pSDC(m), F%level, .false., nvars, F%shape, F%ctx)
             end do
          end if

       end if

       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            F%nodes, F%nflags, F%s0mat, F%qmat)

       call F%sweeper%initialize(F)

    end if

    if ((F%level < pf%nlevels) .and. (.not. associated(F%tau))) then
       allocate(F%tau(nnodes-1))
       do m = 1, nnodes-1
          call F%encap%create(F%tau(m), F%level, .false., nvars, F%shape, F%ctx)
       end do
    else if ((F%level >= pf%nlevels) .and. (associated(F%tau))) then
       do m = 1, nnodes-1
          call F%encap%destroy(F%tau(m))
       end do
       deallocate(F%tau)
       nullify(F%tau)
    end if

    F%allocated = .true.
  end subroutine pf_level_setup

  ! Deallocate PFASST object
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: l

    do l = 1, pf%nlevels
       call pf_level_destroy(pf%levels(l))
    end do

    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)

    ! if (pf%log > 0) then
    !    close(pf%log)
    ! end if
  end subroutine pf_pfasst_destroy

  ! Deallocate PFASST level object
  subroutine pf_level_destroy(F)
    type(pf_level_t), intent(inout) :: F

    integer :: m, p

    if (F%allocated) then
       deallocate(F%q0)
       deallocate(F%send)
       deallocate(F%recv)
       deallocate(F%nodes)
       deallocate(F%nflags)
       deallocate(F%qmat)
       deallocate(F%s0mat)

       if (associated(F%shape)) then
          deallocate(F%shape)
       end if

       if (associated(F%smat)) then
          deallocate(F%smat)
       end if

       call F%encap%destroy(F%qend)
       ! call F%encap%destroy(F%qex)

       do m = 1, F%nnodes
          call F%encap%destroy(F%qSDC(m))
          do p = 1, size(F%fSDC(m,:))
             call F%encap%destroy(F%fSDC(m,p))
          end do
       end do

       if (F%Finterp) then
          do m = 1, F%nnodes
             do p = 1, size(F%fSDC(m,:))
                if (associated(F%pfSDC)) then
                   call F%encap%destroy(F%pfSDC(m,p))
                end if
             end do
          end do
          if (associated(F%pSDC)) then
             call F%encap%destroy(F%pSDC(1))
          end if
       else
          do m = 1, F%nnodes
             if (associated(F%pSDC)) then
                call F%encap%destroy(F%pSDC(m))
             end if
          end do
       end if


       deallocate(F%qSDC)
       deallocate(F%fSDC)
       if (F%Finterp) then
          if (associated(F%pfSDC)) deallocate(F%pfSDC)
       else
          if (associated(F%pSDC)) deallocate(F%pSDC)
       end if

       if (associated(F%tau)) then
          do m = 1, F%nnodes-1
             call F%encap%destroy(F%tau(m))
          end do
          deallocate(F%tau)
          nullify(F%tau)
       end if

       if (associated(F%tmat)) then
          deallocate(F%tmat)
          nullify(F%tmat)
       end if

       if (associated(F%rmat)) then
          deallocate(F%rmat)
          nullify(F%rmat)
       end if

    end if

  end subroutine pf_level_destroy

end module pf_mod_pfasst
