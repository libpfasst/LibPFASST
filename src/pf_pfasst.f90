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

    logical :: read_cmd = .true.
    integer :: l

    if (present(nlevels)) pf%nlevels = nlevels

    ! gather some input from a file and command line
    if (present(nocmd) .and. nocmd)   read_cmd = .false.
    if (present(fname))  then
       call pf_read_opts(pf, read_cmd, fname)
    else
       if (read_cmd)   call pf_read_opts(pf, read_cmd)
    end if

    pf%comm => comm

    allocate(pf%levels(pf%nlevels))
    allocate(pf%hooks(pf%nlevels, PF_MAX_HOOK, PF_MAX_HOOKS))
    allocate(pf%nhooks(pf%nlevels, PF_MAX_HOOK))
    pf%nhooks = 0

    do l = 1, pf%nlevels
       call pf_level_create(pf%levels(l), l)
    end do

    allocate(pf%state)
    pf%state%pstatus = 0
    pf%state%status  = 0

    nullify(pf%cycles%start)
    nullify(pf%cycles%pfasst)
  end subroutine pf_pfasst_create


  !
  ! Create a PFASST level object
  !
  subroutine pf_level_create(level, nlevel)
    type(pf_level_t), intent(inout) :: level
    integer,          intent(in   ) :: nlevel

    level%level = nlevel
    nullify(level%encap)
    nullify(level%interpolate)
    nullify(level%restrict)
    nullify(level%shape)
    nullify(level%tau)
    nullify(level%pF)
    nullify(level%pQ)
    nullify(level%rmat)
    nullify(level%tmat)
    allocate(level%sweeper)
    level%levelctx = c_null_ptr
  end subroutine pf_level_create


  !
  ! Setup PFASST object
  !
  subroutine pf_pfasst_setup(pf)
    use pf_mod_utils
    use pf_mod_cycle

    type(pf_pfasst_t), intent(inout) :: pf

    type(pf_level_t), pointer :: F, G
    integer                   :: l

    if (pf%rank < 0) then
       stop 'invalid PF rank: did you call setup correctly?'
    end if

    do l = 1, pf%nlevels
       call pf_level_setup(pf, pf%levels(l))
    end do

    do l = pf%nlevels, 2, -1
       F => pf%levels(l)
       G => pf%levels(l-1)

       allocate(F%tmat(F%nnodes,G%nnodes))
       call pf_time_interpolation_matrix(F%nodes, F%nnodes, G%nodes, G%nnodes, F%tmat)

       allocate(F%rmat(G%nnodes,F%nnodes))
       call pf_time_interpolation_matrix(G%nodes, G%nnodes, F%nodes, F%nnodes, F%rmat)
    end do

    call pf_cycle_build(pf)

  end subroutine pf_pfasst_setup

  !
  ! Setup (allocate) PFASST level
  !
  ! If the level is already setup, calling this again will allocate
  ! (or deallocate) tau appropriately.
  !
  subroutine pf_level_setup(pf, F)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(in   ) :: pf
    type(pf_level_t),  intent(inout) :: F

    integer :: m, p, nvars, nnodes, npieces

    !
    ! do some sanity checks
    !

    if (F%nvars <= 0) stop "ERROR: Invalid nvars/dofs (pf_pfasst.f90)."
    if (F%nnodes <= 0) stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (F%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."
    if (.not. associated(F%encap)) stop "ERROR: Missing encapsulation (pf_pfasst.f90)."
    if (.not. associated(F%interpolate)) stop "ERROR: Missing spatial interpolation (pf_pfasst.f90)."
    if (.not. associated(F%restrict)) stop "ERROR: Missing spatial restriction (pf_pfasst.f90)."

    nvars   = F%nvars
    nnodes  = F%nnodes
    npieces = F%sweeper%npieces

    F%residual = -1.0_pfdp

    !
    ! (re)allocate tau (may to need create/destroy tau dynamically
    !                   when doing AMR)
    !
    if ((F%level < pf%nlevels) .and. (.not. associated(F%tau))) then
       allocate(F%tau(nnodes-1))
       do m = 1, nnodes-1
          call F%encap%create(F%tau(m), F%level, SDC_KIND_INTEGRAL, &
               nvars, F%shape, F%levelctx, F%encap%encapctx)
       end do
    else if ((F%level >= pf%nlevels) .and. (associated(F%tau))) then
       do m = 1, nnodes-1
          call F%encap%destroy(F%tau(m))
       end do
       deallocate(F%tau)
       nullify(F%tau)
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

    call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
         F%nodes, F%nflags, F%s0mat, F%qmat)

    call F%sweeper%initialize(F)

    !
    ! allocate Q and F
    !
    allocate(F%Q(nnodes))
    allocate(F%F(nnodes,npieces))

    do m = 1, nnodes
       call F%encap%create(F%Q(m), F%level, SDC_KIND_SOL_FEVAL, &
            nvars, F%shape, F%levelctx, F%encap%encapctx)
       do p = 1, npieces
          call F%encap%create(F%F(m,p), F%level, SDC_KIND_FEVAL, &
               nvars, F%shape, F%levelctx, F%encap%encapctx)
       end do
    end do

    !
    ! allocate S, I, and R
    !
    allocate(F%S(nnodes-1))
    allocate(F%I(nnodes-1))
    allocate(F%R(nnodes-1))

    do m = 1, nnodes-1
       call F%encap%create(F%S(m), F%level, SDC_KIND_INTEGRAL, &
            nvars, F%shape, F%levelctx, F%encap%encapctx)
       call F%encap%create(F%I(m), F%level, SDC_KIND_INTEGRAL, &
            nvars, F%shape, F%levelctx, F%encap%encapctx)
       call F%encap%create(F%R(m), F%level, SDC_KIND_INTEGRAL, &
            nvars, F%shape, F%levelctx, F%encap%encapctx)
    end do

    !
    ! allocate pQ and pF
    !
    if (F%level < pf%nlevels) then

       if (F%Finterp) then
          ! store F and Q(1) only
          allocate(F%pF(nnodes,npieces))
          allocate(F%pQ(1))
          do m = 1, nnodes
             do p = 1, npieces
                call F%encap%create(F%pF(m,p), F%level, SDC_KIND_FEVAL, &
                     nvars, F%shape, F%levelctx, F%encap%encapctx)
             end do
          end do
          call F%encap%create(F%pQ(1), F%level, SDC_KIND_SOL_NO_FEVAL, &
               nvars, F%shape, F%levelctx, F%encap%encapctx)
       else
          ! store Q
          allocate(F%pQ(nnodes))
          do m = 1, nnodes
             call F%encap%create(F%pQ(m), F%level, SDC_KIND_SOL_NO_FEVAL, &
                  nvars, F%shape, F%levelctx, F%encap%encapctx)
          end do
       end if

    end if

    !
    ! allocate Qend
    !
    call F%encap%create(F%qend, F%level, SDC_KIND_SOL_NO_FEVAL, &
         nvars, F%shape, F%levelctx, F%encap%encapctx)

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

    if (associated(pf%cycles%start))  deallocate(pf%cycles%start)
    if (associated(pf%cycles%pfasst)) deallocate(pf%cycles%pfasst)

    deallocate(pf%state)
  end subroutine pf_pfasst_destroy


  !
  ! Deallocate PFASST level
  !
  subroutine pf_level_destroy(F)
    type(pf_level_t), intent(inout) :: F

    integer :: m, p

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

    ! Q and F
    do m = 1, F%nnodes
       call F%encap%destroy(F%Q(m))
       do p = 1, size(F%F(m,:))
          call F%encap%destroy(F%F(m,p))
       end do
    end do
    deallocate(F%Q)
    deallocate(F%F)

    ! S, I, and R
    do m = 1, F%nnodes-1
       call F%encap%destroy(F%S(m))
       call F%encap%destroy(F%I(m))
       call F%encap%destroy(F%R(m))
    end do
    deallocate(F%S)
    deallocate(F%I)
    deallocate(F%R)

    ! pQ and pF
    if (associated(F%pQ)) then
       if (F%Finterp) then
          do m = 1, F%nnodes
             do p = 1, size(F%F(m,:))
                call F%encap%destroy(F%pF(m,p))
             end do
          end do
          call F%encap%destroy(F%pQ(1))
          deallocate(F%pF)
          deallocate(F%pQ)
       else
          do m = 1, F%nnodes
             call F%encap%destroy(F%pQ(m))
          end do
          deallocate(F%pQ)
       end if
       nullify(F%pQ)
    end if

    ! qend
    call F%encap%destroy(F%qend)

    ! tau
    if (associated(F%tau)) then
       do m = 1, F%nnodes-1
          call F%encap%destroy(F%tau(m))
       end do
       deallocate(F%tau)
       nullify(F%tau)
    end if

    ! other
    if (associated(F%shape)) then
       deallocate(F%shape)
    end if

    if (associated(F%tmat)) then
       deallocate(F%tmat)
       nullify(F%tmat)
    end if

    if (associated(F%rmat)) then
       deallocate(F%rmat)
       nullify(F%rmat)
    end if

    ! kill the sweeper
    call F%sweeper%destroy(F%sweeper)
    deallocate(F%sweeper)

  end subroutine pf_level_destroy

end module pf_mod_pfasst
