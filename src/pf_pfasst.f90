!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
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
    integer :: l
    
    if (present(nlevels)) pf%nlevels = nlevels

    pf%outdir = ""

    ! gather some input from a file and command line
    read_cmd = .true.
    if (present(nocmd)) then
         if (nocmd)  read_cmd = .false.
    end if
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

  end subroutine pf_pfasst_create


  !
  ! Create a PFASST level object
  !
  subroutine pf_level_create(Lev, nlevel)
    type(pf_level_t), intent(inout) :: Lev
    integer,          intent(in   ) :: nlevel

    Lev%level = nlevel
    nullify(Lev%encap)
    nullify(Lev%interpolate)
    nullify(Lev%restrict)
    nullify(Lev%shape)
    nullify(Lev%tau)
    nullify(Lev%tauQ)
    nullify(Lev%pF)
    nullify(Lev%pQ)
    nullify(Lev%rmat)
    nullify(Lev%tmat)
    allocate(Lev%sweeper)
    Lev%ctx = c_null_ptr
  end subroutine pf_level_create


  !
  ! Setup PFASST object
  !
  subroutine pf_pfasst_setup(pf)
    use pf_mod_utils

    type(pf_pfasst_t), intent(inout) :: pf

    type(pf_level_t), pointer :: Flev, Glev
    integer                   :: l

    if (pf%rank < 0) then
       stop 'invalid PF rank: did you call setup correctly?'
    end if

    do l = 1, pf%nlevels
       call pf_level_setup(pf, pf%levels(l))
    end do

    do l = pf%nlevels, 2, -1
       Flev => pf%levels(l)
       Glev => pf%levels(l-1)

       allocate(Flev%tmat(Flev%nnodes,Glev%nnodes))
       call pf_time_interpolation_matrix(Flev%nodes, Flev%nnodes, Glev%nodes, Glev%nnodes, Flev%tmat)

       allocate(Flev%rmat(Glev%nnodes,Flev%nnodes))
       call pf_time_interpolation_matrix(Glev%nodes, Glev%nnodes, Flev%nodes, Flev%nnodes, Flev%rmat)
    end do

  end subroutine pf_pfasst_setup

  !
  ! Setup (allocate) PFASST level
  !
  ! If the level is already setup, calling this again will allocate
  ! (or deallocate) tau appropriately.
  !
  subroutine pf_level_setup(pf, Lev)
    use pf_mod_quadrature
    type(pf_pfasst_t), intent(in   ) :: pf
    type(pf_level_t),  intent(inout) :: Lev

    integer :: m, p, nvars, nnodes, npieces

    !
    ! do some sanity checks
    !

    if (Lev%nvars <= 0) stop "ERROR: Invalid nvars/dofs (pf_pfasst.f90)."
    if (Lev%nnodes <= 0) stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (Lev%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."
    if (.not. associated(Lev%encap)) stop "ERROR: Missing encapsulation (pf_pfasst.f90)."
    if (.not. associated(Lev%interpolate)) stop "ERROR: Missing spatial interpolation (pf_pfasst.f90)."
    if (.not. associated(Lev%restrict)) stop "ERROR: Missing spatial restriction (pf_pfasst.f90)."

    nvars   = Lev%nvars
    nnodes  = Lev%nnodes
    npieces = Lev%sweeper%npieces

    Lev%residual = -1.0_pfdp

    !
    ! (re)allocate tau (may to need create/destroy tau dynamically
    !                   when doing AMR)
    !
    if ((Lev%level < pf%nlevels) .and. (.not. associated(Lev%tau))) then
       allocate(Lev%tau(nnodes-1))
       do m = 1, nnodes-1
          call Lev%encap%create(Lev%tau(m), Lev%level, SDC_KIND_INTEGRAL, &
               nvars, Lev%shape, Lev%ctx)
       end do
    else if ((Lev%level >= pf%nlevels) .and. (associated(Lev%tau))) then
       do m = 1, nnodes-1
          call Lev%encap%destroy(Lev%tau(m))
       end do
       deallocate(Lev%tau)
       nullify(Lev%tau)
    end if

    !
    ! (re)allocate tauQ (may to need create/destroy tau dynamically
    !                   when doing AMR)
    !
    if ((Lev%level < pf%nlevels) .and. (.not. associated(Lev%tauQ))) then
       allocate(Lev%tauQ(nnodes-1))
       do m = 1, nnodes-1
          call Lev%encap%create(Lev%tauQ(m), Lev%level, SDC_KIND_INTEGRAL, &
               nvars, Lev%shape, Lev%ctx)
       end do
    else if ((Lev%level >= pf%nlevels) .and. (associated(Lev%tauQ))) then
       do m = 1, nnodes-1
          call Lev%encap%destroy(Lev%tauQ(m))
       end do
       deallocate(Lev%tauQ)
       nullify(Lev%tauQ)
    end if

    !
    ! skip the rest if we're already allocated
    !
    if (Lev%allocated) return
    Lev%allocated = .true.

    !
    ! allocate flat buffers (send, and recv)
    !
    allocate(Lev%send(nvars))
    allocate(Lev%recv(nvars))

    !
    ! nodes, flags, and integration matrices
    !
    allocate(Lev%nodes(nnodes))
    allocate(Lev%nflags(nnodes))
    allocate(Lev%s0mat(nnodes-1,nnodes))
    allocate(Lev%qmat(nnodes-1,nnodes))

    if (btest(pf%qtype, 8)) then
       call pf_quadrature(pf%qtype, nnodes, pf%levels(1)%nnodes, &
            Lev%nodes, Lev%nflags, Lev%s0mat, Lev%qmat)
    else
       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            Lev%nodes, Lev%nflags, Lev%s0mat, Lev%qmat)
    end if

    call Lev%sweeper%initialize(Lev)


    !
    ! allocate Q and F
    !
    allocate(Lev%Q(nnodes))
    allocate(Lev%F(nnodes,npieces))

    do m = 1, nnodes
       call Lev%encap%create(Lev%Q(m), Lev%level, SDC_KIND_SOL_FEVAL, &
            nvars, Lev%shape, Lev%ctx)
       do p = 1, npieces
          call Lev%encap%create(Lev%F(m,p), Lev%level, SDC_KIND_FEVAL, &
               nvars, Lev%shape, Lev%ctx)
       end do
    end do

    !
    ! allocate S, I, and R
    !
    allocate(Lev%S(nnodes-1))
    allocate(Lev%I(nnodes-1))
    allocate(Lev%R(nnodes-1))

    do m = 1, nnodes-1
       call Lev%encap%create(Lev%S(m), Lev%level, SDC_KIND_INTEGRAL, &
            nvars, Lev%shape, Lev%ctx)
       call Lev%encap%create(Lev%I(m), Lev%level, SDC_KIND_INTEGRAL, &
            nvars, Lev%shape, Lev%ctx)
       call Lev%encap%create(Lev%R(m), Lev%level, SDC_KIND_INTEGRAL, &
            nvars, Lev%shape, Lev%ctx)
    end do

    !
    ! allocate pQ and pF
    !
    if (Lev%level < pf%nlevels) then

       if (Lev%Finterp) then
          ! store F and Q(1) only
          ! Changed by MM Dec. 20, 2013 to allocate all pQ as well
          !
          allocate(Lev%pF(nnodes,npieces))
          allocate(Lev%pQ(nnodes))
          do m = 1, nnodes
             do p = 1, npieces
                call Lev%encap%create(Lev%pF(m,p), Lev%level, SDC_KIND_FEVAL, &
                     nvars, Lev%shape, Lev%ctx)
             end do
             call Lev%encap%create(Lev%pQ(m), Lev%level, SDC_KIND_SOL_NO_FEVAL, &
                  nvars, Lev%shape, Lev%ctx)
          end do
       else
          ! store Q
          allocate(Lev%pQ(nnodes))
          do m = 1, nnodes
             call Lev%encap%create(Lev%pQ(m), Lev%level, SDC_KIND_SOL_NO_FEVAL, &
                  nvars, Lev%shape, Lev%ctx)
          end do
       end if

    end if

    !
    ! allocate qend and q0
    !
    call Lev%encap%create(Lev%qend, Lev%level, SDC_KIND_SOL_NO_FEVAL, &
         nvars, Lev%shape, Lev%ctx)
    call Lev%encap%create(Lev%q0, Lev%level, SDC_KIND_SOL_NO_FEVAL, &
         nvars, Lev%shape, Lev%ctx)

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
  subroutine pf_level_destroy(Lev)
    type(pf_level_t), intent(inout) :: Lev

    integer :: m, p

    if (.not. Lev%allocated) return

    ! flat buffers
    deallocate(Lev%send)
    deallocate(Lev%recv)

    ! nodes, flags, and integration matrices
    deallocate(Lev%nodes)
    deallocate(Lev%nflags)
    deallocate(Lev%qmat)
    deallocate(Lev%s0mat)

    ! Q and F
    do m = 1, Lev%nnodes
       call Lev%encap%destroy(Lev%Q(m))
       do p = 1, size(Lev%F(m,:))
          call Lev%encap%destroy(Lev%F(m,p))
       end do
    end do
    deallocate(Lev%Q)
    deallocate(Lev%F)

    ! S, I, and R
    do m = 1, Lev%nnodes-1
       call Lev%encap%destroy(Lev%S(m))
       call Lev%encap%destroy(Lev%I(m))
       call Lev%encap%destroy(Lev%R(m))
    end do
    deallocate(Lev%S)
    deallocate(Lev%I)
    deallocate(Lev%R)

    ! pQ and pF
    if (associated(Lev%pQ)) then
       if (Lev%Finterp) then
          do m = 1, Lev%nnodes
             do p = 1, size(Lev%F(m,:))
                call Lev%encap%destroy(Lev%pF(m,p))
             end do
          end do
          call Lev%encap%destroy(Lev%pQ(1))
          deallocate(Lev%pF)
          deallocate(Lev%pQ)
       else
          do m = 1, Lev%nnodes
             call Lev%encap%destroy(Lev%pQ(m))
          end do
          deallocate(Lev%pQ)
       end if
       nullify(Lev%pQ)
    end if

    ! qend
    call Lev%encap%destroy(Lev%qend)
    call Lev%encap%destroy(Lev%q0)


    ! tau
    if (associated(Lev%tau)) then
       do m = 1, Lev%nnodes-1
          call Lev%encap%destroy(Lev%tau(m))
       end do
       deallocate(Lev%tau)
       nullify(Lev%tau)
    end if

    ! tauQ
    if (associated(Lev%tauQ)) then
       do m = 1, Lev%nnodes-1
          call Lev%encap%destroy(Lev%tauQ(m))
       end do
       deallocate(Lev%tauQ)
       nullify(Lev%tauQ)
    end if

    ! other
    if (associated(Lev%shape)) then
       deallocate(Lev%shape)
    end if

    if (associated(Lev%tmat)) then
       deallocate(Lev%tmat)
       nullify(Lev%tmat)
    end if

    if (associated(Lev%rmat)) then
       deallocate(Lev%rmat)
       nullify(Lev%rmat)
    end if

    ! kill the sweeper
    call Lev%sweeper%destroy(Lev%sweeper)
    deallocate(Lev%sweeper)

  end subroutine pf_level_destroy

end module pf_mod_pfasst
