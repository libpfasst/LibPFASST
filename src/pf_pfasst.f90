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
  subroutine pf_pfasst_create(pf, comm, maxlevels)
    use pf_mod_hooks, only: PF_MAX_HOOK
    type(pf_pfasst_t),   intent(inout)         :: pf
    type(pf_comm_t),     intent(inout), target :: comm
    integer,             intent(in)            :: maxlevels

    integer :: l

    pf%comm => comm
    pf%nlevels = 0

    allocate(pf%levels(maxlevels))
    allocate(pf%hooks(maxlevels, PF_MAX_HOOK, PF_MAX_HOOKS))
    allocate(pf%nhooks(maxlevels, PF_MAX_HOOK))

    pf%nhooks = 0

    do l = 1, maxlevels
       call pf_level_create(pf%levels(l), l)
    end do

    nullify(pf%cycles%start)
    nullify(pf%cycles%pfasst)
    nullify(pf%cycles%end)

    pf%nlevels = maxlevels
  end subroutine pf_pfasst_create


  ! Create a PFASST level object
  subroutine pf_level_create(level, nlevel)
    type(pf_level_t), intent(inout) :: level
    integer,          intent(in)    :: nlevel

    level%level = nlevel
    level%shape => null()
    level%tau   => null()
    level%pF    => null()
    level%pQ    => null()
    level%smat  => null()
    level%rmat  => null()
    level%tmat  => null()
  end subroutine pf_level_create


  ! Setup (allocate) PFASST object
  subroutine pf_pfasst_setup(pf)
    use pf_mod_utils

    type(pf_pfasst_t), intent(inout) :: pf

    type(pf_level_t), pointer :: F, G
    integer :: l

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
       allocate(F%Q(nnodes))
       allocate(F%F(nnodes,npieces))
       allocate(F%S(nnodes-1))
       allocate(F%nodes(nnodes))
       allocate(F%nflags(nnodes))
       allocate(F%s0mat(nnodes-1,nnodes))
       allocate(F%qmat(nnodes-1,nnodes))

       if (F%Finterp) then
          if (F%level < pf%nlevels) then
             allocate(F%pF(nnodes,npieces))
             allocate(F%pQ(1))
          end if
       else
          if (F%level < pf%nlevels) then
             allocate(F%pQ(nnodes))
          end if
       end if

       call F%encap%create(F%qend, F%level, SDC_KIND_SOL_NO_FEVAL, &
            nvars, F%shape, F%ctx, F%encap%ctx)

       do m = 1, nnodes
          call F%encap%create(F%Q(m), F%level, SDC_KIND_SOL_FEVAL, &
               nvars, F%shape, F%ctx, F%encap%ctx)
          do p = 1, npieces
             call F%encap%create(F%F(m,p), F%level, SDC_KIND_FEVAL, &
                  nvars, F%shape, F%ctx, F%encap%ctx)
          end do
       end do

       do m = 1, nnodes-1
          call F%encap%create(F%S(m), F%level, SDC_KIND_INTEGRAL, &
               nvars, F%shape, F%ctx, F%encap%ctx)
       end do


       ! create space to store previous iteration info
       if (F%level < pf%nlevels) then

          if (F%Finterp) then  !  Doing store of f and qSDC(1) only
             do m = 1, nnodes
                do p = 1, npieces
                   call F%encap%create(F%pF(m,p), F%level, SDC_KIND_FEVAL, &
                        nvars, F%shape, F%ctx, F%encap%ctx)
                end do
             end do
             call F%encap%create(F%pQ(1), F%level, SDC_KIND_SOL_NO_FEVAL, &
                  nvars, F%shape, F%ctx, F%encap%ctx)
          else   !  Storing all qSDC
             do m = 1, nnodes
                call F%encap%create(F%pQ(m), F%level, SDC_KIND_SOL_NO_FEVAL, &
                     nvars, F%shape, F%ctx, F%encap%ctx)
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
          call F%encap%create(F%tau(m), F%level, SDC_KIND_INTEGRAL, &
               nvars, F%shape, F%ctx, F%encap%ctx)
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

    if (associated(pf%cycles%start)) deallocate(pf%cycles%start)
    if (associated(pf%cycles%pfasst)) deallocate(pf%cycles%pfasst)
    if (associated(pf%cycles%end)) deallocate(pf%cycles%end)

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
          call F%encap%destroy(F%Q(m))
          do p = 1, size(F%F(m,:))
             call F%encap%destroy(F%F(m,p))
          end do
       end do

       do m = 1, F%nnodes-1
          call F%encap%destroy(F%S(m))
       end do

       if (F%Finterp) then
          do m = 1, F%nnodes
             do p = 1, size(F%F(m,:))
                if (associated(F%pF)) then
                   call F%encap%destroy(F%pF(m,p))
                end if
             end do
          end do
          if (associated(F%pQ)) then
             call F%encap%destroy(F%pQ(1))
          end if
       else
          do m = 1, F%nnodes
             if (associated(F%pQ)) then
                call F%encap%destroy(F%pQ(m))
             end if
          end do
       end if


       deallocate(F%Q)
       deallocate(F%F)
       deallocate(F%S)

       if (F%Finterp) then
          if (associated(F%pF)) deallocate(F%pF)
       else
          if (associated(F%pQ)) deallocate(F%pQ)
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
