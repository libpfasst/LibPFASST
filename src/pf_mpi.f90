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


module pf_mod_mpi
  include "mpif.h"
end module pf_mod_mpi

!> Module to implement communication routines in  MPI.
module pf_mod_comm_mpi
  use pf_mod_dtype
  use pf_mod_timer
  implicit none
contains

  !> Subroutine to create an MPI based PFASST communicator using the MPI communicator *mpi_comm*.
  subroutine pf_mpi_create(pf_comm, mpi_comm)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: mpi_comm

    integer :: ierror

    !> assign communicator
    pf_comm%comm = mpi_comm

    !> assign number of processors
    call mpi_comm_size(mpi_comm, pf_comm%nproc, ierror)

    !>  assign procedure pointers
    pf_comm%post => pf_mpi_post
    pf_comm%recv => pf_mpi_recv
    pf_comm%send => pf_mpi_send
    pf_comm%wait => pf_mpi_wait
    pf_comm%broadcast   => pf_mpi_broadcast
    pf_comm%recv_status => pf_mpi_recv_status
    pf_comm%send_status => pf_mpi_send_status
  end subroutine pf_mpi_create


  !> Subroutine to set up the PFASST communicator.
  !! This should be called soon after adding levels to the PFASST controller 
  subroutine pf_mpi_setup(pf_comm, pf,ierror)
    use pf_mod_mpi, only: MPI_REQUEST_NULL

    type(pf_comm_t),   intent(inout) :: pf_comm    !>  communicator 
    type(pf_pfasst_t), intent(inout) :: pf         !>  main pfasst structure
    integer,           intent(inout) :: ierror     !>  error flag

    !>  set the rank
    call mpi_comm_rank(pf_comm%comm, pf%rank, ierror)

    !>  allocate arrarys for and and receive requests
    allocate(pf_comm%recvreq(pf%nlevels))
    allocate(pf_comm%sendreq(pf%nlevels))
    allocate(pf_comm%recvreq_bwd(pf%nlevels))
    allocate(pf_comm%sendreq_bwd(pf%nlevels))

    
    pf_comm%sendreq = MPI_REQUEST_NULL
    pf_comm%sendreq_bwd = MPI_REQUEST_NULL
    pf_comm%statreq = -66
  end subroutine pf_mpi_setup


  !> Subroutine to destroy the PFASST communicator.
  subroutine pf_mpi_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm

    deallocate(pf_comm%recvreq)
    deallocate(pf_comm%sendreq)
    deallocate(pf_comm%recvreq_bwd)
    deallocate(pf_comm%sendreq_bwd)
  end subroutine pf_mpi_destroy

  !>  Subroutine to post receive requests.
  subroutine pf_mpi_post(pf, level, tag, ierror, direction)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(in   ) :: pf
    class(pf_level_t), intent(inout) :: level   !<  level to send from
    integer,           intent(in   ) :: tag     !<  message tag
    integer,           intent(inout) :: ierror  !<  error flag
    integer, optional, intent(in)    :: direction
    integer                          :: dir, source

    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction

!     if(dir==1) source = modulo(pf%rank-1, pf%comm%nproc)
!     if(dir==2) source = modulo(pf%rank+1, pf%comm%nproc)

    if (dir==2 ) then        
        call mpi_irecv(level%recv_bwd, level%mpibuflen, MPI_REAL8, &
                   modulo(pf%rank+1, pf%comm%nproc), tag, pf%comm%comm, pf%comm%recvreq_bwd(level%index), ierror)
    else
        call mpi_irecv(level%recv, level%mpibuflen, MPI_REAL8, &
                   modulo(pf%rank-1, pf%comm%nproc), tag, pf%comm%comm, pf%comm%recvreq(level%index), ierror)
    end if
  end subroutine pf_mpi_post


  !> Subroutine to send convergence status information
  subroutine pf_mpi_send_status(pf, tag,istatus,ierror, direction)
    use pf_mod_mpi, only: MPI_INTEGER4, MPI_STATUS_SIZE, MPI_REQUEST_NULL

    type(pf_pfasst_t), intent(inout) :: pf        !<  main pfasst structure
    integer,           intent(in)    :: tag       !<  message tag
    integer,           intent(in) :: istatus      !<  status flag to send
    integer,           intent(inout) :: ierror    !<  error flag
    integer, optional, intent(in)    :: direction
    integer                          :: dest
    integer    ::  stat(MPI_STATUS_SIZE), dir
    integer(4) :: message!(8)

    message = istatus
    
    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction
    
    if( dir == 2) then
       dest = modulo(pf%rank-1, pf%comm%nproc) 
    else
       dest = modulo(pf%rank+1, pf%comm%nproc) 
    end if


    if (pf%comm%statreq /= -66) then
       call mpi_wait(pf%comm%statreq, stat, ierror)
    end if

    call mpi_issend(message, 1, MPI_INTEGER4, &
                    dest, tag, pf%comm%comm, pf%comm%statreq, ierror)

  end subroutine pf_mpi_send_status

  !> Subroutine to receive convergence status information
  subroutine pf_mpi_recv_status(pf, tag,istatus,ierror, direction)
    use pf_mod_mpi, only: MPI_INTEGER4, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf        !<  main pfasst structure
    integer,           intent(in)    :: tag       !<  message tag
    integer,           intent(inout) :: istatus   !<  status flag to receive
    integer,           intent(inout) :: ierror    !<  error flag
    integer, optional, intent(in)    :: direction
    integer                          :: source
    integer    ::  stat(MPI_STATUS_SIZE), dir

    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction
    
    if(dir == 2) then
       source = modulo(pf%rank+1, pf%comm%nproc) 
    else
       source = modulo(pf%rank-1, pf%comm%nproc) 
    end if
    
    call mpi_recv(istatus, 1, MPI_INTEGER4, &
                  source, tag, pf%comm%comm, stat, ierror)

  end subroutine pf_mpi_recv_status


  !> Subroutine to send solutions
  subroutine pf_mpi_send(pf, level, tag, blocking,ierror, direction)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf       !<  main pfasst structure
    class(pf_level_t), intent(inout) :: level    !<  level to send from
    integer,           intent(in   ) :: tag      !<  message tag
    logical,           intent(in   ) :: blocking !<  true if send is blocking
    integer,           intent(inout) :: ierror   !<  error flag
    integer, optional, intent(in)    :: direction
    integer                          :: dest, dir
    integer ::  stat(MPI_STATUS_SIZE)

    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction
    
    if(dir == 2) then
       dest = modulo(pf%rank-1, pf%comm%nproc) 
    else
       dest = modulo(pf%rank+1, pf%comm%nproc) 
    end if
    
    if (blocking) then
      if(dir == 2) then
          call level%q0%pack(level%send_bwd, 2)
          call mpi_send(level%send_bwd, level%mpibuflen, MPI_REAL8, &
                     dest, tag, pf%comm%comm, stat, ierror)
       else
          call level%qend%pack(level%send, 1)
          call mpi_send(level%send, level%mpibuflen, MPI_REAL8, &
                     dest, tag, pf%comm%comm, stat, ierror)
       end if       
    else
       if(dir == 2) then
          call mpi_wait(pf%comm%sendreq_bwd(level%index), stat, ierror)
          call level%q0%pack(level%send_bwd, 2)
          call mpi_isend(level%send_bwd, level%mpibuflen, MPI_REAL8, &
                      dest, tag, pf%comm%comm, pf%comm%sendreq_bwd(level%index), ierror)
       else
          call mpi_wait(pf%comm%sendreq(level%index), stat, ierror)
          call level%qend%pack(level%send, 1)
          call mpi_isend(level%send, level%mpibuflen, MPI_REAL8, &
                      dest, tag, pf%comm%comm, pf%comm%sendreq(level%index), ierror)
       end if
    end if
  end subroutine pf_mpi_send

  !> Subroutine to receive solutions
  subroutine pf_mpi_recv(pf, level, tag, blocking, ierror, direction)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE
    type(pf_pfasst_t), intent(inout) :: pf     !<  main pfasst structure
    class(pf_level_t), intent(inout) :: level  !<  level to recieve into
    integer,           intent(in   ) :: tag    !<  message tag
    logical,           intent(in   ) :: blocking  !<  true if receive is blocking
    integer,           intent(inout) :: ierror  !<  error flag
    integer, optional, intent(in)    :: direction
    integer                          :: source, dir
    integer ::  stat(MPI_STATUS_SIZE)

    dir = 1 ! default 1: send forward; set to 2 for send backwards
    if(present(direction)) dir = direction
    
    if(dir == 2) then
       source = modulo(pf%rank+1, pf%comm%nproc) 
    else
       source = modulo(pf%rank-1, pf%comm%nproc) 
    end if
    
    if (blocking) then
      if (dir==2 ) then
        call mpi_recv(level%recv_bwd, level%mpibuflen, MPI_REAL8, &
                     source, tag, pf%comm%comm, stat, ierror)
      else
        call mpi_recv(level%recv, level%mpibuflen, MPI_REAL8, &
                     source, tag, pf%comm%comm, stat, ierror)
      end if
    else
      if (dir == 2) then
        call mpi_wait(pf%comm%recvreq_bwd(level%index), stat, ierror)
      else
        call mpi_wait(pf%comm%recvreq(level%index), stat, ierror)
      end if
    end if
  end subroutine pf_mpi_recv

  !
  ! Misc
  !
  subroutine pf_mpi_wait(pf, level, ierror)
    use pf_mod_mpi, only: MPI_STATUS_SIZE
    type(pf_pfasst_t), intent(in) :: pf           !<  main pfasst structure
    integer,           intent(in) :: level        !<  level on which to wait
    integer,           intent(inout) :: ierror    !<  error flag
    integer ::  stat(MPI_STATUS_SIZE)
    call mpi_wait(pf%comm%sendreq(level), stat, ierror)
  end subroutine pf_mpi_wait

  subroutine pf_mpi_broadcast(pf, y, nvar, root,ierror)
    use pf_mod_mpi, only: MPI_REAL8
    type(pf_pfasst_t), intent(inout) :: pf      !<  main pfasst structure
    real(pfdp),        intent(in)    :: y(nvar) !<  data to broadcast
    integer,           intent(in)    :: nvar    !<  size of data to broadcast
    integer,           intent(in)    :: root    !<  rank of broadcaster
    integer,           intent(inout) :: ierror  !<  error flag
    call mpi_bcast(y, nvar, MPI_REAL8, root, pf%comm%comm, ierror)
  end subroutine pf_mpi_broadcast

end module pf_mod_comm_mpi
