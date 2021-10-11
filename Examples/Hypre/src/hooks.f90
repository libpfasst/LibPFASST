!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use encap
  use pf_my_sweeper
  use probin
  implicit none

  interface
     function HypreMaxErr(x, t, init_cond) result(max_err) bind(c, name="HypreMaxErr")
        use iso_c_binding
        type(c_ptr), value :: x
        real(c_double), value :: t
        real(c_double), value :: init_cond
        real(c_double) :: max_err
     end function
  end interface
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: yexact
    real(pfdp) :: maxerr, error 
    real(pfdp) :: residual
    class(hypre_vector_encap), pointer :: y_end
    integer :: nproc, rank, ierr, step

    !> Get the solution at the end of this step
    y_end => cast_as_hypre_vector(pf%levels(level_index)%qend)

    !>  compute error
    error = HypreMaxErr(y_end%c_hypre_vector_ptr, Tfin, init_cond)
    residual = pf%levels(level_index)%residual

    if (solver_type .eq. 1) then
       pf%results%residuals(level_index,1,pf%state%iter+1,1) = residual
       pf%results%errors(level_index,1,pf%state%iter+1,1) = error
    else
       pf%results%residuals(level_index,pf%state%pfblock,pf%state%iter+1,pf%state%sweep) = residual
       pf%results%errors(level_index,pf%state%pfblock,pf%state%iter+1,pf%state%sweep) = error
    end if
   
    !call mpi_comm_rank(pf%comm%comm, rank, ierr)
    !call mpi_comm_size(pf%comm%comm, nproc, ierr)

    if ((pf%state%step .eq. pf%state%nsteps-1) .and. (level_index == pf%nlevels) .and. (pf%state%iter .gt. 0)) then
       print '("error: rank: ", i4.4," step: ",i4.4," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.9e4)', &
            pf%rank,pf%state%step+1, pf%state%iter,level_index, error, residual
       call flush(6)
    end if
  end subroutine echo_error


end module hooks
