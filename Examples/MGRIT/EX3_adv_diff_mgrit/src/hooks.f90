!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pfasst
  use pf_mod_ndarray
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    use my_sweeper, only: exact
    use probin, only: ndim, solver_type, Tfin
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp), allocatable :: yexact_1d(:), yexact_2d(:,:), yexact_3d(:,:,:)
    real(pfdp), pointer :: y_end_1d(:), y_end_2d(:,:), y_end_3d(:,:,:)
    real(pfdp) :: error
    real(pfdp) :: residual
    real(pfdp) :: t

    !if ((solver_type .eq. 1) .or. (solver_type .eq. 4)) then
    !   t = Tfin
    !else
    !   t = Tfin
    !end if

    t = pf%state%t0+pf%state%dt
   
    !>  compute the error at last end point
    if (ndim == 1) then
       allocate(yexact_1d(pf%levels(level_index)%lev_shape(1)))
       y_end_1d => get_array1d(pf%levels(level_index)%qend)
       call exact(t, yexact_1d)
       error = maxval(abs(y_end_1d-yexact_1d))
       deallocate(yexact_1d)
    else if (ndim == 2) then
       allocate(yexact_2d(pf%levels(level_index)%lev_shape(1), &
                          pf%levels(level_index)%lev_shape(2)))
       y_end_2d => get_array2d(pf%levels(level_index)%qend)
       call exact(t, yexact_2d)
       error = maxval(abs(y_end_2d-yexact_2d))
       deallocate(yexact_2d)
    else if (ndim == 3) then
       allocate(yexact_3d(pf%levels(level_index)%lev_shape(1), &
                          pf%levels(level_index)%lev_shape(2), &
                          pf%levels(level_index)%lev_shape(3)))
       y_end_3d => get_array3d(pf%levels(level_index)%qend)
       call exact(t, yexact_3d)
       error = maxval(abs(y_end_3d-yexact_3d))
       deallocate(yexact_3d)
    else
    end if

    residual=pf%levels(level_index)%residual


    !print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
    !     pf%state%step+1, pf%state%iter,level_index, error,residual

    !if (level_index .eq. pf%nlevels) then
    !print '("error: time: ", f8.4," step: ",i8.1," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," resid: ",es14.7," deltaq0: ",es14.7)', &
    !     t,pf%state%step+1, pf%rank, pf%state%iter,level_index,error,pf%levels(level_index)%residual,pf%levels(level_index)%max_delta_q0    
    !end if

    !if ((pf%state%step .eq. pf%state%nsteps-1) .and. (level_index == pf%nlevels)) then
    if ((pf%state%step .eq. pf%state%nsteps-1) .and. &
        (level_index == pf%nlevels) .and. &
        (pf%state%iter .gt. 0) .and. &
        (pf%state%sweep .eq. pf%levels(level_index)%nsweeps)) then
       print '("rank: ", i4.4," step: ",i4.4," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.9e4)', &
            pf%rank,pf%state%step+1, pf%state%iter,level_index, error, residual
       call flush(6)
    end if

    call pf_set_error(pf, level_index, error)
    call pf_set_resid(pf, level_index, residual)
    
  end subroutine echo_error
end module hooks
