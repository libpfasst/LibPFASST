module hooks
  use pf_mod_dtype
  use factory
  use sweeper

  implicit none

contains

  subroutine echo_error(pf, level, state)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    class(imk_sweeper_t), pointer :: imk
    type(zndarray) :: yexact
    type(zndarray), pointer :: qend
    integer :: dim(2)
    real(pfdp) :: error

    imk => cast_as_imk_sweeper(level%ulevel%sweeper)
    qend => cast_as_zndarray(level%qend)
    dim = shape(qend%array)
    call zndarray_build(yexact, dim(1))

    yexact%array = 0.0_pfdp
    call exact(yexact, state%t0+state%dt)
    error = compute_inf_norm(qend%array - yexact%array, dim(1))
    error = error / compute_inf_norm(yexact%array, dim(1))
    print '("time: ", f5.2, " step: ",i3.3," iter: ",i4.3," error (dmat): ",es14.7)', &
         state%t0+state%dt, state%step+1, state%iter, error

    ! stop
    ! print*, '------------------------------------------------------'
    deallocate(yexact%array)

    !if (state%iter == 2) stop
  end subroutine echo_error

  subroutine echo_residual(pf, level, state)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    type(zndarray), pointer :: r, q

    r => cast_as_zndarray(level%R(level%nnodes-1))
    q => cast_as_zndarray(level%Q(level%nnodes-1))

    print '("resid: time: ", f10.5," rank: ",i3.3," step: ",i5.5," iter: ",i4.3," level: ",i1.1," resid: ",es18.10e4)', &
         state%t0+state%dt, pf%rank, state%step+1, state%iter, level%index, maxval(abs(r%array))

    if (any(isnan(real(q%array)))) then
       print*, 'omega is NAN', q%array
       stop
    elseif (any(isnan(real(q%y)))) then
       print*, 'y is NAN', q%y
       stop
    endif
  end subroutine echo_residual


  subroutine save_solution(pf, level, state)
    use probin, only: fbase
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in   ) :: state

    type(zndarray), pointer :: qend
    character(len=256) :: time, filename
    integer :: N,i,j,un

    qend => cast_as_zndarray(level%qend)
    un = 200+pf%rank
    write(time, '(f10.5)') state%t0+state%dt
    write(filename, '("-rank_", i3.3, "-step_",i5.5,"-iter_",i3.3,"-level_",i1.1,"_soln")') &
         pf%rank, state%step+1, state%iter, level%index
    open(unit=un, file=trim(fbase)//'/time_'//trim(adjustl(time))//trim(adjustl(filename)), form='unformatted')

    write(un) qend%y
    
    close(un)

  end subroutine save_solution

end module hooks
