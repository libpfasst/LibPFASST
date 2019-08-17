module hooks
  use pf_mod_dtype
  use mod_zmkpair
  use sweeper

  implicit none

contains

  subroutine echo_error(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(imk_sweeper_t), pointer :: imk
    type(zmkpair) :: yexact
    type(zmkpair), pointer :: qend
    integer :: dim(2)
    real(pfdp) :: error

    
    imk => cast_as_imk_sweeper(pf%levels(level_index)%ulevel%sweeper)
    qend => cast_as_zmkpair(pf%levels(level_index)%qend)
    dim = shape(qend%array)
    call zmkpair_build(yexact, dim(1))

    yexact%array = 0.0_pfdp
    call exact(yexact, pf%state%t0+pf%state%dt)
    error = compute_inf_norm(qend%array - yexact%array, dim(1))
    error = error / compute_inf_norm(yexact%array, dim(1))
    print '("time: ", f5.2, " step: ",i3.3," iter: ",i4.3," error (dmat): ",es14.7)', &
         pf%state%t0+pf%state%dt, pf%state%step+1, pf%state%iter, error

    ! stop
    ! print*, '------------------------------------------------------'
    deallocate(yexact%array)

    !if (pf%state%iter == 2) stop
  end subroutine echo_error

  subroutine echo_residual(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    type(zmkpair), pointer :: r, q

    r => cast_as_zmkpair(pf%levels(level_index)%R(pf%levels(level_index)%nnodes-1))
    q => cast_as_zmkpair(pf%levels(level_index)%Q(pf%levels(level_index)%nnodes-1))

    print '("resid: time: ", f10.5," rank: ",i3.3," step: ",i5.5," iter: ",i4.3," level: ",i1.1," resid: ",es18.10e4)', &
         pf%state%t0+pf%state%dt, pf%rank, pf%state%step+1, pf%state%iter, level_index, maxval(abs(r%array))

    if (any(isnan(real(q%array)))) then
       print*, 'omega is NAN', q%array
       stop
    elseif (any(isnan(real(q%y)))) then
       print*, 'y is NAN', q%y
       stop
    endif
  end subroutine echo_residual


  subroutine save_solution(pf, level_index)
    use probin, only: fbase
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    type(zmkpair), pointer :: qend
    character(len=256) :: time, filename
    integer :: N,i,j,un

    qend => cast_as_zmkpair(pf%levels(level_index)%qend)
    un = 200+pf%rank
    write(time, '(f10.5)') pf%state%t0+pf%state%dt
    write(filename, '("-rank_", i3.3, "-step_",i5.5,"-iter_",i3.3,"-level_",i1.1,"_soln")') &
         pf%rank, pf%state%step+1, pf%state%iter, level_index
    open(unit=un, file=trim(fbase)//'/time_'//trim(adjustl(time))//trim(adjustl(filename)), form='unformatted')

    write(un) qend%y
    
    close(un)

  end subroutine save_solution

end module hooks
