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
    integer :: dim(2),j
    real(pfdp) :: error,t_end,errd

    real(pfdp)   :: q_ex(20)
    real(pfdp)   :: q_diag(20)
    
    
    !  20 particles run with 512
    q_ex =(/ 1.49416313027551D+00, 7.07931923619350D-01, 9.18515291508953D-01, 7.82054090221508D-01, 2.69391225356707D-01, 2.07255594164204D-01, 2.88290898242198D-01, &
         3.95894491436281D-01, 4.87366427184367D-01, 4.08042164173400D-01, 4.08042164173065D-01, 4.87366427184291D-01, 3.95894491436375D-01, 2.88290898242375D-01, &
         2.07255594164250D-01, 2.69391225356918D-01, 7.82054090222874D-01, 9.18515291509728D-01, 7.07931923618189D-01, 1.49416313027296D+00 /)
    
    t_end=pf%state%t0+pf%state%dt
    
    imk => cast_as_imk_sweeper(pf%levels(level_index)%ulevel%sweeper)
    qend => cast_as_zmkpair(pf%levels(level_index)%qend)
    dim = shape(qend%array)
    call zmkpair_build(yexact, dim(1))

!    yexact%array = 0.0_pfdp
!    call exact(yexact, pf%state%t0+pf%state%dt)
!    error = compute_inf_norm(qend%array - yexact%array, dim(1))
!    error = error / compute_inf_norm(yexact%array, dim(1))
!    print '("time: ", f5.2, " step: ",i3.3," iter: ",i4.3," error (dmat): ",es14.7)', &
!         pf%state%t0+pf%state%dt, pf%state%step+1, pf%state%iter, error
    if (abs(t_end-0.05_pfdp) .lt. 1e-12) then
       do j=1,20
          q_diag(j) = real(qend%y(j,j))
       end do
       print *,q_diag
       print *,q_ex
       errd = maxval(abs(q_diag-q_ex))
       
       print *,'Rank ',pf%rank,'Nsteps ',pf%state%step+1,' Iter ',pf%state%iter,' error at end=',errd, 'Nnodes=',pf%nnodes(level_index), 'qtype=',pf%qtype
    endif
    
    deallocate(yexact%array)

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
