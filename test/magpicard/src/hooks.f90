module hooks
  use pf_mod_dtype
  use pf_mod_zndarray
  use sweeper
  implicit none

contains

  subroutine echo_error(pf, level_index)
    use probin, only:  magnus_order, nnodes
    type(pf_pfasst_t), intent(inout) :: pf
    integer,  intent(in   ) :: level_index

    class(magpicard_sweeper_t), pointer :: magpicard
    type(zndarray) :: yexact
    type(zndarray), pointer :: qend
    integer :: dim(2)
    real(pfdp) :: errd,t_end
    complex(pfdp),      pointer :: q_array(:,:)
    real(pfdp)   :: q_ex(20)
    real(pfdp)   :: q_diag(20)

    
   !  20 particles run with 512
   q_ex =(/ 1.49416313027551D+00, 7.07931923619350D-01, 9.18515291508953D-01, 7.82054090221508D-01, 2.69391225356707D-01, 2.07255594164204D-01, 2.88290898242198D-01, &
            3.95894491436281D-01, 4.87366427184367D-01, 4.08042164173400D-01, 4.08042164173065D-01, 4.87366427184291D-01, 3.95894491436375D-01, 2.88290898242375D-01, &
            2.07255594164250D-01, 2.69391225356918D-01, 7.82054090222874D-01, 9.18515291509728D-01, 7.07931923618189D-01, 1.49416313027296D+00 /)
    magpicard => cast_as_magpicard_sweeper(pf%levels(level_index)%ulevel%sweeper)
    qend => cast_as_zndarray(pf%levels(level_index)%qend)
    q_array=>get_array2d(pf%levels(level_index)%qend)

    t_end=pf%state%t0+pf%state%dt

    if (abs(t_end-0.05_pfdp) .lt. 1e-12) then
       q_diag = real(qend%flatarray(1:400:21))
       errd = maxval(abs(q_diag-q_ex))
      
       print *,'Rank ',pf%rank,'Nsteps ',pf%state%step+1,' Iter ',pf%state%iter,' error at end=',errd, 'totTime=',pf%runtimes(1),'Ord=',magnus_order, 'Nnodes=',nnodes(1), 'qtype=',pf%qtype
    endif


!!$    print '("time: ", f5.2, " step: ",i3.3," iter: ",i4.3," error (dmat): ",es14.7)', &
!!$         state%t0+state%dt, state%step+1, state%iter, error

  end subroutine echo_error



  subroutine save_solution(pf, level_index)
    use probin, only: fbase
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    type(zndarray), pointer :: qend,Fend
    character(len=256) :: time, filename
    integer :: un,istat,system
    complex(pfdp),      pointer :: q_array(:,:)

    !  Solution at the end

    qend => cast_as_zndarray(pf%levels(level_index)%qend)
    q_array=>get_array2d(pf%levels(level_index)%qend)
    

    istat= system('mkdir -p ' // trim(fbase))
    un = 200+pf%rank
    write(time, '(f10.5)') pf%state%t0+pf%state%dt
    write(filename, '("-rank_", i3.3, "-step_",i5.5,"-iter_",i3.3,"-level_",i1.1,"_soln")') &
         pf%rank, pf%state%step+1, pf%state%iter, level_index
    open(unit=un, file=trim(fbase)//'/time_'//trim(adjustl(time))//trim(adjustl(filename)), form='unformatted')
    write(un) q_array
    
    close(un)
    
    write(filename, '("P_step_",i5.5)') &
          pf%state%step+1

    open(unit=20, file=trim(fbase)//'/'//trim(adjustl(filename)), form='formatted')    
    write(20,*) real(qend%flatarray)
    write(20,*) ' '
    write(20,*) aimag(qend%flatarray)    
    close(20)


    Fend => cast_as_zndarray(pf%levels(level_index)%F(pf%levels(level_index)%nnodes,1))
    write(filename, '("F_step_",i5.5)') &
          pf%state%step+1
    open(unit=21, file=trim(fbase)//'/'//trim(adjustl(filename)), form='formatted')    
    write(21,*) real(Fend%flatarray)
    write(21,*) ' '
    write(21,*) aimag(Fend%flatarray)    
    close(21)
 
  end subroutine save_solution

end module hooks
