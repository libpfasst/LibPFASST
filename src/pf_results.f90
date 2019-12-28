!! Module for storing results for eventual output
!
! This file is part of LIBPFASST.
!
!>  Module for the storing results for eventual output
module pf_mod_results
  use pf_mod_dtype
  use pf_mod_utils
  implicit none
  
contains
  subroutine initialize_results(this, nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,level_index,outdir,save_residuals)

    class(pf_results_t), intent(inout) :: this
    integer, intent(in) :: nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,level_index
    character(len=*), intent(in) :: outdir
    logical, intent(in) :: save_residuals

    character(len = 128) :: fname  !!  output file name for residuals
    character(len = 128) :: datpath  !!  path to output files
    character(len = 256) :: fullname  !!  output file name for residuals
    integer :: istat,system,ierr
    
    !  Set up the directory to dump results
    istat= system('mkdir -p dat')
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in initialize_results")       
    istat= system('mkdir -p dat/' // trim(outdir))       
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in initialize_results")
    this%datpath= 'dat/' // trim(outdir) // '/'
!    this%datpath=this%datpath //   '/'
    
    if (save_residuals) then
       
       write (fname, "(A16,I0.1,A4)") 'residuals_size_L',level_index,'.dat'
       fullname = trim(this%datpath) // trim(fname)
       
       if (rank_in == 0) then
          open(unit=123, file=trim(fullname), form='formatted')
          write(123,'(I5, I5, I5, I5)') nsteps_in, niters_in, nprocs_in, nsweeps_in
          close(unit=123)
       end if
    end if
       
!    this%dump => dump_results
    this%destroy => destroy_results

    this%nsteps=nsteps_in
    this%nblocks=nsteps_in/nprocs_in
    this%niters=niters_in
    this%nprocs=nprocs_in
    this%nsweeps=nsweeps_in
    this%rank=rank_in
    this%level_index=level_index    

    ierr=0
    if(.not.allocated(this%errors)) allocate(this%errors(niters_in, this%nblocks, nsweeps_in),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)               
    if(.not.allocated(this%residuals)) allocate(this%residuals(niters_in, this%nblocks, nsweeps_in),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                   
    if(.not.allocated(this%delta_q0)) allocate(this%delta_q0(niters_in, this%nblocks, nsweeps_in),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                   

    this%errors = -1.0_pfdp
    this%residuals = -1.0_pfdp
    this%delta_q0 = -1.0_pfdp
  end subroutine initialize_results

  subroutine dump_resids(this)
    type(pf_results_t), intent(inout) :: this
    integer :: i, j, k, istat,system,istream,nstep
    character(len = 128) :: fname  !!  output file name for residuals
    character(len = 256) :: fullname  !!  output file name for residuals
    character(len = 128) :: datpath  !!  directory path
    character(len = 128) :: dirname  !!  directory name

    ! Create directory for residual output if necessary
    datpath = trim(this%datpath) // 'residuals'
    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_resids")
    ! Create directory for this processor
    write (dirname, "(A6,I0.3)") '/Proc_',this%rank
    datpath=trim(datpath) // trim(dirname) 
    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_resids")

    !  output residuals for each sweep and block and iteration
    write (fname, "(A5,I0.1,A4)") '/Lev_',this%level_index,'.dat'
    fullname = trim(datpath) // trim(fname)
    istream=5000+this%rank
    open(istream, file=trim(fullname), form='formatted')

    do j = 1, this%nblocks
       nstep=(j-1)*this%nprocs+this%rank+1
       do i = 1 , this%niters
          do k = 1, this%nsweeps
             if (this%residuals(i, j, k) .ge. 0.0) then
                write(istream, '(I5, I4,I4, I4, e22.14)') nstep,j,i,k,this%residuals(i, j, k)
             end if
          end do
       end do
    enddo
    close(istream)

    !  output residuals only for last iteration
    write (fname, "(A5,I0.1,A9)") '/Lev_',this%level_index,'_iter.dat'
    fullname = trim(datpath) // trim(fname)
    istream=5000+this%rank
    open(istream, file=trim(fullname), form='formatted')

    do j = 1, this%nblocks
       nstep=(j-1)*this%nprocs+this%rank+1
       do i = this%niters,1,-1
          if (this%residuals(i, j, this%nsweeps) .ge. 0.0) then
             write(istream, '(I5,I4, I4, e22.14)') nstep,j,i,this%residuals(i, j, this%nsweeps)
             exit
          end if
       end do
    enddo
    close(istream)
    
  end subroutine dump_resids
  subroutine dump_delta_q0(this)
    type(pf_results_t), intent(inout) :: this
    integer :: i, j, k, istat,system,istream,nstep
    character(len = 128) :: fname  !!  output file name for residuals
    character(len = 256) :: fullname  !!  output file name for residuals
    character(len = 128) :: datpath  !!  directory path
    character(len = 128) :: dirname  !!  directory name

    ! Create directory for residual output if necessary
    datpath = trim(this%datpath) // 'delta_q0'
    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_delta_q0")
    ! Create directory for this processor
    write (dirname, "(A6,I0.3)") '/Proc_',this%rank
    datpath=trim(datpath) // trim(dirname) 
    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_delta_q0")

    !  output delta_q0 for each sweep and block and iteration
    write (fname, "(A5,I0.1,A4)") '/Lev_',this%level_index,'.dat'
    fullname = trim(datpath) // trim(fname)
    istream=5000+this%rank
    open(istream, file=trim(fullname), form='formatted')

    do j = 1, this%nblocks
       nstep=(j-1)*this%nprocs+this%rank+1
       do i = 1 , this%niters
          do k = 1, this%nsweeps
             if (this%delta_q0(i, j, k) .ge. 0.0) then
                write(istream, '(I5, I4,I4, I4, e22.14)') nstep,j,i,k,this%delta_q0(i, j, k)
             end if
          end do
       end do
    enddo
    close(istream)

    !  output delta_q0 only for last iteration
    write (fname, "(A5,I0.1,A9)") '/Lev_',this%level_index,'_iter.dat'
    fullname = trim(datpath) // trim(fname)
    istream=5000+this%rank
    open(istream, file=trim(fullname), form='formatted')

    do j = 1, this%nblocks
       nstep=(j-1)*this%nprocs+this%rank+1
       do i = this%niters,1,-1
          if (this%delta_q0(i, j, this%nsweeps) .ge. 0.0) then
             write(istream, '(I5,I4, I4, e22.14)') nstep,j,i,this%delta_q0(i, j, this%nsweeps)
             exit
          end if
       end do
    enddo
    close(istream)
    
  end subroutine dump_delta_q0
  subroutine dump_errors(this)
    type(pf_results_t), intent(inout) :: this

    integer :: i, j, k, istat,system,istream,nstep
    character(len = 128   ) :: fname  !!  output file name for residuals
    character(len = 256   ) :: fullname  !!  output file name for residuals
    character(len = 128   ) :: datpath  !!  directory path
    character(len = 128) :: dirname  !!  directory name
    
    !  Create directory for output if necessary
    datpath = trim(this%datpath) // 'errors'
    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_errors")

    ! Create directory for this processor
    write (dirname, "(A6,I0.3)") '/Proc_',this%rank
    datpath=trim(datpath) // trim(dirname) 
    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_errors")    

    !  Build file name for output per sweep
    write (fname, "(A5,I0.1,A4)") '/Lev_',this%level_index,'.dat'
    fullname = trim(datpath) // trim(fname)               !  full path to file
    istream=2000+this%rank                                !  define unit number
    open(istream, file=trim(fullname), form='formatted')  !  open file

    !  output errors  per sweep
    do j = 1, this%nblocks
       nstep=(j-1)*this%nprocs+this%rank+1
       do i = 1 , this%niters
          do k = 1, this%nsweeps
             if (this%errors(i, j, k) .ge. 0.0) then
                write(istream, '(I5, I4, I4, I4, e22.14)') nstep,j,i,k,this%errors(i, j, k)
             end if
          end do
       end do
    enddo
    close(istream)

    !  Build file name for output at end of step
    write (fname, "(A5,I0.1,A9)") '/Lev_',this%level_index,'_iter.dat'
    fullname = trim(datpath) // trim(fname)               !  full path to file
    istream=2000+this%rank                               !  define unit number
    open(istream, file=trim(fullname), form='formatted')  !  open file

    !  output errors  per sweep
    do j = 1, this%nblocks
       nstep=(j-1)*this%nprocs+this%rank+1
       do i = this%niters,1,-1
          if (this%errors(i, j, this%nsweeps) .ge. 0.0) then
             write(istream, '(I5,I4, I4, e22.14)') nstep,j,i,this%errors(i, j, this%nsweeps)
             exit
          end if
       end do
    enddo
    flush(istream)
    close(istream)

  end subroutine dump_errors

  subroutine dump_timings(this,pf)
    type(pf_results_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    character(len = 128   ) :: pname     !!  processor name
    character(len = 256   ) :: fullname  !!  output file name for runtimes
    character(len = 128   ) :: datpath  !!  directory path
    character(len = 128   ) :: strng      !  used for string conversion
    integer :: istat,j, iout,system,nlev

    datpath = trim(this%datpath) // '/runtimes/'
    istat= system('mkdir -p '// trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_timings")

    ! Create directory for this processor
    write (pname, "(A5,I0.3)") 'Proc_',this%rank
    datpath=trim(datpath) // trim(pname) 

    !  Write a json file with timer numbers and times
    fullname = trim(datpath) // '_times.json'
    iout = 4000+pf%rank !  Use processor dependent file number
    nlev=pf%nlevels
    !  output timings
    open(iout, file=trim(fullname), form='formatted')
    write(iout,*) '{'
    write(iout,"(A24,e12.6,A1)")  '"total" :',       pf%pf_runtimes%t_total, ','
    if (pf%save_timings > 1) then
       write(iout,"(A24,e12.6,A1)")  '"predictor" :',   pf%pf_runtimes%t_predictor, ','
       write(iout,"(A24,e12.6,A1)")  '"iteration" :',   pf%pf_runtimes%t_iteration, ','
       write(iout,"(A24,e12.6,A1)")  '"broadcast" :',   pf%pf_runtimes%t_broadcast, ','
       write(iout,"(A24,e12.6,A1)")  '"step" :',        pf%pf_runtimes%t_step, ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_hooks(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"hooks" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_sweeps(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"sweeps" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_interpolate(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"interpolate" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_restrict(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"restrict" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_send(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"wait" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_wait(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"send" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_receive(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"receive" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_residual(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"residual" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_feval(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"feval" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_fcomp(1:nlev),nlev))
       write(iout,"(A24,A60,A1)")  '"fcomp" :', adjustl(strng), ','
       strng=trim(convert_real_array(pf%pf_runtimes%t_aux(1:nlev),nlev))
       write(iout,"(A24,A60)")  '"aux" :', adjustl(strng)
    end if
    
    write(iout,*) '}'
    
    close(iout)
    
  end subroutine dump_timings

  subroutine destroy_results(this)
    type(pf_results_t), intent(inout) :: this
    
    if(allocated(this%errors))  deallocate(this%errors)
    if(allocated(this%residuals))  deallocate(this%residuals)
    if(allocated(this%delta_q0))  deallocate(this%delta_q0)
  end subroutine destroy_results

end module pf_mod_results
