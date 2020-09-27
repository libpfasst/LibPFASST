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
  subroutine initialize_results(pf)
    type(pf_pfasst_t), intent(inout) :: pf    

    character(len = 128) :: datpath  !!  path to output files
    character(len = 128) :: dirname  !!  output file name for residuals
    integer :: istat,system
    integer nlevs,nsteps,nblocks,nprocs,niters,max_nsweeps    !  local variables for code brevity
    nlevs=pf%nlevels
    nblocks=pf%nlevels
    nsteps=pf%state%nsteps
    nprocs=pf%comm%nproc
    niters=pf%niters
    nblocks=nsteps/nprocs
    max_nsweeps=max(maxval(pf%nsweeps),maxval(pf%nsweeps_pred))
    

    !  Load up the results structure
    pf%results%nlevs=nlevs
    pf%results%nsteps=nsteps
    pf%results%nprocs=nprocs
    pf%results%niters=niters
    pf%results%nblocks=nblocks
    pf%results%max_nsweeps=max_nsweeps
    pf%results%rank=pf%rank
    pf%results%save_residuals=pf%save_residuals
    pf%results%save_errors=pf%save_errors
    pf%results%save_delta_q0=pf%save_delta_q0
    istat=0
    if(.not.allocated(pf%results%errors) .and. pf%results%save_errors) then
       allocate(pf%results%errors(nlevs, nblocks,niters+1, max_nsweeps),stat=istat)
       if (istat /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',istat)
       pf%results%errors = -1.0_pfdp
    end if
    
    if(.not.allocated(pf%results%residuals) .and. pf%results%save_residuals) then
       allocate(pf%results%residuals(nlevs, nblocks,niters+1, max_nsweeps),stat=istat)
       if (istat /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',istat)
       pf%results%residuals = -1.0_pfdp
    endif
    

    if(.not.allocated(pf%results%delta_q0) .and. pf%results%save_delta_q0) then
       allocate(pf%results%delta_q0(nlevs,nblocks,niters+1, max_nsweeps),stat=istat)
       if (istat /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',istat)
       pf%results%delta_q0 = -1.0_pfdp
    end if
    
    if(.not.allocated(pf%results%iters)) allocate(pf%results%iters(nblocks),stat=istat)
    if (istat /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',istat)                   
    pf%results%iters = niters  

    datpath= 'dat/' // trim(pf%outdir) // '/'

    ! Create directory for this processor
    write (dirname, "(A5,I0.4)") 'Proc_',pf%results%rank
    datpath=trim(datpath) // trim(dirname)

    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make Proc directory")
    !  Final path for all the stat files
    pf%results%datpath= trim(datpath)

  end subroutine initialize_results

  subroutine dump_results(this)
    type(pf_results_t), intent(inout) :: this
    integer :: kblock, kiter, ksweep,nstep,klevel
    integer :: system
    integer :: estream,rstream,qstream,istream
    character(len = 256) :: iname  !!  output file name for residuals
    character(len = 256) :: resname  !!  output file name for residuals
    character(len = 256) :: errname  !!  output file name for errors
    character(len = 256) :: q0name   !!  output file name for delta_q0


    !  open files for output
    if (this%save_residuals) then
       resname = trim(this%datpath) // '/residual.dat'
       rstream=5000+this%rank
       open(rstream, file=trim(resname), form='formatted',err=999)
    end if
    if (this%save_errors) then
       errname = trim(this%datpath) // '/error.dat'
       estream=6000+this%rank
       open(estream, file=trim(errname), form='formatted',err=999)
    end if
    if (this%save_delta_q0) then
       q0name = trim(this%datpath) // '/delta_q0.dat'
       qstream=7000+this%rank
       open(qstream, file=trim(q0name), form='formatted',err=999)
    end if
    if (this%save_residuals .or. this%save_errors .or. this%save_delta_q0) then
       if (this%save_residuals) write(rstream,*)'#level   step     block  iter  sweep   residual'
       if (this%save_errors) write(estream,*)'#level   step     block  iter  sweep   error'
       if (this%save_delta_q0) write(qstream,*)'#level   step     block  iter  sweep   delta q0'
       do klevel=1,this%nlevs
          do kblock = 1, this%nblocks
             nstep=(kblock-1)*this%nprocs+this%rank+1
             do kiter = 0 , this%niters
                do ksweep = 1, this%max_nsweeps
                   if (this%save_residuals) write(rstream,101 ) klevel,nstep,kblock,kiter,ksweep,this%residuals(klevel, kblock,kiter+1, ksweep)
                   if (this%save_errors) write(estream,101) klevel,nstep,kblock,kiter,ksweep,this%errors(klevel,kblock,kiter+1,  ksweep)
                   if (this%save_delta_q0) write(qstream,101) klevel,nstep,kblock,kiter,ksweep,this%delta_q0(klevel, kblock,kiter+1, ksweep)
                end do
             end do
          enddo
       enddo
101    format(I3,I10, I10,I6, I6, e22.14)
       if (this%save_residuals) close(rstream)
       if (this%save_errors) close(estream)
       if (this%save_delta_q0) close(qstream)
    end if

    !  output file for iters
    iname = trim(this%datpath) // '/iter.dat'
    istream=8000+this%rank
    open(istream, file=trim(iname), form='formatted',err=999)
    do kblock = 1, this%nblocks
       nstep=(kblock-1)*this%nprocs+this%rank+1
       write(istream, '(I10,I10, e22.14)') nstep,kblock,this%iters(kblock)
    enddo
    close(istream)
    return
    999 call pf_stop(__FILE__,__LINE__, "Error opening file")    
  end subroutine dump_results

  subroutine dump_timingsl(this,pf)
    type(pf_results_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    character(len = 128   ) :: pname     !!  processor name
    character(len = 256   ) :: fullname  !!  output file name for runtimes
    character(len = 128   ) :: datpath  !!  directory path
    character(len = 128   ) :: strng      !  used for string conversion
    integer :: j, iout,system,nlev,k,kmax
    real(pfdp) :: qarr(pf%nlevels)

    !  Write a json file with timer numbers and times
    fullname = trim(this%datpath) // '/runtime.json'
    iout = 4000+pf%rank !  Use processor dependent file number
    nlev=pf%nlevels
    !  output timings
    open(iout, file=trim(fullname), form='formatted',err=998)
    write(iout,*) '{'
    kmax=1
    if (pf%save_timings > 1) kmax=4

    do k=1,kmax
       write(iout,"(A24,A1,e14.6,A1)") wrap_timer_name(timer_names(k)), ': ', pf%pf_timers%runtimes(k,1), ','
    end do
    if (pf%save_timings > 1) then
       do k=kmax+1,PF_NUM_TIMERS
          if (nlev .eq. 1) then
             write(iout,"(A24,A1,e14.6,A1)") wrap_timer_name(timer_names(k)), ':', pf%pf_timers%runtimes(k,1), ','
          else
             qarr=pf%pf_timers%runtimes(k,1:nlev)
             strng=trim(convert_real_array(qarr,nlev))       
             write(iout,"(A24,A1,A60,A1)")  wrap_timer_name(timer_names(k)),':', adjustl(strng), ','
          end if
       end do
    end if
    write(iout,*) '"foo":"0"'
    write(iout,*) '}'
    
    close(iout)
    return
998 call pf_stop(__FILE__,__LINE__, "Error opening json file")    

  end subroutine dump_timingsl

  
  subroutine destroy_results(this)
    type(pf_results_t), intent(inout) :: this
    integer :: istat
    istat=0
    if(allocated(this%errors))  deallocate(this%errors,stat=istat)
    if(allocated(this%residuals))  deallocate(this%residuals,stat=istat)
    if(allocated(this%delta_q0))  deallocate(this%delta_q0,stat=istat)
    if(allocated(this%iters))  deallocate(this%iters,stat=istat)
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Error deallocating arrays")        
  end subroutine destroy_results

end module pf_mod_results
