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
    integer :: istat,system,ierr

    pf%results%nlevs=pf%nlevels
    pf%results%nsteps=pf%state%nsteps
    pf%results%nprocs=pf%comm%nproc
    pf%results%niters=pf%niters
    pf%results%nsweeps=pf%nsweeps
    pf%results%rank=pf%rank
    pf%results%nblocks=pf%results%nsteps/pf%results%nprocs
    pf%results%max_nsweeps=maxval(pf%results%nsweeps)

    ierr=0
    if(.not.allocated(pf%results%errors)) allocate(pf%results%errors(pf%results%nlevs, pf%results%nblocks,pf%results%niters+1, pf%results%max_nsweeps),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)               
    if(.not.allocated(pf%results%residuals)) allocate(pf%results%residuals(pf%results%nlevs, pf%results%nblocks,pf%results%niters+1, pf%results%max_nsweeps),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                   
    if(.not.allocated(pf%results%delta_q0)) allocate(pf%results%delta_q0(pf%results%nlevs,pf%results%nblocks,pf%results%niters+1, pf%results%max_nsweeps),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                   
    if(.not.allocated(pf%results%iters)) allocate(pf%results%iters(pf%results%nblocks),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                   

    pf%results%errors = -1.0_pfdp
    pf%results%residuals = -1.0_pfdp
    pf%results%delta_q0 = -1.0_pfdp
    pf%results%iters = pf%results%niters

    !  Set up the directory to dump results
    istat= system('mkdir -p dat')  !  May be there already
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make dat directory dat in initialize_results")       
    istat= system('mkdir -p dat/' // trim(pf%outdir))       
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make base directory in initialize_results")

    pf%results%datpath= 'dat/' // trim(pf%outdir) // '/'

    ! Create directory for this processor
    write (dirname, "(A5,I0.3)") 'Proc_',pf%results%rank
    datpath=trim(pf%results%datpath) // trim(dirname) 
    istat= system('mkdir -p ' // trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make Proc directory in dump_resids")
    !  Pf%Results is the final path for all the stat files
    pf%results%datpath= trim(datpath)

  end subroutine initialize_results

  subroutine dump_results(this)
    type(pf_results_t), intent(inout) :: this
    integer :: kblock, kiter, ksweep,nstep,klevel
    integer :: istat
    integer :: system
    integer :: estream,rstream,qstream,istream
    character(len = 256) :: iname  !!  output file name for residuals
    character(len = 256) :: resname  !!  output file name for residuals
    character(len = 256) :: errname  !!  output file name for errors
    character(len = 256) :: q0name   !!  output file name for delta_q0
    character(len = 128) :: dirname  !!  directory name


    !  output residuals for each sweep and block and iteration
    resname = trim(this%datpath) // '/residual.dat'
    errname = trim(this%datpath) // '/error.dat'
    q0name = trim(this%datpath) // '/delta_q0.dat'
    rstream=5000+this%rank
    estream=6000+this%rank
    qstream=7000+this%rank

    open(rstream, file=trim(resname), form='formatted')
    open(estream, file=trim(errname), form='formatted')
    open(qstream, file=trim(q0name), form='formatted')
    do klevel=1,this%nlevs
       do kblock = 1, this%nblocks
          nstep=(kblock-1)*this%nprocs+this%rank+1
          do kiter = 0 , this%niters
             do ksweep = 1, this%max_nsweeps
                write(rstream,101 ) klevel,nstep,kblock,kiter,ksweep,this%residuals(klevel, kblock,kiter+1, ksweep)
                write(estream,101) klevel,nstep,kblock,kiter,ksweep,this%errors(klevel,kblock,kiter+1,  ksweep)
                write(qstream,101) klevel,nstep,kblock,kiter,ksweep,this%delta_q0(klevel, kblock,kiter+1, ksweep)
             end do
          end do
       enddo
    enddo
    101 format(I3,I10, I10,I5, I4, e22.14)
    close(rstream)
    close(estream)
    close(qstream)

    !  output residuals only for last iteration
    iname = trim(this%datpath) // '/iter.dat'
    istream=8000+this%rank
    open(istream, file=trim(iname), form='formatted')

    do kblock = 1, this%nblocks
       nstep=(kblock-1)*this%nprocs+this%rank+1
       write(istream, '(I10,I10, e22.14)') nstep,kblock,this%iters(kblock)
    enddo
    close(istream)
    
  end subroutine dump_results

  subroutine dump_timingsl(this,pf)
    type(pf_results_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    character(len = 128   ) :: pname     !!  processor name
    character(len = 256   ) :: fullname  !!  output file name for runtimes
    character(len = 128   ) :: datpath  !!  directory path
    character(len = 128   ) :: strng      !  used for string conversion
    integer :: istat,j, iout,system,nlev,k,kmax
    real(pfdp) :: qarr(pf%nlevels)

    !  Write a json file with timer numbers and times
    fullname = trim(this%datpath) // '/runtime.json'
    iout = 4000+pf%rank !  Use processor dependent file number
    nlev=pf%nlevels
    !  output timings
    open(iout, file=trim(fullname), form='formatted')
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
    
  end subroutine dump_timingsl

  
  subroutine destroy_results(this)
    type(pf_results_t), intent(inout) :: this
    
    if(allocated(this%errors))  deallocate(this%errors)
    if(allocated(this%residuals))  deallocate(this%residuals)
    if(allocated(this%delta_q0))  deallocate(this%delta_q0)
    if(allocated(this%iters))  deallocate(this%iters)
  end subroutine destroy_results

end module pf_mod_results
