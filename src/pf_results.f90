!! Storing results for eventual output
!
! This file is part of LIBPFASST.
!
!>  Module for the storing results for eventual output
module pf_mod_results
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  
contains
    subroutine initialize_results(this, nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,level_index,datpath)
    class(pf_results_t), intent(inout) :: this
    integer, intent(in) :: nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,level_index
    character(len=*), intent(in) :: datpath

    character(len = 50   ) :: fname  !!  output file name for residuals
    character(len = 100   ) :: fullname  !!  output file name for residuals
    integer :: istat
    
    !  Set up the directory to dump results
    istat= system('mkdir -p ' // trim(datpath))
    this%datpath= trim(datpath) // '/'
    
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in initialize_results")

    write (fname, "(A17,I0.1,A4)") 'residuals_size_L',level_index,'.dat'
    fullname = trim(this%datpath) // trim(fname)
    
    if (rank_in == 0) then
       open(unit=123, file=trim(fullname), form='formatted')
       write(123,'(I5, I5, I5, I5)') nsteps_in, niters_in, nprocs_in, nsweeps_in
       close(unit=123)
    end if
    
    this%dump => dump_results
    this%destroy => destroy_results

    this%nsteps=nsteps_in
    this%nblocks=nsteps_in/nprocs_in
    this%niters=niters_in
    this%nprocs=nprocs_in
    this%nsweeps=nsweeps_in
    this%rank=rank_in
    this%level=level_index    
    
    if(.not.allocated(this%errors)) allocate(this%errors(niters_in, this%nblocks, nsweeps_in))
    if(.not.allocated(this%residuals)) allocate(this%residuals(niters_in, this%nblocks, nsweeps_in))

    this%errors = 0.0_pfdp
    this%residuals = 0.0_pfdp
  end subroutine initialize_results

  subroutine dump_results(this)
    type(pf_results_t), intent(inout) :: this
    integer :: i, j, k, istat
    character(len = 50   ) :: fname  !!  output file name for residuals
    character(len = 100   ) :: fullname  !!  output file name for residuals
    character(len = 50   ) :: datpath  !!  directory path

    
    datpath = trim(this%datpath) // 'residuals'
    istat= system('mkdir -p ' // trim(datpath))
    
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_results")

    write (fname, "(A6,I0.3,A5,I0.1,A4)") '/Proc_',this%rank,'_Lev_',this%level,'.dat'
    fullname = trim(datpath) // trim(fname)
    !  output residuals
    open(100+this%rank, file=trim(fullname), form='formatted')
    do j = 1, this%nblocks
       do i = 1 , this%niters
          do k = 1, this%nsweeps
             write(100+this%rank, '(I4, I4, I4, e21.14)') j,i,k,this%residuals(i, j, k)
          end do
       end do
    enddo
    close(100+this%rank)

  end subroutine dump_results

    subroutine dump_timings(pf)
    type(pf_pfasst_t), intent(inout) :: pf
    character(len = 25   ) :: fname  !!  output file name for runtimes
    character(len = 50   ) :: fullname  !!  output file name for runtimes
    character(len = 50   ) :: datpath  !!  directory path
    integer :: istat,j, istream

    datpath = trim(pf%outdir) // 'runtimes'

    istat= system('mkdir -p '// trim(datpath))
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_timings")

    write (fname, "(A6,I0.3,A4)")  '/Proc_',pf%rank,'.dat'    
    fullname = trim(datpath) // trim(fname)
    istream = 200+pf%rank !  Use processor dependent file number
    !  output timings
    open(istream, file=trim(fullname), form='formatted')
    do j = 1, 100
       if (pf%runtimes(j) > 0.0d0) then
          write(istream, '(a16,  f23.8)') timer_names(j),pf%runtimes(j)
       end if
    end do
       
    close(istream)

  end subroutine dump_timings

  subroutine destroy_results(this)
    type(pf_results_t), intent(inout) :: this

    if(allocated(this%errors))  deallocate(this%errors)
    if(allocated(this%residuals))  deallocate(this%residuals)
  end subroutine destroy_results

end module pf_mod_results
