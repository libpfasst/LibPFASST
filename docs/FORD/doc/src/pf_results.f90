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
    subroutine initialize_results(this, nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,lev_ind)
    class(pf_results_t), intent(inout) :: this
    integer, intent(in) :: nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,lev_ind
    character(len = 25   ) :: fname  !!  output file name for residuals
    integer :: istat
    
    !  Set up the directory to dump results
    istat= system('mkdir -p dat')
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in initialize_results")

    if (rank_in == 0) then
       write (fname, "(A20,I0.1,A4)") 'dat/residuals_size_L',lev_ind,'.dat'
       open(unit=123, file=fname, form='formatted')
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
    this%level=lev_ind    

    if(.not.allocated(this%errors)) allocate(this%errors(niters_in, this%nblocks, nsweeps_in))
    if(.not.allocated(this%residuals)) allocate(this%residuals(niters_in, this%nblocks, nsweeps_in))

    this%errors = 0.0_pfdp
    this%residuals = 0.0_pfdp
  end subroutine initialize_results

  subroutine dump_results(this)
    type(pf_results_t), intent(inout) :: this
    integer :: i, j, k
    character(len = 32   ) :: fname_r  !!  output file name for residuals
    character(len = 25   ) :: fname_t  !!  output file name for runtimes
    character(len = 21) :: fname_e     !!  output file name errors
    integer :: istat
    
    istat= system('mkdir -p dat/residuals')
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in initialize_results")

    write (fname_r, "(A19,I0.3,A5,I0.1,A4)") 'dat/residuals/Proc_',this%rank,'_Lev_',this%level,'.dat'
!    write (fname_e, "(A10,I0.3,A1,I0.3,A4)") 'dat/errors_',  this%rank,'_',this%level,'.dat'

    !  output residuals
    open(100+this%rank, file=fname_r, form='formatted')
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
    integer :: istat,j, istream
    
    istat= system('mkdir -p dat/runtimes')
    if (istat .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory in dump_timings")

    write (fname, "(A18,I0.3,A4)")         'dat/runtimes/Proc_',pf%rank,'.dat'    

    istream = 200+pf%rank !  Use processor dependent file number
    !  output timings
    open(istream, file=fname, form='formatted')
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
