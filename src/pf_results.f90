!! Storing results for eventual output
!
! This file is part of LIBPFASST.
!
!>  Module for the storing results for eventual output
module pf_mod_results
  use pf_mod_dtype
  implicit none

  

  
contains
    subroutine initialize_results(this, nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,lev_ind)
    class(pf_results_t), intent(inout) :: this
    integer, intent(in) :: nsteps_in, niters_in, nprocs_in, nsweeps_in,rank_in,lev_ind
    character(len = 24   ) :: fname  !!  output file name for residuals
    
    if (rank_in == 0) then
       write (fname, "(A17,I0.3,A4)") 'dat/results_size_',lev_ind,'.dat'
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

    print *,'initialize results',niters_in, this%nblocks, nsweeps_in,lev_ind

    if(.not.allocated(this%errors)) allocate(this%errors(niters_in, this%nblocks, nsweeps_in))
    if(.not.allocated(this%residuals)) allocate(this%residuals(niters_in, this%nblocks, nsweeps_in))

    this%errors = 0.0_pfdp
    this%residuals = 0.0_pfdp
  end subroutine initialize_results

  subroutine dump_results(this)
    type(pf_results_t), intent(inout) :: this
    integer :: i, j, k
    character(len = 24   ) :: fname_r  !!  output file name for residuals
    character(len = 21) :: fname_e     !!  output file name errors
    

    write (fname_r, "(A13,I0.3,A1,I0.3,A4)") 'dat/residual_',this%rank,'_',this%level,'.dat'
    write (fname_e, "(A10,I0.3,A1,I0.3,A4)") 'dat/errors_',  this%rank,'_',this%level,'.dat'

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

  subroutine destroy_results(this)
    type(pf_results_t), intent(inout) :: this

    if(allocated(this%errors))  deallocate(this%errors)
    if(allocated(this%residuals))  deallocate(this%residuals)
  end subroutine destroy_results

end module pf_mod_results
