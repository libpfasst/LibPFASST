!
! This file is part of LIBPFASST.
!
!>  Module for the calling of user defined routines from various places in the pfasst algorithm
module pf_mod_results
  use pf_mod_dtype
  implicit none

  

  
contains
    subroutine initialize_results(this, nsteps_in, niters_in, nprocs_in, nlevels_in,rank_in)
    class(pf_results_t), intent(inout) :: this
    integer, intent(in) :: nsteps_in, niters_in, nprocs_in, nlevels_in,rank_in

    if (rank_in == 0) then
       open(unit=123, file='result-size.dat', form='formatted')
       write(123,'(I5, I5, I5, I5)') nsteps_in, niters_in, nprocs_in, nlevels_in
       close(unit=123)
    end if


    this%dump => dump_results
    this%destroy => destroy_results

    this%nsteps=nsteps_in
    this%nblocks=nsteps_in/nprocs_in
    this%niters=niters_in
    this%nprocs=nprocs_in
    this%nlevels=nlevels_in
    this%p_index=rank_in+100

    write (this%fname_r, "(A13,I0.3,A4)") 'dat/residual_',rank_in,'.dat'
    write (this%fname_e, "(A10,I0.3,A4)") 'dat/errors_',rank_in,'.dat'

    allocate(this%errors(niters_in, this%nblocks, nlevels_in), &
         this%residuals(niters_in, this%nblocks, nlevels_in))

    this%errors = 0.0_pfdp
    this%residuals = 0.0_pfdp
  end subroutine initialize_results

  subroutine dump_results(this)
    type(pf_results_t), intent(inout) :: this
    integer :: i, j, k
    
    open(unit=this%p_index, file=this%fname_r, form='formatted')
    do k = 1, this%nlevels
       do j = 1, this%nblocks
          do i = 1 , this%niters
             write(this%p_index, '(I4, I4, I4, e21.14)') i,j,k,this%residuals(i, j, k)
          end do
       end do
    enddo
    close(this%p_index)

  end subroutine dump_results

  subroutine destroy_results(this)
    type(pf_results_t), intent(inout) :: this

    if(allocated(this%errors))  deallocate(this%errors)
    if(allocated(this%residuals))  deallocate(this%residuals)
  end subroutine destroy_results

end module pf_mod_results
