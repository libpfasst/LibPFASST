!! File holding the module to exit gracefully and informatively
!
! This file is part of LIBPFASST.
!
!> Module to exit gracefully and informatively
module pf_mod_stop
contains  
  subroutine pf_stop(pf_file,Nline,msg, val, rank)
    character(len=*), intent(in) :: pf_file
    integer, intent(in):: Nline
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: val
    integer, intent(in), optional :: rank    

    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    if (present(rank))   print *,'Rank = ',rank
    print *,'Stopping in File: ', pf_file    
    print *,'Line number: ', Nline
    print *,msg
    if (present(val))   print *,' = ',val
    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    call EXIT(-1)
  end subroutine pf_stop
end module pf_mod_stop
