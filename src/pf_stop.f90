!
! This file is part of LIBPFASST.
!
!> Module to exit gracefully and informatively
module pf_mod_stop
contains  
  subroutine pf_stop(pf_file,Nline,msg, N)
    character(len=*), intent(in) :: pf_file
    integer, intent(in):: Nline
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: N

    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *,'Stopping in File: ', pf_file    
    print *,'Line number: ', Nline
    print *,msg
    if (present(N))   print *,'value=',N
    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    call EXIT(-1)
  end subroutine pf_stop
end module pf_mod_stop
