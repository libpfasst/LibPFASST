!
! This file is part of LIBPFASST.
!
!>  Module containing a collection of "use" statements to simplify
!!  including the common main modules in writing applications that use libpfasst
module pfasst
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_parallel
  use pf_mod_pfasst
#ifndef NOMPI
  use pf_mod_comm_mpi
#endif
  use pf_mod_imexQ
end module pfasst

