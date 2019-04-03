!! Some convenient use statements
!
! This file is part of LIBPFASST.
!
!>  Module containing a collection of "use" statements to simplify
!!  including the common main modules in writing applications that use libpfasst
module pfasst
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_results
  use pf_mod_parallel
  use pf_mod_pfasst
  use pf_mod_utils
#ifndef NOMPI
  use pf_mod_comm_mpi
#endif
  
end module pfasst

