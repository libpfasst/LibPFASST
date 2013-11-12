!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module pmain
  use pfasst
  use feval
  implicit none
contains

  subroutine pth_main(ptr) bind(c)
    type(c_ptr), intent(in   ), value :: ptr

    type(pf_pfasst_t), pointer :: pf
    type(ndarray),     pointer :: q0

    call c_f_pointer(ptr, pf)
    call pf_pfasst_setup(pf)

    allocate(q0)
    call ndarray_create_simple(q0, [ pf%levels(pf%nlevels)%nvars ])
    call initial(q0)
    call pf_pfasst_run(pf, c_loc(q0), 0.01d0, 0.0d0, pf%comm%nproc)
    call ndarray_destroy(c_loc(q0))
    call pthread_exit(c_null_ptr)
  end subroutine pth_main

end module pmain
