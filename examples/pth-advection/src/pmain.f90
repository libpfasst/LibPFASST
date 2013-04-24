!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module pmain
  use pfasst
  implicit none
contains

  subroutine pth_main(ptr) bind(c)
    use encap_array1d
    use feval

    type(c_ptr), intent(in), value :: ptr
    type(pf_pfasst_t), pointer     :: pf

    type(array1d), target :: q0
    integer :: l

    call c_f_pointer(ptr, pf)
    call pf_pfasst_setup(pf)

    allocate(q0%array(pf%levels(pf%nlevels)%nvars))
    call initial(q0)

    call pf_pfasst_run(pf, c_loc(q0), 0.01d0, 0.0d0, pf%comm%nproc)

    deallocate(q0%array)

    call pthread_exit(c_null_ptr)
  end subroutine pth_main

end module pmain
