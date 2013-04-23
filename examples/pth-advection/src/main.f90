!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

program fpfasst
  use pfasst
  use pmain
  use feval
  use hooks
  use transfer
  use encap_array1d
  implicit none

  type(pf_pfasst_t),  pointer :: pf
  type(pf_comm_t),    target  :: tcomm
  type(pf_sweeper_t), target  :: sweeper
  type(pf_encap_t),   target  :: encap


  integer(c_long), pointer :: tid(:)
  integer(c_long) :: ret

  integer :: l, t, nthreads
  integer, target :: nvars(3), nnodes(3)

  nvars  = [ 16, 32, 64 ]
  nnodes = [ 3, 5, 9 ]

  nthreads = 4


  !
  ! initialize pfasst and launch
  !
  call pf_pthreads_create(tcomm, nthreads, size(nvars))

  call array1d_encap_create(encap)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)

  allocate(tid(nthreads))
  do t = 1, nthreads
     allocate(pf)
     call pf_pfasst_create(pf, tcomm, size(nvars))

     pf%rank   = t - 1
     pf%niters = 12
     pf%qtype  = SDC_GAUSS_LOBATTO
     pf%echo_timings = .false.
     pf%levels(1)%nsweeps = 2

     do l = 1, size(nvars)
        pf%levels(l)%nvars  = nvars(l)
        pf%levels(l)%nnodes = nnodes(l)

        call feval_create_workspace(pf%levels(l)%ctx, nvars(l))

        pf%levels(l)%interpolate => interpolate
        pf%levels(l)%restrict    => restrict
        pf%levels(l)%encap       => encap
        pf%levels(l)%sweeper     => sweeper
     end do

     call pf_pthreads_setup(tcomm, pf)
     call pf_add_hook(pf, 3, PF_POST_ITERATION, echo_error)

     ret = pthread_create(c_loc(tid(t)), c_null_ptr, &
          c_funloc(pth_main), c_loc(pf))
  end do

  do t = 1, nthreads
     ret = pthread_join(tid(t), c_null_ptr)
  end do

  do t = 1, nthreads
     call c_f_pointer(tcomm%pfs(t-1), pf)
     do l = 1, size(nvars)
        call feval_destroy_workspace(pf%levels(l)%ctx)
     end do
  end do

  call pf_pthreads_destroy(tcomm)
  call pf_imex_destroy(sweeper)
  call fftw_cleanup()

  deallocate(tid)

  call pthread_exit(c_null_ptr)
end program fpfasst
