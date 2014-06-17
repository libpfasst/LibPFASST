module hooks
  use pfasst
  use encap
  use initial
  implicit none
contains

  subroutine project_hook(pf, level, state, ctx)
    use feval, only: project
    implicit none

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    integer :: n
    complex(c_double_complex), pointer :: ustar(:,:,:,:)
    type(carray4), pointer :: qend

    n = level%shape(1)
    allocate(ustar(n,n,n,3))
    call c_f_pointer(level%qend, qend)

    !$omp parallel workshare
    ustar = qend%array
    !$omp end parallel workshare

    call project(ctx, n, n, n, ustar, qend%array)

    deallocate(ustar)
  end subroutine project_hook

  subroutine echo_error_hook(pf, level, state, levelctx)
    use feval
    implicit none

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: levelctx

    type(carray4)          :: qex
    type(carray4), pointer :: q
    type(feval_t), pointer :: ctx

    real(pfdp) :: e0, e1, r

    call c_f_pointer(levelctx, ctx)
    call carray4_create(qex, level%shape)

    call c_f_pointer(level%q(1), q)
    call exact(qex, ctx%nu, state%t0)
    e0 = maxval(abs(qex%array-q%array))

    call c_f_pointer(level%qend, q)
    call exact(qex, ctx%nu, state%t0+state%dt)
    e1 = maxval(abs(qex%array-q%array))

    call c_f_pointer(level%R(level%nnodes-1), q)
    r = maxval(abs(q%array))

    print *, 'err0, err1, res: ', e0, e1, r

    deallocate(qex%array)
  end subroutine echo_error_hook

  subroutine dump(dname, fname, q)
    use pf_mod_ndarray
    character(len=*), intent(in   ) :: dname, fname
    type(carray4),    intent(in   ) :: q

    real(c_double), pointer :: rarray(:)
    type(c_ptr) :: tmp

    tmp = c_loc(q%array(1,1,1,1))
    call c_f_pointer(tmp, rarray, [ 2*product(q%shape) ])

    call ndarray_dump_numpy(trim(dname)//c_null_char, trim(fname)//c_null_char, '<c16'//c_null_char, &
         4, q%shape, 2*product(q%shape), rarray)

  end subroutine dump

  subroutine dump_hook(pf, level, state, levelctx)
    use pf_mod_ndarray

    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: levelctx

    character(len=256)     :: fname
    type(carray4), pointer :: qend

    if (state%iter < 0) return

    call c_f_pointer(level%qend, qend)

    write(fname, "('s',i0.5,'i',i0.3,'l',i0.2,'.npy')") &
         state%step, state%iter, level%level

    call dump(pf%outdir, fname, qend)

  end subroutine dump_hook

end module hooks
