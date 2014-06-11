module hooks
  use pfasst
  use encap
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

  subroutine dump(dname, fname, q)
    use pf_mod_ndarray
    character(len=*), intent(in   ) :: dname, fname
    type(carray4),    intent(in   ) :: q

    real(c_double), pointer :: rarray(:)
    type(c_ptr) :: tmp

    tmp = c_loc(q%array(1,1,1,1))
    call c_f_pointer(tmp, rarray, [ product(2*q%shape) ])

    call ndarray_dump_numpy(trim(dname)//c_null_char, trim(fname)//c_null_char, '<c16'//c_null_char, &
         4, q%shape, product(2*q%shape), rarray)

  end subroutine dump



  subroutine dump_hook(pf, level, state, levelctx)
    use pf_mod_ndarray

    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: levelctx

    character(len=256)     :: fname
    type(carray4), pointer :: qend

    call c_f_pointer(level%qend, qend)

    write(fname, "('s',i0.5,'i',i0.3,'l',i0.2,'.npy')") &
         state%step, state%iter, level%level

    call dump(pf%outdir, fname, qend)

  end subroutine dump_hook

end module hooks
