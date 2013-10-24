

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

end module hooks
