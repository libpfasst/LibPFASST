!!  Restriction operators
!
! This file is part of LIBPFASST.
!
module pf_mod_restrict
  !!  Module to restrict solutions between pfasst levels and create the FAS tau correction
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_hooks
  implicit none
contains



  subroutine restrict_time_space_fas(pf, t0, dt, level_index, flags, mystep)
    !! Restrict (in time and space) fine level to coarse and set coarse level FAS correction.
    !!
    !! The coarse function values are re-evaluated after restriction.
    !! Note that even if the number of variables and nodes is the same,
    !! we should still compute the FAS correction since the function
    !! evaluations may be different.
    type(pf_pfasst_t), intent(inout),target :: pf
    real(pfdp),        intent(in)    :: t0            !!  time at beginning of step
    real(pfdp),        intent(in)    :: dt            !!  time step
    integer,           intent(in)    :: level_index   !! defines which level to restrict
    integer, optional, intent(in)    :: flags, mystep    

    !>  Local variables
    class(pf_level_t), pointer :: c_lev_ptr    
    class(pf_level_t), pointer :: f_lev_ptr

    integer    :: m, step

    real(pfdp), allocatable :: c_times(:)  !!  Simulation time at coarse nodes  
    real(pfdp), allocatable :: f_times(:)  !!  Simulation time at fine nodes
    class(pf_encap_t), allocatable :: &
         c_tmp_array(:), &    ! coarse integral of coarse function values
         f_int_array(:), &    ! fine integral of fine function values
         f_int_arrayr(:)      ! coarse integral of restricted fine function values
    
    f_lev_ptr => pf%levels(level_index);
    c_lev_ptr => pf%levels(level_index-1)

    step = pf%state%step+1
    if(present(mystep)) step = mystep
    
    call call_hooks(pf, level_index, PF_PRE_RESTRICT_ALL)
    call start_timer(pf, TRESTRICT + level_index - 1)
    
    !> create workspaces
    call c_lev_ptr%ulevel%factory%create_array(c_tmp_array, c_lev_ptr%nnodes, &
      c_lev_ptr%index,   c_lev_ptr%shape)
    call c_lev_ptr%ulevel%factory%create_array(f_int_arrayr, c_lev_ptr%nnodes, &
      c_lev_ptr%index,   c_lev_ptr%shape)
    call c_lev_ptr%ulevel%factory%create_array(f_int_array, f_lev_ptr%nnodes, &
      f_lev_ptr%index,   f_lev_ptr%shape)
    allocate(c_times(c_lev_ptr%nnodes))
    allocate(f_times(f_lev_ptr%nnodes))

    !> restrict q's and recompute f's
    c_times = t0 + dt*c_lev_ptr%nodes
    f_times = t0 + dt*f_lev_ptr%nodes

    call restrict_sdc(f_lev_ptr, c_lev_ptr, f_lev_ptr%Q, c_lev_ptr%Q, .false., f_times, flags)

    !>  Recompute the functions
    call c_lev_ptr%ulevel%sweeper%evaluate_all(c_lev_ptr, c_times, flags=flags, step=step)


    !>  Compute  FAS correction
    do m = 1, c_lev_ptr%nnodes-1
       call c_lev_ptr%tauQ(m)%setval(0.0_pfdp, flags)
    end do
    if (pf%state%iter >= pf%taui0)  then
       ! compute '0 to node' integral on the coarse level
      call c_lev_ptr%ulevel%sweeper%integrate(c_lev_ptr, c_lev_ptr%Q, &
        c_lev_ptr%F, dt, c_tmp_array, flags)
       ! compute '0 to node' integral on the fine level
      call f_lev_ptr%ulevel%sweeper%integrate(f_lev_ptr, f_lev_ptr%Q, &
        f_lev_ptr%F, dt, f_lev_ptr%I, flags)
       !  put tau in on fine level
      if (level_index < pf%state%finest_level) then
          do m = 1, f_lev_ptr%nnodes-1
             call f_lev_ptr%I(m)%axpy(1.0_pfdp, f_lev_ptr%tauQ(m), flags)
          end do
       end if
  
       ! restrict '0 to node' integral on the fine level  in time and space
       call restrict_sdc(f_lev_ptr, c_lev_ptr, f_lev_ptr%I, f_int_arrayr, .true.,f_times, flags)

      ! compute '0 to node' tau correction
       do m = 1, c_lev_ptr%nnodes-1
          call c_lev_ptr%tauQ(m)%axpy(1.0_pfdp, f_int_arrayr(m), flags)
          call c_lev_ptr%tauQ(m)%axpy(-1.0_pfdp, c_tmp_array(m), flags)
       end do
    end if

    call end_timer(pf, TRESTRICT + level_index - 1)
    call call_hooks(pf, level_index, PF_POST_RESTRICT_ALL)

    !>  Clean up
    call c_lev_ptr%ulevel%factory%destroy_array(c_tmp_array, c_lev_ptr%nnodes, &
      c_lev_ptr%index,   c_lev_ptr%shape)
    call c_lev_ptr%ulevel%factory%destroy_array(f_int_arrayr, c_lev_ptr%nnodes, &
      c_lev_ptr%index,  c_lev_ptr%shape)
    call f_lev_ptr%ulevel%factory%destroy_array(f_int_array, f_lev_ptr%nnodes, &
      f_lev_ptr%index,   f_lev_ptr%shape)

    deallocate(c_times)
    deallocate(f_times)
  end subroutine restrict_time_space_fas


  subroutine restrict_sdc(f_lev_ptr, c_lev_ptr, f_encap_array, c_encap_array, IS_INTEGRAL,f_time, flags)

    !! Restrict (in time and space) f_sol_array  to c_sol_array
    !! Depending on the flag INTEGRAL, we may be restricting solutions, or integrals of F
    
    class(pf_level_t),  intent(inout) :: f_lev_ptr   !!   pointer to fine level
    class(pf_level_t),  intent(inout) :: c_lev_ptr   !!   pointer to coarse level
    class(pf_encap_t),  intent(inout) :: f_encap_array(:)   !! array of fine level data to be restricted
    class(pf_encap_t),  intent(inout) :: c_encap_array(:)   !! array of coarse level data to be computed
    logical,            intent(in)    :: IS_INTEGRAL       !! flag determines if it is integral data being restricted
    real(pfdp),         intent(in) :: f_time(:)             !! time at the fine nodes
    integer, optional, intent(in)    :: flags    

    class(pf_encap_t), allocatable :: f_encap_array_c(:)  !!  fine solution restricted in space only
    integer :: m,j
    integer :: f_nnodes,c_nnodes


    f_nnodes = f_lev_ptr%nnodes
    c_nnodes = c_lev_ptr%nnodes

    !!  do the restriction
    if (IS_INTEGRAL) then   ! Restriction of integrals
       call c_lev_ptr%ulevel%factory%create_array(f_encap_array_c, f_nnodes-1, c_lev_ptr%index, c_lev_ptr%shape)
       !  spatial restriction
       do m = 1, f_nnodes-1
          call f_lev_ptr%ulevel%restrict(f_lev_ptr, c_lev_ptr, f_encap_array(m), f_encap_array_c(m), f_time(m), flags)
       end do

       ! temporal restriction
       ! when restricting '0 to node' integral terms, skip the first entry since it is zero
       if (present(flags)) then
          if ((flags .eq. 0) .or. (flags .eq. 1)) &
            call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat(2:,2:), f_encap_array_c, .true., flags=1)
          if ((flags .eq. 0) .or. (flags .eq. 2)) &
            call pf_apply_mat_backward(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat(2:,2:), f_encap_array_c, .true., flags=2)
       else
          call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat(2:,2:), f_encap_array_c, .true.)
       end if
       call c_lev_ptr%ulevel%factory%destroy_array(f_encap_array_c, f_nnodes-1, c_lev_ptr%index, c_lev_ptr%shape)
    else
       call c_lev_ptr%ulevel%factory%create_array(f_encap_array_c, f_nnodes, c_lev_ptr%index, c_lev_ptr%shape)
       !  spatial restriction
       do m = 1, f_nnodes
          call f_lev_ptr%ulevel%restrict(f_lev_ptr, c_lev_ptr, f_encap_array(m), f_encap_array_c(m), f_time(m), flags)
       end do
       ! temporal restriction
       if (present(flags)) then
          if ((flags .eq. 0) .or. (flags .eq. 1)) &
            call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat, f_encap_array_c, .true., flags)
          if ((flags .eq. 0) .or. (flags .eq. 2)) &
            call pf_apply_mat_backward(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat, f_encap_array_c, .true., flags=2)
        else
           call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_ptr%rmat, f_encap_array_c, .true.)
        end if
       call c_lev_ptr%ulevel%factory%destroy_array(f_encap_array_c, f_nnodes, c_lev_ptr%index, c_lev_ptr%shape)
    end if

  end subroutine restrict_sdc

  subroutine pf_apply_mat(dst, a, mat, src, zero, flags)
    !! Apply a matrix (tmat or rmat) to src and add to dst.
    !! Mathematically this is 
    !!     dst= dst + a*mat*src
    !!  Where dst and src are vectors, mat is a matrix, and a is a scalar
    !!  If the optional variable "zero" is provided and is true, then we compute
    !!     dst=  a*mat*src
    class(pf_encap_t), intent(inout) :: dst(:)       !!  destination vector
    real(pfdp),        intent(in)    :: a            !!  scalar
    real(pfdp),        intent(in)    :: mat(:, :)    !!  matrix
    class(pf_encap_t), intent(in)    :: src(:)       !!  src vector
    logical,           intent(in), optional :: zero   !! If false, don't zero out the the dst variable before computing 
    integer,           intent(in), optional :: flags  !! Used for choosing which variable to operate on 
    
    !!  Local variables
    logical :: lzero   !!  local version of input parameter zero
    integer :: which   !!  local version of flags
    integer :: n, m    !!  size of mat   
    integer :: i, j    !!  loop variables

    lzero = .true.; if (present(zero)) lzero = zero    
    which = 1;      if(present(flags)) which = flags
        
    n = size(mat, dim=1)
    m = size(mat, dim=2)
        
    do i = 1, n
      if (lzero) call dst(i)%setval(0.0_pfdp, flags)
      do j = 1, m
         if (a*mat(i, j) /= 0.0_pfdp)  call dst(i)%axpy(a * mat(i, j), src(j), flags)
      end do
    end do
  end subroutine pf_apply_mat
  

  subroutine pf_apply_mat_backward(dst, a, mat, src, zero, flags)
    !! Apply a matrix (tmat or rmat) to src and add to dst.
    class(pf_encap_t), intent(inout) :: dst(:)       !!  destination vector
    real(pfdp),        intent(in)    :: a            !!  scalar
    real(pfdp),        intent(in)    :: mat(:, :)    !!  matrix
    class(pf_encap_t), intent(in)    :: src(:)       !!  src vector
    logical,           intent(in), optional :: zero   !! If false, don't zero out the the dst variable before computing 
    integer,           intent(in), optional :: flags  !! Used for choosing which variable to operate on 

    
    !!  Local variables
    logical :: lzero   !!  local version of input parameter zero
    integer :: which   !!  local version of flags
    integer :: n, m    !!  size of mat   
    integer :: i, j    !!  loop variables

    lzero = .true.; if (present(zero)) lzero = zero    
    which = 2;      if(present(flags)) which = flags
    
    if( which /= 2 ) &
      stop "pf_apply_mat_backward can only be used for restricting the backward integrals with which==2"

    n = size(mat, dim=1)
    m = size(mat, dim=2)
        
    do i = 1, n
      if (lzero) call dst(n+1-i)%setval(0.0_pfdp, 2)
      do j = 1, m
        if (a*mat(i, j) /= 0.0_pfdp)  call dst(n+1-i)%axpy(a * mat(i, j), src(m+1-j), 2)
      end do
    end do
  end subroutine pf_apply_mat_backward

end module pf_mod_restrict
