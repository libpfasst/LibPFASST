!!  Restriction operators
!
! This file is part of LIBPFASST.
!
module pf_mod_restrict
  !!  Module to restrict solutions between pfasst levels and create the FAS tau correction
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_hooks
  use pf_mod_utils
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
    class(pf_level_t), pointer :: c_lev_p    
    class(pf_level_t), pointer :: f_lev_p

    integer    :: m, step

    real(pfdp), allocatable :: c_times(:)  !!  Simulation time at coarse nodes  
    real(pfdp), allocatable :: f_times(:)  !!  Simulation time at fine nodes
    
    
    f_lev_p => pf%levels(level_index);
    c_lev_p => pf%levels(level_index-1)
    
    step = pf%state%step+1
    if(present(mystep)) step = mystep
    
    call call_hooks(pf, level_index, PF_PRE_RESTRICT_ALL)
    call start_timer(pf, TRESTRICT + level_index - 1)
    

    allocate(c_times(c_lev_p%nnodes))
    allocate(f_times(f_lev_p%nnodes))

    !> restrict q's and recompute f's
    c_times = t0 + dt*c_lev_p%nodes
    f_times = t0 + dt*f_lev_p%nodes

    call restrict_ts(f_lev_p, c_lev_p, f_lev_p%Q, c_lev_p%Q, f_times, flags)

    !>  Recompute the functions
    call c_lev_p%ulevel%sweeper%evaluate_all(pf,level_index-1, c_times, flags=flags, step=step)

    !>  Compute  FAS correction
    do m = 1, c_lev_p%nnodes-1
       call c_lev_p%tauQ(m)%setval(0.0_pfdp, flags)
    end do
    if (pf%state%iter >= pf%taui0)  then
       ! compute '0 to node' integral on the coarse level
      call c_lev_p%ulevel%sweeper%integrate(pf,level_index-1,  c_lev_p%Q, &
           c_lev_p%F, dt, c_lev_p%I, flags)
       ! compute '0 to node' integral on the fine level
      call f_lev_p%ulevel%sweeper%integrate(pf,level_index, f_lev_p%Q, &
        f_lev_p%F, dt, f_lev_p%I, flags)
       !  put tau in on fine level
      if (level_index < pf%state%finest_level) then
          do m = 1, f_lev_p%nnodes-1
             call f_lev_p%I(m)%axpy(1.0_pfdp, f_lev_p%tauQ(m), flags)
          end do
       end if
       !  Subtract coarse integral
       do m = 1, c_lev_p%nnodes-1
          call c_lev_p%tauQ(m)%axpy(-1.0_pfdp, c_lev_p%I(m), flags)
       end do

       ! restrict '0 to node' integral on the fine level  in time and space
       call restrict_ts_integral(f_lev_p, c_lev_p, f_lev_p%I, c_lev_p%I, f_times, flags)

       ! Add fine restriction of fine integral (stored on coarse)
       do m = 1, c_lev_p%nnodes-1
          call c_lev_p%tauQ(m)%axpy(1.0_pfdp, c_lev_p%I(m), flags)
       end do

!!$       if (pf%use_Sform) then
!!$          do m = c_lev_p%nnodes-1,2,-1
!!$!             call c_lev_p%tauQ(m)%axpy(-1.0_pfdp, c_lev_p%tauQ(m-1), flags)
!!$          end do
!!$       end if
       
    end if

    call end_timer(pf, TRESTRICT + level_index - 1)
    call call_hooks(pf, level_index, PF_POST_RESTRICT_ALL)

    deallocate(c_times)
    deallocate(f_times)
  end subroutine restrict_time_space_fas


  subroutine restrict_ts(f_lev_p, c_lev_p, f_encap_array, c_encap_array, f_time, flags)

    !! Restrict (in time and space) f_sol_array  to c_sol_array
    !! This version is for point values (either functions or solutions)
    
    class(pf_level_t),  intent(inout) :: f_lev_p   !!   pointer to fine level
    class(pf_level_t),  intent(inout) :: c_lev_p   !!   pointer to coarse level
    class(pf_encap_t),  intent(inout) :: f_encap_array(:)   !! array of fine level data to be restricted
    class(pf_encap_t),  intent(inout) :: c_encap_array(:)   !! array of coarse level data to be computed
    real(pfdp),         intent(in) :: f_time(:)             !! time at the fine nodes
    integer, optional, intent(in)    :: flags    

    class(pf_encap_t), allocatable :: f_encap_array_c(:)  !!  fine solution restricted in space only
    integer :: m,j
    integer :: f_nnodes,c_nnodes


    f_nnodes = f_lev_p%nnodes
    c_nnodes = c_lev_p%nnodes

    !!  Create a temp array for the spatial restriction
    call c_lev_p%ulevel%factory%create_array(f_encap_array_c, f_nnodes, c_lev_p%index, c_lev_p%shape)

    !  spatial restriction
    do m = 1, f_nnodes
       call f_lev_p%ulevel%restrict(f_lev_p, c_lev_p, f_encap_array(m), f_encap_array_c(m), f_time(m), flags)
    end do

    ! temporal restriction
    if (present(flags)) then
       if ((flags .eq. 0) .or. (flags .eq. 1)) &
            call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_p%rmat, f_encap_array_c, .true., flags)
       if ((flags .eq. 0) .or. (flags .eq. 2)) &
            call pf_apply_mat_backward(c_encap_array, 1.0_pfdp, f_lev_p%rmat, f_encap_array_c, .true., flags=2)
    else
       call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_p%rmat, f_encap_array_c, .true.)
    end if
    call c_lev_p%ulevel%factory%destroy_array(f_encap_array_c)


  end subroutine restrict_ts

  subroutine restrict_ts_integral(f_lev_p, c_lev_p, f_encap_array, c_encap_array,f_time, flags)

    !! Restrict (in time and space) f_sol_array  to c_sol_array
    !! This version is for integrals
    
    class(pf_level_t),  intent(inout) :: f_lev_p   !!   pointer to fine level
    class(pf_level_t),  intent(inout) :: c_lev_p   !!   pointer to coarse level
    class(pf_encap_t),  intent(inout) :: f_encap_array(:)   !! array of fine level data to be restricted
    class(pf_encap_t),  intent(inout) :: c_encap_array(:)   !! array of coarse level data to be computed
    real(pfdp),         intent(in) :: f_time(:)             !! time at the fine nodes
    integer, optional, intent(in)    :: flags    

    class(pf_encap_t), allocatable :: f_encap_array_c(:)  !!  fine solution restricted in space only
    integer :: m,j
    integer :: f_nnodes,c_nnodes


    f_nnodes = f_lev_p%nnodes
    c_nnodes = c_lev_p%nnodes

    !!  Create a temp array for the spatial restriction
    call c_lev_p%ulevel%factory%create_array(f_encap_array_c, f_nnodes-1, c_lev_p%index, c_lev_p%shape)

    !  spatial restriction
    do m = 1, f_nnodes-1
       call f_lev_p%ulevel%restrict(f_lev_p, c_lev_p, f_encap_array(m), f_encap_array_c(m), f_time(m), flags)
    end do
    
    ! temporal restriction
    ! when restricting '0 to node' integral terms, skip the first entry since it is zero
    if (present(flags)) then
       if ((flags .eq. 0) .or. (flags .eq. 1)) &
            call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_p%rmat(2:,2:), f_encap_array_c, .true., flags=1)
       if ((flags .eq. 0) .or. (flags .eq. 2)) &
            call pf_apply_mat_backward(c_encap_array, 1.0_pfdp, f_lev_p%rmat(2:,2:), f_encap_array_c, .true., flags=2)
    else
       call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev_p%rmat(2:,2:), f_encap_array_c, .true.)
    end if
    call c_lev_p%ulevel%factory%destroy_array(f_encap_array_c)

  end subroutine restrict_ts_integral


end module pf_mod_restrict
