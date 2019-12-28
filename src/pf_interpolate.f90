!! Interpolation operators
!
! This file is part of LIBPFASST.
!
!> Module to do interpolation between pfasst levels
module pf_mod_interpolate
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_hooks
  use pf_mod_utils
  implicit none
contains

  !> Subroutine to interpolate (in time and space) level_index-1 to level_index
  !! Interpolation is done by interpolating increments.  
  !! The fine function values are re-evaluated after interpolation.
  subroutine interpolate_time_space(pf, t0, dt, level_index, F_INTERP, flags)
    type(pf_pfasst_t), intent(inout),target :: pf      !! main pfasst structure
    real(pfdp),        intent(in)    :: t0             !! time at beginning of time interval
    real(pfdp),        intent(in)    :: dt             !! time step
    integer,           intent(in)    :: level_index    !! defines which level to interpolate to
    logical,           intent(in)    :: F_INTERP !!  Flag, if true, then do interp on f not sol
    integer, optional, intent(in)    :: flags

    !  Local variables
    class(pf_level_t), pointer :: c_lev   !  Pointer to coarse level
    class(pf_level_t), pointer :: f_lev   !  Pointer to fine level

    integer    :: m, p, step,ierr
    real(pfdp), allocatable :: c_times(:)   ! coarse level node times
    real(pfdp), allocatable :: f_times(:)   ! fine level node times


    f_lev => pf%levels(level_index)   ! fine level
    c_lev => pf%levels(level_index-1) ! coarse level

    call call_hooks(pf, level_index, PF_PRE_INTERP_ALL)
    if (pf%save_timings > 1) call pf_start_timer(pf, T_INTERPOLATE,level_index)
    
    
    step = pf%state%step+1

    !> create workspaces
    if (f_lev%interp_workspace_allocated   .eqv. .false.) then  
       call c_lev%ulevel%factory%create_array(f_lev%c_delta,  c_lev%nnodes, c_lev%index,  c_lev%lev_shape)
       call f_lev%ulevel%factory%create_array(f_lev%cf_delta, c_lev%nnodes, f_lev%index,  f_lev%lev_shape)
       f_lev%interp_workspace_allocated  = .true.     
    end if
    !> set time at coarse and fine nodes
    allocate(c_times(c_lev%nnodes),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)
    allocate(f_times(f_lev%nnodes),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    

    c_times = t0 + dt*c_lev%nodes
    f_times = t0 + dt*f_lev%nodes

    do m = 1, c_lev%nnodes
       call f_lev%c_delta(m)%setval(0.0_pfdp,flags)
       call f_lev%cf_delta(m)%setval(0.0_pfdp,flags)
    end do

    !>  interpolate coarse level correction in space only
    do m = 1, c_lev%nnodes
       call f_lev%c_delta(m)%copy(c_lev%Q(m), flags)
       call f_lev%c_delta(m)%axpy(-1.0_pfdp, c_lev%pQ(m), flags)
       call f_lev%ulevel%interpolate(f_lev,c_lev, f_lev%cf_delta(m), f_lev%c_delta(m), c_times(m), flags)
    end do

    !> interpolate corrections in time
    call pf_apply_mat(f_lev%Q, 1.0_pfdp, f_lev%tmat, f_lev%cf_delta, .false., flags)

    !> either interpolate function values or recompute them
    if (F_INTERP) then         !  Interpolating F
      do p = 1,SIZE(c_lev%F(1,:))
          do m = 1, c_lev%nnodes
             call f_lev%c_delta(m)%setval(0.0_pfdp, flags)
             call f_lev%cf_delta(m)%setval(0.0_pfdp, flags)
          end do
          ! interpolate coarse corrections  in space
          do m = 1, c_lev%nnodes
            call f_lev%c_delta(m)%copy(c_lev%F(m,p), flags)
            call f_lev%c_delta(m)%axpy(-1.0_pfdp, c_lev%pF(m,p), flags)
            call f_lev%ulevel%interpolate(f_lev, c_lev, f_lev%cf_delta(m), f_lev%c_delta(m), c_times(m), flags)
         end do

         ! interpolate corrections  in time
          call pf_apply_mat(f_lev%F(:,p), 1.0_pfdp, f_lev%tmat, f_lev%cf_delta, .false., flags)

       end do !  Loop on npieces
    else    ! recompute function values
       call f_lev%ulevel%sweeper%evaluate_all(pf,level_index, f_times, flags=flags, step=step)
    end if  !  Feval

    !> destroy local data structures
    deallocate(c_times,f_times)


    if (pf%save_timings > 1) call pf_stop_timer(pf, T_INTERPOLATE,f_lev%index)
    call call_hooks(pf, f_lev%index, PF_POST_INTERP_ALL)
  end subroutine interpolate_time_space

  !>  Subroutine to update the fine initial condition from coarse increment by spatial interpolation
  subroutine interpolate_q0(pf, f_lev, c_lev, flags)

    type(pf_pfasst_t), intent(inout) :: pf          !!  main pfasst structure
    class(pf_level_t),  intent(inout) :: f_lev  !!  fine level
    class(pf_level_t),  intent(inout) :: c_lev  !!  coarse level
    integer, optional, intent(in)    :: flags       !!  optional: specify component on which to operate
                                                    !   here flags more or less is logical, if it is present we operate on component 1
                                                    !   of the ndarray-type


    call call_hooks(pf, f_lev%index, PF_PRE_INTERP_Q0)
    if (pf%save_timings > 1) call pf_start_timer(pf, T_INTERPOLATE, f_lev%index )


    call c_lev%delta_q0%setval(0.0_pfdp,flags)
    call f_lev%delta_q0%setval(0.0_pfdp,flags)
    

    !>  restrict fine initial data to coarse
    call f_lev%ulevel%restrict(f_lev, c_lev, f_lev%q0, c_lev%delta_q0, pf%state%t0, flags)
    !>  get coarse level correction
    call c_lev%delta_q0%axpy(-1.0_pfdp, c_lev%q0, flags)    
    !>  interpolate correction in space
    call f_lev%ulevel%interpolate(f_lev, c_lev, f_lev%delta_q0, c_lev%delta_q0, pf%state%t0, flags)
    !> update fine inital condition

    call f_lev%q0%axpy(-1.0_pfdp, f_lev%delta_q0, flags)
    if (pf%save_timings > 1) call pf_stop_timer(pf, T_INTERPOLATE,f_lev%index)
    call call_hooks(pf, f_lev%index, PF_POST_INTERP_Q0)

  end subroutine interpolate_q0
  
    !>  Subroutine to update the fine terminal condition from coarse increment by spatial interpolation
    !>  used for adjoint solver
  subroutine interpolate_qend(pf, f_lev, c_lev)
  
    type(pf_pfasst_t), intent(inout) :: pf          !!  main pfasst structure
    class(pf_level_t),  intent(inout) :: f_lev  !!  fine level
    class(pf_level_t),  intent(inout) :: c_lev  !!  coarse level
    
    call call_hooks(pf, f_lev%index, PF_PRE_INTERP_Q0)
    if (pf%save_timings > 1) call pf_start_timer(pf, T_INTERPOLATE, f_lev%index)
    
    call c_lev%delta_q0%setval(0.0_pfdp)
    call f_lev%delta_q0%setval(0.0_pfdp)

    !>  restrict fine initial data to coarse
    call f_lev%ulevel%restrict(f_lev, c_lev, f_lev%qend, c_lev%delta_q0, pf%state%t0, flags=2)
    !>  get coarse level correction
    call c_lev%delta_q0%axpy(-1.0_pfdp, c_lev%qend, flags=2)    

    !>  interpolate correction in space
    call f_lev%ulevel%interpolate(f_lev, c_lev, f_lev%delta_q0, c_lev%delta_q0, pf%state%t0, flags=2)

    !> update fine inital condition
    call f_lev%qend%axpy(-1.0_pfdp, f_lev%delta_q0, flags=2)

    if (pf%save_timings > 1) call pf_stop_timer(pf, T_INTERPOLATE,f_lev%index)
    call call_hooks(pf, f_lev%index, PF_POST_INTERP_Q0)


  end subroutine interpolate_qend

end module pf_mod_interpolate
