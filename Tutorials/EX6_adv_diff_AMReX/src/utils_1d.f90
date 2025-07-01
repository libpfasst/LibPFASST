module pf_mod_rutils
  use pf_mod_solutions
  use pf_mod_dtype
  use pf_mod_AMReX_mfab
  implicit none
  
contains

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines currently only for adv-diff   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Routine to return the exact solution
  subroutine exact(t, y_exact)  
    real(pfdp), intent(in)  :: t
    type(pf_amrex_mfab_t), intent(inout) :: y_exact
    !> local helper
    real(pfdp), allocatable :: yex(:)
    allocate(yex(y_exact%arr_shape(1))) ! use either y_exact%ndof or y_exact%arr_shape - y_exact%pack_size contains ghost cells 

    ! compute and unpack flatarray into mfab
    call exact_realspace(t,yex)
    call y_exact%unpack(yex)
  end subroutine exact

  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, yex)
    use probin, only: nu, a, kfreq, Lx, ic_type
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)
    
    yex=0.0_pfdp
    
    if (ic_type .eq. 1) then
      call exact_ad_cos(t,yex,nu,a,kfreq,Lx)
    else
      call exact_ad_exp(t,yex,nu,a,Lx)
    end if
   end subroutine exact_realspace  
  
end module pf_mod_rutils
