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
    use probin, only: nu, vx, vy, kfreqx, kfreqy, ic_type
    real(pfdp), intent(in)  :: t
    type(pf_amrex_mfab_t), intent(inout) :: y_exact
    !> local helper
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx    
    real(pfdp), allocatable :: y_exact_flat(:)
    integer :: lo(3), hi(3), ix, iy, k
    real(pfdp) :: Lx, Ly, x, y, uy

    ! allocate
    allocate(y_exact_flat(y_exact%max_dof_proc)) 
    call amrex_mfiter_build(mfi, y_exact%mfab, tiling=.true.)
    Lx = amrex_probhi(1) - amrex_problo(1)
    Ly = amrex_probhi(2) - amrex_problo(2)
    k = 0
    ! loop over boxes
    do while (mfi%next())
      bx = mfi%tilebox()
      lo = bx%lo          ! min and max index per dimension in c++ style
      hi = bx%hi          ! min and max index per dimension in c++ style
      do iy = lo(2), hi(2)
        y = iy * y_exact%geom%dx(2)
        uy=ad_cos_ex(t, y,nu,vy,kfreqy,Ly)
        do ix = lo(1), hi(1)
          k = k + 1
          x = ix * y_exact%geom%dx(1) 
          y_exact_flat(k) = ad_cos_ex(t, x,nu,vx,kfreqx,Lx)*uy
        end do
      end do
    end do  
    ! unpack into mfab & de-alloc
    call y_exact%unpack(y_exact_flat)
    deallocate(y_exact_flat)
    
  end subroutine exact

end module pf_mod_rutils
