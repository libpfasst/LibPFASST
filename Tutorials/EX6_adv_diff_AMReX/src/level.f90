!
! This file is part of LIBPFASST.
!
! Level specification for  amrex_mfab encapsulations

module pf_my_level
  use pf_mod_pfasst
  use pf_mod_AMReX_mfab
  use amrex_multifabutil_module     ! contains average_down, needed for restrict
  use amrex_fillpatch_module        ! contains fill_patch, needed for interpolate 
  use iso_c_binding

  implicit none

  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: my_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type my_level_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  These are the transfer functions that must be  provided for the level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>  Interpolate from coarse  level to fine
  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use probin, only : geom_bc_lo, geom_bc_hi
    use amrex_interpolater_module

    class(my_level_t), intent(inout)  :: this
    class(pf_level_t), intent(inout)  :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t), intent(inout)  :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   )  :: t
    integer, intent(in), optional     :: flags

    !> local 
    class(pf_amrex_mfab_t), pointer         :: mfab_f,mfab_c        ! pointer to encaps
    class(pf_AMReX_mfab_factory_t), pointer :: factory_f, factory_c
    integer     :: refine_ratio
    integer     :: nc
    integer, allocatable :: lo_bc(:,:), hi_bc(:,:)

    !> set pointers
    factory_c => cast_as_AMReX_factory(c_lev%ulevel%factory)
    factory_f => cast_as_AMReX_factory(f_lev%ulevel%factory)
    mfab_f => cast_as_amrex_mfab(f_vec)
    mfab_c => cast_as_amrex_mfab(c_vec)

    !> compute refinement ratio
    ! checked during initalization that refinement ratio is of integer-type and same for all dimensions 
    refine_ratio = mfab_f%arr_shape(1) / mfab_c%arr_shape(1)
    
    !> interpolate using amrex_fillcoarsepatch 
    if (mfab_c%ncomp .eq. 1) then       ! single component
      nc = 1
      allocate(lo_bc(amrex_spacedim, nc))
      allocate(hi_bc(amrex_spacedim, nc))
      lo_bc(:,nc) = geom_bc_lo
      hi_bc(:,nc) = geom_bc_hi
      ! input is: mfab_f, t_old_c, mfab_old_c, t_new_c, mfab_new_c, geom_c, fill_physbc, geom_f, fill_physbc, &
      !             t, src_comp(old), dst_comp(new), num_comp, refine_ratio, interpolation_mode, bc_lo, bc_hi
      ! since we only interpolate in space, we use mfab_c for mfab_old and mfab_new
      ! interpolation_mode are defined in amrex/Src/F_Interfaces/AmrCore/AMReX_interpolater_mod.f90
      !   summary of interpolation modes:
      !           amrex_interp_pc            = 0
      !           amrex_interp_node_bilinear = 1
      !           amrex_interp_cell_bilinear = 2
      !           amrex_interp_quadratic     = 3
      !           amrex_interp_lincc         = 4
      !           amrex_interp_cell_cons     = 5
      !           amrex_interp_protected     = 6
      !           amrex_interp_quartic       = 7
      !           amrex_interp_face_divfree  = 8
      !           amrex_interp_face_linear   = 9
      call amrex_fillcoarsepatch(mfab_f%mfab, t, mfab_c%mfab, t, mfab_c%mfab, &
            & factory_c%geom, fill_physbc, factory_f%geom, fill_physbc, &
            & t, nc, nc, nc, refine_ratio, amrex_interp_cell_cons, lo_bc, hi_bc)
      deallocate(lo_bc)
      deallocate(hi_bc)
    else                                ! multiple components
      allocate(lo_bc(amrex_spacedim, mfab_f%ncomp))
      allocate(hi_bc(amrex_spacedim, mfab_f%ncomp))
      do nc = 1, mfab_f%ncomp
        lo_bc(:,nc) = geom_bc_lo
        hi_bc(:,nc) = geom_bc_hi
      end do                                
      ! same as single just loop over components
      do nc = 1, mfab_f%ncomp
        call amrex_fillcoarsepatch(mfab_f%mfab, t, mfab_c%mfab, t, mfab_c%mfab, &
              & factory_c%geom, fill_physbc, factory_f%geom, fill_physbc, &
              & t, nc, nc, mfab_f%ncomp, refine_ratio, amrex_interp_cell_cons, lo_bc, hi_bc)
      end do
      deallocate(lo_bc)
      deallocate(hi_bc)
    end if
    
  end subroutine interpolate
  
  !> function needed for interpolation - adapted from amrex-tutorials/ExampleCodes/FortranInterface/Advection_F
  subroutine fill_physbc (pmf, scomp, ncomp, time, pgeom) bind(c)
    use amrex_geometry_module, only : amrex_is_all_periodic
    use amrex_filcc_module, only : amrex_filcc
    use probin, only : geom_bc_lo, geom_bc_hi
    type(c_ptr), value :: pmf, pgeom
    integer(c_int), value :: scomp, ncomp
    real(amrex_real), value :: time

    type(amrex_geometry) :: geom
    type(amrex_multifab) :: mf
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
    integer :: plo(4), phi(4), nc
    integer, allocatable :: lo_bc(:,:), hi_bc(:,:)

    if (.not. amrex_is_all_periodic()) then
       geom = pgeom
       mf = pmf

       allocate(lo_bc(amrex_spacedim, mf%nc))
       allocate(hi_bc(amrex_spacedim, mf%nc))
       do nc = 1, mf%nc
        lo_bc(:,nc) = geom_bc_lo
        hi_bc(:,nc) = geom_bc_hi
       end do

       call amrex_mfiter_build(mfi, mf, tiling=.false.)
       do while(mfi%next())
          p => mf%dataptr(mfi)
          if (.not. geom%domain%contains(p)) then ! part of this box is outside the domain
             plo = lbound(p)
             phi = ubound(p)
             call amrex_filcc(p, plo, phi,         & ! fortran array and bounds
                  geom%domain%lo, geom%domain%hi,  & ! index extent of whole problem domain
                  geom%dx,                         & ! cell size in real
                  geom%get_physical_location(plo), & ! physical location of lower left corner
                  lo_bc, hi_bc)            ! bc types for each component

             ! amrex_filcc doesn't fill EXT_DIR (see amrex_bc_types_module for a list of bc types
             ! In that case, the user needs to fill it.
          end if
       end do
       call amrex_mfiter_destroy(mfi)
       deallocate(lo_bc)
       deallocate(hi_bc)
    end if
  end subroutine fill_physbc

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use pf_my_sweeper, only: my_sweeper_t

    class(my_level_t), intent(inout)  :: this
    class(pf_level_t), intent(inout)  :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t), intent(inout)  :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   )  :: t             !<  time of solution
    integer, intent(in), optional     :: flags

    !> local 
    class(pf_amrex_mfab_t), pointer         :: mfab_f,mfab_c        ! pointer to encaps
    class(pf_AMReX_mfab_factory_t), pointer :: factory_f, factory_c
    integer     :: refine_ratio, nc
    
    !> set pointers
    factory_c => cast_as_AMReX_factory(c_lev%ulevel%factory)
    factory_f => cast_as_AMReX_factory(f_lev%ulevel%factory)
    mfab_f => cast_as_amrex_mfab(f_vec)
    mfab_c => cast_as_amrex_mfab(c_vec) 
    
    !> compute refinement ratio
    ! checked during initalization that refinement ratio is of integer-type and same for all dimensions 
    refine_ratio = mfab_f%arr_shape(1) / mfab_c%arr_shape(2)

    !> restrict using amrex_average_down
    if (mfab_c%ncomp .eq. 1) then       ! single component
      nc = mfab_c%ncomp
      call amrex_average_down(mfab_f%mfab, mfab_c%mfab, factory_f%geom, factory_c%geom, nc, nc, refine_ratio)
    else                                ! loop over components
      do nc = 1, mfab_c%ncomp
        call amrex_average_down(mfab_f%mfab, mfab_c%mfab, factory_f%geom, factory_c%geom, nc, mfab_c%ncomp, refine_ratio)
      end do
    end if
    
  end subroutine restrict

end module pf_my_level
