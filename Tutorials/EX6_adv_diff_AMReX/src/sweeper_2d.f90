! MODULE: pf_my_sweeper
! IMEX Sweeper for 2-D  example in spectral space:
module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_AMReX_mfab
  use pf_mod_imex_sweeper
  use pf_mod_solutions
  use pf_mod_rutils
  use amrex_abeclaplacian_module
  use amrex_linop_module
  use amrex_multigrid_module
  
  implicit none

  !>  extend the IMEX sweeper
  type, extends(pf_imex_sweeper_t) :: my_sweeper_t
     integer ::     nx, ny
     !  Useful for making f_eval and f_comp generic
     class(pf_amrex_mfab_t), pointer :: f_encap,rhs_encap,y_encap 
     real(pfdp),        pointer :: yvec(:), rhsvec(:), fvec(:)
     
     ! AMReX MLMG Laplacian, solver and LinOP 
     type(amrex_abeclaplacian)  :: abeclap
     type(amrex_multigrid)      :: mlmg

   contains
     procedure :: f_eval       !  Computes the advection and diffusion terms - explicit part
     procedure :: f_comp       !  Computes the advection and diffusion terms - implicit part
     procedure :: initialize
     procedure :: destroy
  end type my_sweeper_t

contains
  

   subroutine advect_2d (lo, hi, y, f, plo, phi, dx)
      ! evaluate RHS for single patch
      use amrex_fort_module,  only : amrex_real
      use probin,             only: a, b        ! advective velocities
      ! in- / outputs
      integer, intent(in) :: lo(2), hi(2), plo(2), phi(2)
      real(amrex_real), intent(in   ) :: y(plo(1):phi(1), plo(2):phi(2))
      real(amrex_real), intent(inout) :: f(plo(1):phi(1), plo(2):phi(2))
      real(amrex_real), intent(in) :: dx(2)
      ! locals
      integer :: ix,iy
      real(amrex_real) :: dxinv, dyinv
      real(amrex_real) :: fx(lo(1):hi(1)+1,lo(2):hi(2)  )
      real(amrex_real) :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1)
      
      ! compute dx
      dxinv = 1.d0/dx(1)
      dyinv = 1.d0/dx(2)
    
      ! x-fluxes
      do iy=lo(2),hi(2)
         do ix=lo(1),hi(1)+1
            !fx(ix,iy) = ( y(ix-1,iy) - y(ix,iy) ) * 0.5_pfdp
            fx(ix,iy) = (-y(ix+1,iy) + 7.0_pfdp*(y(ix,iy) + y(ix-1,iy)) - y(ix-2,iy)) / 12.0_pfdp
         end do
      end do
    
      ! y-fluxes
      do iy=lo(2),hi(2)+1
         do ix=lo(1),hi(1)
            !fy(ix,iy) = ( y(ix,iy-1) - y(ix,iy) ) * 0.5_pfdp
            fy(ix,iy) = (-y(ix,iy+1) + 7.0_pfdp*(y(ix,iy) + y(ix,iy-1)) - y(ix,iy-2)) / 12.0_pfdp
         end do
      end do

      ! eval f
      do iy=lo(2),hi(2)
         do ix=lo(1),hi(1)
            f(ix,iy) = -a * dxinv * (fx(ix+1,iy  )-fx(ix,iy)) &
                     - b * dyinv * (fy(ix  ,iy+1)-fy(ix,iy))
         end do
      end do

   end subroutine advect_2d
  
   subroutine diff_2d (lo, hi, y, f, plo, phi, dx)
      ! evaluate RHS for single patch
      use amrex_fort_module,  only : amrex_real
      use probin,             only: nu
      ! in- / outputs
      integer, intent(in) :: lo(2), hi(2), plo(2), phi(2)
      real(amrex_real), intent(in   ) :: y(plo(1):phi(1), plo(2):phi(2))
      real(amrex_real), intent(inout) :: f(plo(1):phi(1), plo(2):phi(2))
      real(amrex_real), intent(in) :: dx(2)
      ! local
      integer :: ix, iy
      real(amrex_real) :: dxinv, dyinv
      real(amrex_real) :: fx(lo(1):hi(1)+1,lo(2):hi(2)  )
      real(amrex_real) :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1)
    
      !
      dxinv = 1.d0/dx(1)
      dyinv = 1.d0/dx(2)
    
      ! x-fluxes
      do iy=lo(2),hi(2)
         do ix=lo(1),hi(1)+1
            fx(ix,iy) = ( y(ix,iy) - y(ix-1,iy) ) * dxinv
            !fx(ix,iy) = (15.0_pfdp*(y(ix,iy) - y(ix-1,iy)) + y(ix-2,iy) - y(ix+1,iy)) * dxinv / 12.0_pfdp 
         end do
      end do

      ! y-fluxes
      do iy=lo(2),hi(2)+1
         do ix=lo(1),hi(1)
            fy(ix,iy) = ( y(ix,iy) - y(ix,iy-1) ) * dyinv
            !fy(ix,iy) = (15.0_pfdp*(y(ix,iy) - y(ix,iy-1)) + y(ix,iy-2) - y(ix,iy+1)) * dyinv / 12.0_pfdp 
         end do
      end do

      ! eval f - value is divergence of flux
      do iy=lo(2),hi(2)
         do ix=lo(1),hi(1)      
            f(ix,iy) = nu * ( (fx(ix+1,iy) - fx(ix,iy)) * dxinv &
                            + (fy(ix,iy+1) - fy(ix,iy)) * dyinv)
         end do
      end do

   end subroutine diff_2d

  ! The rest of the stuff is dimension independent
  ! include 'sweeper_include.f90' 
  ! somehow this does not work and since f_eval is dimension-dependent
  ! omit sweeper_include.f90 by adding functionaliy into each sweeper_$(Ndim)d.f90
  !

   !>  Helper function to return sweeper pointer
   function as_my_sweeper(sweeper) result(r)
      class(pf_sweeper_t), intent(inout), target :: sweeper
      class(my_sweeper_t), pointer :: r
      select type(sweeper)
      type is (my_sweeper_t)
         r => sweeper
      class default
         stop
      end select
   end function as_my_sweeper

   !>  Routine to set up sweeper variables and operators
   subroutine initialize(this, pf,level_index)
      use probin, only:  Ndim, splitting, mlmg_max_iter, linop_maxorder, linop_bc_lo, linop_bc_hi
      class(my_sweeper_t), intent(inout) :: this
      type(pf_pfasst_t),   intent(inout),target :: pf
      integer, intent(in) :: level_index
      ! local
      class(pf_level_t), pointer    :: lev
      real(8)                       :: a_scalar, b_scalar
      type(amrex_multifab)          :: alpha, beta(Ndim)
      
      !>  get level + shape
      lev => pf%levels(level_index)
      this%nx=lev%lev_shape(1)
      this%ny=lev%lev_shape(2)
      
      !> call superclass initialize
      call this%imex_initialize(pf,level_index)
      ! set IMEX options
      if (splitting .eq. 0 ) then     ! do implicit only
         this%implicit=.TRUE.
         this%explicit=.FALSE.
         print *, 'sweeper_2d.f90: fully implicit detected'
      elseif (splitting .eq. 1 ) then ! IMEX
         this%implicit=.TRUE.
         this%explicit=.TRUE.
         print *, 'sweeper_2d.f90: IMEX detected'
      elseif (splitting .eq. 2 ) then ! do explicit only 
         this%implicit=.FALSE.
         this%explicit=.TRUE.
         print *, 'sweeper_2d.f90: fully explicit detected'
      else                            ! do IMEX otherwise
         this%implicit=.TRUE.
         this%explicit=.TRUE.
      end if
      
      
      !> Initialize AMReX Laplacian
      ! inputs: abeclap, geom, ba, dm, metric_term, agglomeration, consolidation, max_coarsening_level 
      ! optional:   - metric_term           :: logical - default:  0 (in path_to_amrex/Src/LinearSolvers/MLMG/AMReX_MLLinOp.H = true)
      !             - agglomeration         :: logical - default:  0 (in path_to_amrex/Src/LinearSolvers/MLMG/AMReX_MLLinOp.H = true)
      !             - consolidation         :: logical - default:  0 (in path_to_amrex/Src/LinearSolvers/MLMG/AMReX_MLLinOp.H = true)
      !             - max_coarsening_level  :: integer - default: 30 (in path_to_amrex/Src/LinearSolvers/MLMG/AMReX_MLLinOp.H =   30)
      print *, 'Initialize AMReX MLMG Laplacian'
      select type (fac => lev%ulevel%factory)
        type is (pf_amrex_mfab_factory_t)
        call amrex_abeclaplacian_build(this%abeclap, [fac%geom], [fac%ba], [fac%dm])
      end select
      print *, 'Finished initalizing AMReX MLMG Laplacian'

      !> set order of stencil
      !linop_maxorder = 2
      call this%abeclap%set_maxorder(linop_maxorder) 
      print *, 'Setted maxorder of Laplacian stencil to: ', linop_maxorder

      !> define and set boundary conditions
      ! build array of boundary conditions needed by MLABecLaplacian (LinOpBCType)
      ! Definition of LinOpBCType from path_to_amrex/Src/Boundary/AMReX_LO_BCTYPES.H:
      ! LinOpBCType :: Integer
      !     interior         =    0
      !     Dirichlet        =  101
      !     Neumann          =  102
      !     reflect_odd      =  103
      !     Marshak          =  104
      !     SanchezPomraning =  105
      !     inflow           =  106
      !     inhomogNeumann   =  107
      !     Robin            =  108
      !     symmetry         =  109
      !     Periodic         =  200
      !     bogus            = 1729
      !
      !linop_bc_lo(Ndim), linop_bc_hi(Ndim) defined in inputfile, available via probin
      print *, 'LinOP BC lo, hi:', linop_bc_lo,',', linop_bc_hi
      call this%abeclap%set_domain_bc(linop_bc_lo, linop_bc_hi)  ! set boundary conditions
      print *, 'Set boundary conditions of Laplacian'

      !> set coefficients for equation (a_scalar * \alpha + b_scalar\nabla  * \beta\nabla) u = f
      ! a_scalar, b_scalar are scalar constants, \alpha and \beta are scalar fields, u is the unknown, f the rhs
      ! set a_scalar, b_scalar
      a_scalar = 1.0       ! 1.0 as value just for initialization        
      b_scalar = 1.0       ! 1.0 as value just for initialization  
      call this%abeclap%set_scalars(a_scalar, b_scalar)
      print *, 'Initialized A, B coefficients of Laplacian'
      ! set \alpha - \alpha is cell centered value -> single mfab
      select type (fac => lev%ulevel%factory)
        type is (pf_amrex_mfab_factory_t)
        call amrex_multifab_build(alpha, fac%ba, fac%dm, 1, 0)
      end select
      call alpha%setval(1.0_pfdp)                 ! 1.0 as value just for initialization
      call this%abeclap%set_acoeffs(0, alpha)     ! inputs: (level, alpha)
      print *, 'Initialized alpha coefficients of Laplacian'
      
      !> build beta - array pf face centered multifab
      call amrex_multifab_build(beta(1), alpha%ba, alpha%dm, 1, 0, [.true., .false.])  ! build multifab with nodalized boxarray
      call beta(1)%setval(1.0_pfdp)               ! 1.0 as value just for initialization   
      call amrex_multifab_build(beta(2), alpha%ba, alpha%dm, 1, 0, [.false., .true.])  ! build multifab with nodalized boxarray
      call beta(2)%setval(1.0_pfdp)               ! 1.0 as value just for initialization   
      call this%abeclap%set_bcoeffs(0, beta)      ! inputs: (level, beta)
      print *, 'beta(1)'
      call amrex_print(beta(1)%ba)
      print *, 'beta(2)'
      call amrex_print(beta(2)%ba)
      print *, 'Initialized beta coefficients of Laplacian'
      
      !> build MLMG solver
      call amrex_multigrid_build(this%mlmg, this%abeclap)
      ! mlmg options 
      call this%mlmg%set_max_iter(mlmg_max_iter)  ! max number of iterations
      !call this%mlmg%set_max_fmg_iter()          ! max number of multigrid cycles before switching to V-cycle -> not used here
      !call this%mlmg%set_fixed_iter()            ! sets fixed number of iterations
      call this%mlmg%set_bottom_solver(1)         ! s :: int 
      print *, 'Use default biconjugate gradient stabilized method (bicgstab) as AMReX MLMG solver!'
        ! Defined in path_to_amrex/Src/F_Interfaces/LinearSolvers/AMReX_multigrid_fi.cpp
        ! s = [ 0 -> smoother (such as Gauss-Seidel)
        !       1 -> bicgstab (default)
        !       2 -> cg       (Conjugate Gradient)
        !       3 -> hypre    (Hypre solver - external solver)
        !       4 -> petsc    (only for cell centered multifabs)]
        ! currently AMReX F-Interfaces does not support bicgcg and cgbicg
      call this%mlmg%set_verbose(0)           ! [0, 1] - set verbosity of solver
      call this%mlmg%set_bottom_verbose(0)    ! [0, 1] - set verbosity of bottom-solver

      !    allocate(this%tmp(this%nx))

   end subroutine initialize

   ! DESTROY: calls superclass destroy and deallocates all array type parameters
   subroutine destroy(this,pf,level_index)
      class(my_sweeper_t), intent(inout) :: this
      type(pf_pfasst_t),   intent(inout),target :: pf
      integer,             intent(in)    :: level_index
      
      !> Destroy AMReX MLMG & Laplacian
      call amrex_abeclaplacian_destroy(this%abeclap)
      call amrex_multigrid_destroy(this%mlmg)
      
      !>  Call base sweeper destroy
      call this%imex_destroy(pf,level_index)
      
      !   deallocate(this%tmp)
   
   end subroutine destroy

   !F_EVAL: evaluates the rhs function f(t,y(t)) at y, t.
   subroutine f_eval(this, y, t, level_index, f,piece)
      use probin, only: eq_type, splitting
      use amrex_fort_module, only : amrex_real
      ! arguments
      class(my_sweeper_t), intent(inout) :: this
      class(pf_encap_t),   intent(in)    :: y
      real(pfdp),          intent(in)    :: t             ! not used here - passed for compatibility
      integer,             intent(in)    :: level_index   ! not used here - passed for compatibility
      class(pf_encap_t),   intent(inout) :: f
      integer,             intent(in   ) :: piece  !> 1: advection, 2: diffusion

      ! local variables
      integer :: plo(4), phi(4)
      type(amrex_box) :: bx
      type(amrex_mfiter) :: mfi
      real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: y_arr, f_arr
    
      !> encaps
      this%y_encap => cast_as_amrex_mfab(y) 
      this%f_encap => cast_as_amrex_mfab(f) 
      
      !> solves
      if (piece .eq. 1) then  !  Explicit solve of advection term
         ! This fills periodic ghost cells and ghost cells from neighboring grids
         call this%y_encap%mfab%fill_boundary(this%f_encap%geom)

         ! use AMReX MFIter to iterate over Multifabs (=subdomains)
         call amrex_mfiter_build(mfi, this%y_encap%mfab, tiling=.true.)

         do while (mfi%next())
            ! get box
            bx = mfi%tilebox()
        
            ! pointer to solution and rhs
            y_arr => this%y_encap%mfab%dataptr(mfi)
            f_arr => this%f_encap%mfab%dataptr(mfi)
      
            ! get array bounds for solution and RHS from total array
            plo = lbound(y_arr)
            phi = ubound(y_arr)
        
            ! eval advect_2d
            call advect_2d(bx%lo, bx%hi, y_arr, f_arr, plo, phi, this%y_encap%geom%dx)
         end do

         ! destroy
         call amrex_mfiter_destroy(mfi)

      else if (piece .eq. 2) then  !  Explicit solve of diffusion term
         ! This fills periodic ghost cells and ghost cells from neighboring grids
         call this%y_encap%mfab%fill_boundary(this%f_encap%geom)

         ! use AMReX MFIter to iterate over Multifabs (=subdomains)
         call amrex_mfiter_build(mfi, this%y_encap%mfab, tiling=.true.)

         do while (mfi%next())
            ! get box
            bx = mfi%tilebox()
        
            ! pointer to solution and rhs
            y_arr => this%y_encap%mfab%dataptr(mfi)
            f_arr => this%f_encap%mfab%dataptr(mfi)
      
            ! get array bounds for solution and RHS from total array
            plo = lbound(y_arr)
            phi = ubound(y_arr)
        
            ! eval diff_1d
            call diff_2d(bx%lo, bx%hi, y_arr, f_arr, plo, phi, this%y_encap%geom%dx)
         end do

         ! destroy
         call amrex_mfiter_destroy(mfi)

      else  
        call pf_stop(__FILE__,__LINE__,'Unknown case for piece in f_eval ')  
      endif

   end subroutine f_eval


   ! Solve for y and return f2 also.
   subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
      use probin, only: splitting, nu, mlmg_tol_rel, mlmg_tol_abs
      class(my_sweeper_t), intent(inout) :: this
      class(pf_encap_t),   intent(inout) :: y
      real(pfdp),          intent(in   ) :: t
      real(pfdp),          intent(in   ) :: dtq
      class(pf_encap_t),   intent(in   ) :: rhs
      integer,             intent(in   ) :: level_index     ! unused here - passed for compatibility (?)
      class(pf_encap_t),   intent(inout) :: f
      integer,             intent(in   ) :: piece !> 1: advection, 2: diffusion

      !>  local variables
      real(pfdp)        :: qij, invdtq
      real(amrex_real)  :: a_scalar
      real(amrex_real)  :: y_sol
      
      !> encaps
      this%y_encap => cast_as_amrex_mfab(y) 
      this%f_encap => cast_as_amrex_mfab(f) 
      this%rhs_encap => cast_as_amrex_mfab(rhs)

      !> sanity check - should not call f_comp in explicit only
      if (splitting .eq. 2)  then
         print *,'WARNING: We should not be calling fcomp for fully explicit!'
         call this%y_encap%copy(rhs)
         call this%f_encap%setval(0.0_pfdp)
         return
      end if

      !> define case to solve
      if (piece == 1) then        ! implicit advection
        call pf_stop(__FILE__,__LINE__,'f_comp for advection term not implemented')
      else if (piece == 2) then   ! implicit diffusion
        ! fill ghost cells of each grid - includes periodic domain boundaries
        call this%rhs_encap%mfab%fill_boundary(this%rhs_encap%geom)     ! does this take care of all BCs?
        ! fill_boundary final ref to path_to_amrex/Src/Base/AMReX_FabArray.H
        ! seems to fill all boundaries, including periodic ones if periodic is provided
        ! strange thing SDC example in AMReX calls rhs.FillBoundary(geom.periodicity()); and afterwards FillDomainBoundary(rhs, geom, bc);
        !
        if (amrex_is_all_periodic() .eqv. .false.) then
          ! Fill non-periodic physical boundaries
          print *, 'WARNING: non-periodic BC'
          ! in C++ code:
          ! FillDomainBoundary(rhs, geom, bc);
          ! bc = [bc_x_lo, bc_y_lo, bc_x_hi, bc_y_hi] for 2D with bc_x/y_lo/hi :: Integer see below
          ! 
          ! added fill_domain_boundary to interface in path_to_amrex/Src/F_Interfaces/Base/AMReX_multifab_mod.f90
          ! should work like this:
          !bc = [0, 0]
          !call this%rhs_encap%mfab%fill_domain_boundary(this%y_encap%geom, bc) 
          !
          ! AMReX_linop_mod.f90 also offers set_domain_bc(this::amrex_linop, lobc::integer(Ndim), hibc::integer(Ndim))
          ! amrex_abeclaplacian is child of amrex_linop -> are set during initialization
        end if

         !> set diffusion parameters
         qij = nu * dtq
         a_scalar = 1.0_pfdp
         call this%abeclap%set_scalars(a_scalar, qij)
        
         !> set boundary conditions
         call this%abeclap%set_level_bc(0, this%rhs_encap%mfab)
         call this%abeclap%set_level_bc(0, this%y_encap%mfab)

         !> mlmg solve
         ! relative error tolerance mlmg_tol_rel is hardcoded to at least 1e-16
         !   see -> https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html# 
         y_sol = this%mlmg%solve([this%y_encap%mfab], [this%rhs_encap%mfab], mlmg_tol_rel, mlmg_tol_abs)
         !call this%mlmg%get_fluxes(this%f_encap%mfab)
         !print *, 'AMReX MLMG fluxes computed'
         ! further functionality of mlmg:
         !call this%mlmg%get_grad_solutions(grad_sol)    ! init grad_sol before as type(amrex_multifab)
         !call this%mlmg%comp_residual(residual,this%y_encap%mfab,this%rhs_encap%mfab)   ! init residual before as type(amrex_multifab)

         !> compute f
         ! y-dtq*f(y,t) = rhs -> f(y,t) = 1/dtq * y + (-1/dtq) * rhs
         invdtq = 1 / dtq
         ! dstmfab%lincomb(): dstmfab -> a * mfab1 + b * mfab2
         ! inputs: lincomb(a,mfab1,fab1_comp,b,mfab2,mfab2_comp,dstmfab_comp,ncomp,nghost)
         call this%f_encap%mfab%lincomb(invdtq,this%y_encap%mfab,1,-invdtq,this%rhs_encap%mfab,1, &
              1, this%f_encap%ncomp,this%f_encap%nghost)
      else  
        call pf_stop(__FILE__,__LINE__,'Unknown case for piece in f_comp')
      end if

   end subroutine f_comp

end module pf_my_sweeper


