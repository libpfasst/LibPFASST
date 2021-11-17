! MODULE: pf_my_sweeper_ad
! !> @author
!> Tommaso Buvoli
!
! Last Modified: Dec 28, 2018
!
! Description!
! Exponential Sweeper for 1-D advection/diffusion example: u_t + v*u_x = nu*u_xx

module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_exp_sweeper
  use phi_mod
  use pf_mod_fftpackage
  implicit none


  !>  extend the exponential sweeper
  type, extends(pf_exp_t) :: ad_sweeper_t
     integer ::     nx
     ! phi storage and scratch space
     complex(pfdp), allocatable :: phi(:,:)
     complex(pfdp), allocatable :: swp_phi(:,:,:)
     complex(pfdp), allocatable :: res_phi(:,:,:)
     complex(pfdp), allocatable :: yhat(:) ! used in applyPhi to store Fourier transform of final result
     complex(pfdp), allocatable :: bhat(:) ! used in applyPhi to store Fourier transform of final result
     real(pfdp) :: h_swpPhi = real(0.0, pfdp)
     real(pfdp) :: h_resPhi = real(0.0, pfdp)
     ! fft object and differentiaton matrices
     type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: lap(:) ! Laplacian operators
     complex(pfdp), allocatable :: ddx(:) ! First derivative operators

   contains

     procedure :: f_eval       !  Computes the advection and diffusion terms
     procedure :: phib         !  Computes phi functions
     procedure :: swpPhib      !  Computes phi functions
     procedure :: resPhib      !  Computes phi functions
     procedure, private :: storePhi
     procedure, private :: applyPhi
     procedure :: initialize
     procedure :: destroy
  end type ad_sweeper_t

contains

    !>  Helper function to return sweeper pointer
    function as_ad_sweeper(sweeper) result(r)
        class(pf_sweeper_t), intent(inout), target :: sweeper
        class(ad_sweeper_t), pointer :: r
        select type(sweeper)
        type is (ad_sweeper_t)
        r => sweeper
        class default
        stop
        end select
    end function as_ad_sweeper

    !>  Routine to set up sweeper variables and operators
    subroutine initialize(this, pf,level_index)
      class(ad_sweeper_t), intent(inout) :: this
      type(pf_pfasst_t),   intent(inout),target :: pf
      integer, intent(in) :: level_index
      
      integer ::  nx,nnodes
      class(pf_level_t), pointer :: lev
        

      lev => pf%levels(level_index)
      this%nx=lev%lev_shape(1)
      nx=this%nx
      nnodes=lev%nnodes

      ! call superclass initialize
      call this%exp_initialize(pf,level_index)
      
      ! Decide if phi will be recomputed
      this%use_phib=.true.   !  False means they will not be recomputed 

      
      ! allocate space for phi functions
      allocate(this%phi(lev%nnodes+1, nx))
      allocate(this%swp_phi(lev%nnodes - 1, lev%nnodes+1, nx))
      allocate(this%res_phi(lev%nnodes - 1, lev%nnodes+1, nx))
      allocate(this%yhat(nx))
      allocate(this%bhat(nx))
      
      ! allocate fft & differentiation matrices
      allocate(this%fft_tool)
      call this%fft_tool%fft_setup([nx],1)
      allocate(this%lap(nx))
      allocate(this%ddx(nx))
      call this%fft_tool%make_lap(this%lap)
      call this%fft_tool%make_deriv(this%ddx)
    end subroutine initialize
      

     ! DESTROY: calls superclass destroy and deallocates all array type parameters
     subroutine destroy(this,pf,level_index)
       class(ad_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index

       call this%exp_destroy(pf,level_index)
       
       ! deallocate arrays
       deallocate(this%phi)
       deallocate(this%swp_phi)
       deallocate(this%res_phi)
       deallocate(this%lap)
       deallocate(this%ddx)
       deallocate(this%yhat)
       deallocate(this%bhat)
       call this%fft_tool%fft_destroy()
       deallocate(this%fft_tool)
       
    end subroutine destroy

  ! ==================================================================
  ! Required subroutines for any sweeper that inherits from pf_exp_t
  ! ==================================================================

     !F_EVAL: evaluates the nonlinear function N(t,y(t)) at y, t.
     subroutine f_eval(this, y, t, level_index, n)
        use probin, only:  splitting, v,nu

        ! arguments
        class(ad_sweeper_t), intent(inout) :: this
        class(pf_encap_t),   intent(in)    :: y
        real(pfdp),          intent(in)    :: t
        integer,             intent(in)    :: level_index
        class(pf_encap_t),   intent(inout) :: n

        ! local variables
        real(pfdp),         pointer :: yvec(:), nvec(:)
        type(pf_fft_t),     pointer :: fft

        yvec  => get_array1d(y)
        nvec  => get_array1d(n)
        fft   => this%fft_tool

        select case (splitting)
            case (1) ! u_xx exponential, u_x explicit
                call fft%conv(yvec,-v * this%ddx,nvec)
            case (2) ! Fully Exponential
                call n%setval(real(0.0, pfdp))
            case (3) ! u_x exponential, u_xx explicit
                call fft%conv(yvec,nu * this%lap,nvec)
        end select

     end subroutine f_eval

     ! PHIB: compute phi products without storing result (useful for variable timestepping)
     subroutine phib(this, t, h, b, y)
         use probin, only: nu, v, splitting
         ! arguments
         class(ad_sweeper_t), intent(inout) :: this
         real(pfdp),          intent(in)    :: t
         real(pfdp),          intent(in)    :: h
         class(pf_encap_t),   intent(in)    :: b(:)
         class(pf_encap_t),   intent(inout) :: y

         select case (splitting)
             case (1)
                 call phi(t * h * (nu * this%lap), size(b) - 1, this%phi) ! initialize phi
             case (2)
                 call phi(t * h * (nu * this%lap - v * this%ddx), size(b) - 1, this%phi) ! initialize phi
             case (3)
                 call phi(t * h * (- v * this%ddx), size(b) - 1, this%phi) ! initialize phi
         end select

         call this%applyPhi(t, h, this%phi, b, y)      ! compute phi products

     end subroutine phib

     ! SWPPHIB: compute phi products for sweeper and stores result
     subroutine swpPhib(this, j, h, b, y) 
         
         ! arguments
         class(ad_sweeper_t), intent(inout) :: this
         integer,             intent(in)    :: j
         real(pfdp),          intent(in)    :: h
         class(pf_encap_t),   intent(in)    :: b(:)
         class(pf_encap_t),   intent(inout) :: y

         ! compute and store phi
         if(h .NE. this%h_swpPhi) then
             call storePhi(this, this%eta, h, this%swp_phi)
             this%h_swpPhi = h
         end if

         ! compute phi products
         call this%applyPhi(this%eta(j), h, this%swp_phi(j,:,:), b, y)

     end subroutine swpPhib

     ! RESPHIB: compute phi products for residual computation and stores result
     subroutine resPhib(this, j, h, b, y)

          ! arguments
          class(ad_sweeper_t), intent(inout) :: this
          integer,             intent(in)    :: j
          real(pfdp),          intent(in)    :: h
          class(pf_encap_t),   intent(in)    :: b(:)
          class(pf_encap_t),   intent(inout) :: y

          ! compute and store phi
          if(h .NE. this%h_resPhi) then
              call storePhi(this, this%nodes(2:size(this%nodes)), h, this%res_phi)
              this%h_resPhi = h
          end if

          ! compute phi products
          call this%applyPhi(this%nodes(j+1), h, this%res_phi(j,:,:), b, y)

      end subroutine resPhib

  ! ==================================================================
  ! helper subroutines for resPhib and swpPhib
  ! ==================================================================

      ! STOREPHI: computes and stores phi functions phi_j(t(i) h L)
      subroutine storePhi(this, ts, h, phi_storage)

          use probin, only: nu, v, splitting
          ! arguments
          class(ad_sweeper_t), intent(inout) :: this
          real(pfdp), intent(in)  :: ts(:)
          real(pfdp), intent(in)  :: h
          complex(pfdp), intent(out) :: phi_storage(:,:,:)
          ! local variables
          integer :: i, n, m

          m = size(ts)         ! number of t to compute
          n = size(this%nodes) ! highest order phi function
          do i = 1, m
             select case (splitting)
                 case (1)
                     call phi_zvector(ts(i) * h * (nu * this%lap), n, phi_storage(i,:,:)) ! initialize phi
                 case (2)
                     call phi_zvector(ts(i) * h * (nu * this%lap - v * this%ddx), n, phi_storage(i,:,:)) ! initialize phi
                 case (3)
                     call phi_zvector(ts(i) * h * (- v * this%ddx), n, phi_storage(i,:,:)) ! initialize phi
             end select
          end do

      end subroutine storePhi

      ! APPLYPHI: compute phi_0(t h L ) b_1 + h \sum_{i = 1} phi_i(t h L) b_i
      subroutine applyPhi(this, t, h, phis, b, y)

          ! arguments
          class(ad_sweeper_t), intent(inout) :: this
          real(pfdp),          intent(in)    :: t
          real(pfdp),          intent(in)    :: h
          complex(pfdp),       intent(in)    :: phis(:,:)
          class(pf_encap_t),   intent(in)    :: b(:)
          class(pf_encap_t),   intent(inout) :: y
          ! local variables
          integer :: i, n
          real(pfdp), pointer :: bvec(:), yvec(:)
          type(pf_fft_t),     pointer :: fft

          n = size(b)
          yvec => get_array1d(y)

          ! prepare fft variables
          fft    => this%fft_tool

          this%yhat = CMPLX(0.0, 0.0, pfdp)
          do i = 1, n
             bvec => get_array1d(b(i))

             call fft%fft(bvec,this%bhat)
             if(i .eq. 1) then
                 this%yhat = phis(1,:) * this%bhat                              ! phi_0(hL) b(1)           [no h factor]
             else
                this%yhat =this%yhat+  t**(i-1) * h * phis(i,:) * this%bhat   ! t^{n-1} phi_i(hL) h b(n) [h factor included]
             end if
         end do
         call fft%ifft(this%yhat,yvec)

      end subroutine applyPhi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial(y_0)
    type(pf_ndarray_t), intent(inout) :: y_0
    call exact(0.0_pfdp, y_0%flatarray)
  end subroutine initial

  !> Routine to return the exact solution
  subroutine exact(t, yex)
    use probin, only: nprob,nu, v, t00, kfreq
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    integer    :: nx, i, ii, k,nbox
    real(pfdp) :: tol, x, t0,Dx, omega

    nx = size(yex)
    Dx = 1.0d0/dble(nx)

    if (nprob .eq. 1) then
       !  Using sin wave initial condition
       omega = two_pi*kfreq
       do i = 1, nx
          x = dble(i-1-nx/2)/dble(nx) - t*v 
          yex(i) = cos(omega*x)*exp(-omega*omega*nu*t)
       end do
    else  !  Use periodic image of Gaussians
       yex=0
       if (nu .gt. 0.0) then
          nbox = ceiling(sqrt(4.0*nu*(t00+t)*37.0d0))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = dble(i-1-nx/2)/dble(nx) - t*v + dble(k)
                yex(i) = yex(i) + sqrt(t00)/sqrt(t00+t)*exp(-x*x/(4.0*nu*(t00+t)))
             end do
          end do
       else
          nbox = ceiling(sqrt(37.0d0))  !  Decide how many periodic images
          do k = -nbox,nbox
             do i = 1, nx
                x = dble(i-1-nx/2)/dble(nx) - t*v  + dble(k)
                yex(i) = yex(i) + exp(-x*x)
             end do
          end do
       end if  ! nbox

    end if
  end subroutine exact


end module pf_my_sweeper


