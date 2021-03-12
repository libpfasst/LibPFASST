!!  Exponential RK Integrator Module
! =====================================================================================
! MODULE: pf_mod_exprkstepper
! !> @author
! Tommaso Buvoli
!
! Last Modified: Nov 4, 2019
!
!> Exponential RK module
!!
!! This module implements an exponential RK method of the form:
!!
!!      Y_1 = F(y_n)
!!      Y_i = exp(h c(i) L) y_n + \sum_{j=1}^{i-1} \alpha_{ij}(hL) N(c_j, Y_j)    i = 2, ..., s
!!
!!      y_{n+1} = exp(h L) y_n + \sum_{j=1}^s \beta_{j}(hL) N(c_j, Y_j)
!!
!! where the functions
!!
!!      alpha_{ij}(L) = \sum_{k=1}^m A(1,k,i,j) * \varphi_k( A(2,k,j,i) * L )
!!      beta_{j}(L)   = \sum_{k=1}^m b(1,k,j) * \varphi_k( b(2,k,j) * L)
!!
!! The information for the EXPRK A matrix is stored as a 4D array A(:,:,:,:) where A(:, :, j, i) specifies the linear
!! combination of phi functions that that comprise the (i+1,j)th entry in the A matrix. In particular:
!!
!!          1. size(A, 3) = size(A, 4) is equal to (s-1) - first stage is explicit and not represented
!!
!!          2. size(A, 2) returns the highest-order phi-function used by the method
!!
!!          3. A(:, k, j, i) is a vector of dimension 2 that cooresponds to the term A(1, k, j, i) * phi_{k-1}(A(2, k, j, i) hL)        
!!
!! Similarly, the b vector is stored a 3D array b(:,:,:) where b(:,:,i) contains the information for the ith entry of b      

module pf_mod_erkstepper
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  type, extends(pf_stepper_t), abstract :: pf_erk_stepper_t
    
     real(pfdp), allocatable :: A(:,:,:,:)     
     integer,   allocatable  :: AF(:,:) ! AF(j,i) > 0 iff A(:,:,j,i) is non empty. AF(i,j) represents the number of nonzero entries coefficients in Tablaeu for (i < a) and (i=a, j <=b)   
     real(pfdp), allocatable :: b(:,:,:)
     real(pfdp), allocatable :: c(:)
     real(pfdp), allocatable :: d(:, :)     
     integer :: nstages
     integer :: nnz_A

     class(pf_encap_t), allocatable  :: Y_stage       !!  Stage value
     class(pf_encap_t), allocatable  :: PFY           !!  temp storage variable for storing phi products with stage derivatives         
     class(pf_encap_t), allocatable  :: y_n           !!  Local y0
     class(pf_encap_t), allocatable  :: y_np1         !!  Local yend

     class(pf_encap_t), pointer      :: F(:,:)        !!  Pointer to F

   contains
     procedure(pf_f_eval), deferred :: f_eval
     procedure(pf_compA),  deferred :: compA
     procedure(pf_compb),  deferred :: compB
     procedure(pf_compd),  deferred :: compD
     procedure :: do_n_steps  => erk_do_n_steps
     procedure :: initialize  => erk_initialize
     procedure :: destroy     => erk_destroy
     procedure :: erk_initialize 
     procedure :: erk_destroy    
  end type pf_erk_stepper_t

  interface
 
    subroutine pf_compA(this, dt, i, j, F, val)
       import pf_erk_stepper_t, pf_encap_t, pfdp
        class(pf_erk_stepper_t), intent(inout) :: this
        real(pfdp), intent(in) :: dt
        integer, intent(in) :: i, j
        class(pf_encap_t), intent(in)    :: F(:,:)
        class(pf_encap_t), intent(inout) :: val
    end subroutine pf_compA

    subroutine pf_compB(this, dt, i, F, val)
       import pf_erk_stepper_t, pf_encap_t, pfdp
        class(pf_erk_stepper_t), intent(inout) :: this
        real(pfdp), intent(in) :: dt
        integer, intent(in) :: i
        class(pf_encap_t), intent(in) :: F(:,:)
        class(pf_encap_t), intent(inout) :: val
    end subroutine pf_compB

    subroutine pf_compD(this, dt, i, y_n, val)
       import pf_erk_stepper_t, pf_encap_t, pfdp
        class(pf_erk_stepper_t), intent(inout) :: this
        real(pfdp), intent(in) :: dt
        integer, intent(in) :: i
        class(pf_encap_t), intent(in) :: y_n
        class(pf_encap_t), intent(inout) :: val
    end subroutine pf_compD
 
    subroutine pf_f_eval(this,y, t, level_index, f)
       import pf_erk_stepper_t, pf_encap_t, pfdp
       class(pf_erk_stepper_t),   intent(inout) :: this
       class(pf_encap_t), intent(in   ) :: y
       real(pfdp),        intent(in   ) :: t
       integer,           intent(in   ) :: level_index
       class(pf_encap_t), intent(inout) :: f
    end subroutine pf_f_eval

  end interface

contains

  subroutine erk_initialize(this, pf, level_index)
    class(pf_erk_stepper_t), intent(inout)  :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    
    integer    :: nstages,npieces   !  Local copies for convenience
    integer    :: i, j   
    real(pfdp) :: gamma, delta
    type(pf_level_t), pointer  :: lev    !  Current level
    integer :: s, m ! stages, and highest order phi

    lev => pf%levels(level_index)   !  Assign level pointer

    select case (this%order)
        case(1) ! Exponential Euler
                s = 1 ! number of stages (including explicit first stage)
                m = 2 ! number of phi functions (including \varphi_0)

                allocate(this%A(2, m, s - 1, s - 1))
                allocate(this%b(2, m, s))
                allocate(this%d(2, s))
                allocate(this%c(s - 1))

                this%A = 0.0_pfdp
                this%b = 0.0_pfdp
                this%d = 0.0_pfdp
                this%c = 0.0_pfdp

                ! -- A Matrix ------------------------------------------------------------------------
               
                ! -- b Vector ------------------------------------------------------------------------
                ! \varphi_1(hL) * F(Y1)
                this%b(:, 2, 1)    = [1.0_pfdp,  1.0_pfdp]

                ! -- d Vector ------------------------------------------------------------------------
                this%d(:, 1)       = [1.0_pfdp,  1.0_pfdp] ! exp(h L) y_n

                ! -- c Vector ------------------------------------------------------------------------

        case(2) ! Eqn (22) from S. M. Cox and P. C. Matthews, "Exponential time differencing for stiff systems." 2002 

                s = 2 ! number of stages (including explicit first stage)
                m = 3 ! number of phi functions (including \varphi_0)
                
                allocate(this%A(2, m, s - 1, s - 1))
                allocate(this%b(2, m, s))
                allocate(this%d(2, s))
                allocate(this%c(s - 1)) 
                
                this%A = 0.0_pfdp
                this%b = 0.0_pfdp
                this%d = 0.0_pfdp
                this%c = 0.0_pfdp

                ! -- A Matrix ------------------------------------------------------------------------
                this%A(:, 2, 1, 1) = [1.0_pfdp,  1.0_pfdp] ! \varphi_1(hL)
                
                ! -- b Vector ------------------------------------------------------------------------
                ! (\varphi_1(hL) - \varphi_2(hL)) * F(Y1)
                this%b(:, 3, 1)    = [-1.0_pfdp, 1.0_pfdp]
                this%b(:, 2, 1)    = [1.0_pfdp,  1.0_pfdp]  
                
                ! \varphi_2(h L) * F(Y2)
                this%b(:, 3, 2)    = [1.0_pfdp,  1.0_pfdp]

                ! -- d Vector ------------------------------------------------------------------------
                this%d(:, 1)       = [1.0_pfdp,  1.0_pfdp] ! exp(h L) y_n
                this%d(:, 2)       = [1.0_pfdp,  1.0_pfdp] ! exp(h L) y_n

                ! -- c Vector ------------------------------------------------------------------------
                this%c(:)          = [1.0_pfdp]

        case(3) ! Eqn (23)-(25) from S. M. Cox and P. C. Matthews, "Exponential time differencing for stiff systems." 2002 

                s = 3 ! number of stages (including explicit first stage)
                m = 4 ! number of phi functions (including \varphi_0)

                allocate(this%A(2, m, s - 1, s - 1))
                allocate(this%b(2, m, s))
                allocate(this%d(2, s))
                allocate(this%c(s - 1)) 
                
                this%A = 0.0_pfdp
                this%b = 0.0_pfdp
                this%d = 0.0_pfdp
                this%c = 0.0_pfdp

                ! -- A Matrix ------------------------------------------------------------------------
                this%A(:, 2, 1, 1) = [0.5_pfdp,  0.5_pfdp] ! (1/2) \varphi_1(hL/2)
                this%A(:, 2, 1, 2) = [-1.0_pfdp, 1.0_pfdp] ! -1 \varphi_1(hL) 
                this%A(:, 2, 2, 2) = [2.0_pfdp,  1.0_pfdp] ! 2  \varphi_1(hL)

                ! -- b Vector ------------------------------------------------------------------------
                ! (\varphi_1(hL) - 3 \varphi_2(hL) + 4 \varphi_3(hL)) * F(Y1)
                this%b(:, 2, 1)    = [1.0_pfdp, 1.0_pfdp]
                this%b(:, 3, 1)    = [-3.0_pfdp, 1.0_pfdp] 
                this%b(:, 4, 1)    = [4.0_pfdp, 1.0_pfdp] 
                ! (4 \varphi_2(hL) - 8 \varphi_3(hL)) * F(Y2)
                this%b(:, 3, 2)    = [4.0_pfdp, 1.0_pfdp] 
                this%b(:, 4, 2)    = [-8.0_pfdp, 1.0_pfdp]
                ! (- \varphi_2(hL) + 4 \varphi_3(hL)) * F(Y3) 
                this%b(:, 3, 3)    = [-1.0_pfdp, 1.0_pfdp] 
                this%b(:, 4, 3)    = [4.0_pfdp, 1.0_pfdp]

                ! -- d Vector ------------------------------------------------------------------------
                this%d(:, 1)       = [1.0_pfdp,  0.5_pfdp] ! exp(h L / 2) y_n
                this%d(:, 2)       = [1.0_pfdp,  1.0_pfdp] ! exp(h L) y_n
                this%d(:, 3)       = [1.0_pfdp,  1.0_pfdp] ! exp(h L) y_n

                ! -- c Vector ------------------------------------------------------------------------
                this%c(:)          = [0.5_pfdp, 1.0_pfdp]



        case (4) ! Eqn (51) from S. Krogstad. "Generalized integrating factor methods for stiff PDEs", 2005

                s = 4 ! number of stages 
                m = 4 ! number of phi functions
                
                allocate(this%A(2, m, s - 1, s - 1))
                allocate(this%b(2, m, s))
                allocate(this%d(2, s))
                allocate(this%c(s)) 
                
                this%A = 0.0_pfdp
                this%b = 0.0_pfdp
                this%d = 0.0_pfdp
                this%c = 0.0_pfdp

                ! -- A Matrix ------------------------------------------------------------------------ 
                ! A(1,1) = (1/2) \varphi_1(h/2 L)
                this%A(:, 2, 1, 1) = [0.5_pfdp, 0.5_pfdp] 

                ! A(2,1) = (1/2) \varphi_1(h/2 L) - \varphi_2(h/2 L)
                this%A(:, 2, 1, 2) = [0.5_pfdp, 0.5_pfdp]
                this%A(:, 3, 1, 2) = [-1.0_pfdp,  0.5_pfdp] 
                
                ! A(2,2) = (1/2) \varphi_2(h/2 L)
                this%A(:, 3, 2, 2) = [1.0_pfdp, 0.5_pfdp]
                
                ! A(3, 1) = \varphi_1(h L) - 2 \varphi_2(h L)
                this%A(:, 2, 1, 3) = [1.0_pfdp, 1.0_pfdp]
                this%A(:, 3, 1, 3) = [-2.0_pfdp, 1.0_pfdp]

                ! A(3,3) = 2 \varphi_2(h L)
                this%A(:, 3, 3, 3) = [2.0_pfdp, 1.0_pfdp]

                ! -- b Vector ------------------------------------------------------------------------
                ! (\varphi_1(hL) - 3 \varphi_2(hL) + 4 \varphi_3(hL) * F(y_n)
                this%b(:, 2, 1) = (/ 1.0_pfdp,  1.0_pfdp /) 
                this%b(:, 3, 1) = (/ -3.0_pfdp, 1.0_pfdp /)
                this%b(:, 4, 1) = (/ 4.0_pfdp,  1.0_pfdp /) 

                ! (2 \varphi_2(hL) - 4 \varphi_3(hL)) * F(Y_1)
                this%b(:, 3, 2) = (/ 2.0_pfdp,  1.0_pfdp /)
                this%b(:, 4, 2) = (/ -4.0_pfdp, 1.0_pfdp /)

                ! (2 \varphi_2(hL) - 4 \varphi_3(hL)) * F(Y_2)
                this%b(:, 3, 3) = (/ 2.0_pfdp,  1.0_pfdp /)
                this%b(:, 4, 3) = (/ -4.0_pfdp, 1.0_pfdp /)

                ! (-1 \varphi_2(hL) - 4 \varphi_3(hL)) * F(Y_3)
                this%b(:, 3, 4) = (/ -1.0_pfdp, 1.0_pfdp /)
                this%b(:, 4, 4) = (/ 4.0_pfdp, 1.0_pfdp /)

                ! -- d Vector ------------------------------------------------------------------------
                this%d(:, 1) = [1.0_pfdp,  0.5_pfdp] ! exp(h/2 L) y_n
                this%d(:, 2) = [1.0_pfdp,  0.5_pfdp] ! exp(h/2 L) y_n
                this%d(:, 3) = [1.0_pfdp,  1.0_pfdp] ! exp(h L) y_n 
                this%d(:, 4) = [1.0_pfdp,  1.0_pfdp] ! exp(h L) y_n

                ! -- c Vector ------------------------------------------------------------------------
                this%c(:) = [0.5_pfdp, 0.5_pfdp, 1.0_pfdp]
    
        end select

        ! Form the Matrix A_flag
        allocate(this%AF(s - 1, s - 1))
        this%nnz_A = 0
        this%AF = 0
        do i = 1, s - 1
                do  j = 1, s - 1
                        if(count(pack(this%A(:,:,j,i), .true.) .ne. 0.0_pfdp) .ne. 0) then
                                this%nnz_A   = this%nnz_A + 1
                                this%AF(j,i) = this%nnz_A
                        endif 
                enddo
        enddo

!!$        do i = 1, s - 1
!!$         write (*,*) this%AF(:,i)
!!$        enddo



        npieces = 1
        nstages = s
        this%nstages = nstages
        this%npieces = npieces

        !  Store the info in the pf structure
        pf%rk_order(level_index)=this%order
        pf%rk_nstages(level_index)=this%nstages

        ! Allocate space for local variables
        call lev%ulevel%factory%create_single(this%y_n,    level_index,  lev%lev_shape)
        call lev%ulevel%factory%create_single(this%y_np1,  level_index,  lev%lev_shape)
        call lev%ulevel%factory%create_single(this%Y_stage,level_index,  lev%lev_shape)
        call lev%ulevel%factory%create_single(this%PFY,    level_index,  lev%lev_shape)
        call lev%ulevel%factory%create_array(lev%Frkflt,   nstages*npieces, level_index,  lev%lev_shape)
        do i = 1, nstages*npieces
                call lev%Frkflt(i)%setval(0.0_pfdp, 0)
        end do    
        this%F(1:nstages,1:npieces) => lev%Frkflt

  end subroutine erk_initialize

  subroutine erk_destroy(this, pf,level_index)
    class(pf_erk_stepper_t),   intent(inout) :: this
    type(pf_pfasst_t),  target,  intent(inout) :: pf
    integer,              intent(in)    :: level_index 

    type(pf_level_t), pointer  :: lev        !  Current level
    lev => pf%levels(level_index)   !  Assign level pointer

   call lev%ulevel%factory%destroy_single(this%y_np1)
   call lev%ulevel%factory%destroy_single(this%y_n)
   call lev%ulevel%factory%destroy_single(this%Y_stage)
   call lev%ulevel%factory%destroy_single(this%PFY)
   call lev%ulevel%factory%destroy_array(lev%Frkflt)
  end subroutine erk_destroy


  !> Perform N steps of ark on level level_index and set yend appropriately.
  subroutine erk_do_n_steps(this, pf, level_index, t0, y0,yend,big_dt, nsteps_rk)
    use pf_mod_timer
    use pf_mod_hooks
    
    class(pf_erk_stepper_t),   intent(inout) :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in   )         :: t0           !!  Time at start of time interval
    class(pf_encap_t), intent(in   )         :: y0           !!  Starting value
    class(pf_encap_t), intent(inout)         :: yend         !!  Final value
    real(pfdp),        intent(in   )         :: big_dt       !!  Size of time interval to integrato on
    integer,           intent(in)            :: level_index  !!  Level of the index to step on
    integer,           intent(in)            :: nsteps_rk    !!  Number of steps to use
    
    class(pf_level_t), pointer               :: lev          !!  Pointer to level level_index
    integer                                  :: i, j, n      !!  Loop counters
    real(pfdp)                               :: tn           !!  Time at beginning of RKstep
    real(pfdp)                               :: dt           !!  Size of each ark step
    real(pfdp)                               :: tc           !!  stage time
    
    lev => pf%levels(level_index)   !! Assign pointer to appropriate level
    dt = big_dt/real(nsteps_rk, pfdp)   ! Set the internal time step size based on the number of rk steps
    
    call this%y_n%copy(y0)
    tn = t0

!    print *,'rank: ', pf%rank,' doing ',nsteps_rk,' steps on level ',level_index
    do n = 1, nsteps_rk      ! Loop over time steps
       ! Reset initial condition
       if (n > 1) then
          call this%y_n%copy(this%y_np1)
       endif
       ! Loop over stage values
       call pf_start_timer(pf,T_FEVAL,level_index)       
       call this%f_eval(this%y_n, tn, level_index, this%F(1,1))
       call pf_stop_timer(pf,T_FEVAL,level_index)       
       do i = 1, this%nstages - 1
          call this%compD(dt, i, this%y_n, this%y_stage)
          do j = 1, i
             if(this%AF(j,i) .gt. 0) then
                call this%compA(dt, i, j, this%F, this%PFY)
                call this%y_stage%axpy(1.0_pfdp, this%PFY)
             endif
          enddo
          call pf_start_timer(pf,T_FEVAL,level_index)       
          call this%f_eval(this%y_stage, tn + dt * this%c(i), level_index, this%F(i+1,1))
          call pf_stop_timer(pf,T_FEVAL,level_index)       
       enddo
       call this%compD(dt, this%nstages, this%y_n, this%y_np1)
       do i = 1, this%nstages
          call this%compB(dt, i, this%F, this%PFY)
          call this%y_np1%axpy(1.0_pfdp, this%PFY)
       enddo
       tn = t0 + dt             
    end do ! End Loop over time steps

    call yend%copy(this%y_np1)

  end subroutine erk_do_n_steps
  
end module pf_mod_erkstepper
