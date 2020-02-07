!> Definition of the local sweeper extending the IMEX sweeper
module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_zndsysarray
  use pf_mod_imex_sweeper
  use pf_mod_fftpackage
  use pf_mod_solutions
  implicit none

  !>  extend the exponential sweeper
  type, extends(pf_imex_sweeper_t) :: my_sweeper_t
     integer ::     nx

     type(pf_fft_t), pointer :: fft_tool     ! fft routines
     complex(pfdp), allocatable :: lap(:,:)  ! Laplacian operator
     complex(pfdp), allocatable :: ilap(:,:) ! inverse Laplacian operator
     complex(pfdp), allocatable :: ddx(:,:)  ! First derivative operator in x
     complex(pfdp), allocatable :: ddy(:,:)  ! First derivative operator in y
     complex(pfdp), allocatable :: tmp(:,:)  ! used for operator in real space
     complex(pfdp), allocatable :: tmphat(:,:) ! used for operator in spectral space
     complex(pfdp), allocatable :: u(:,:) !  u-velocity in real space
     complex(pfdp), allocatable :: v(:,:) !  v-velocity
    
   contains

     procedure :: f_eval       !  Computes the advection and diffusion terms
     procedure :: f_comp       !  Computes the advection and diffusion terms
     procedure :: initialize
     procedure :: destroy
     procedure :: exact
     procedure :: set_ic
  end type my_sweeper_t

contains

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
      use probin, only:   v0,eq_type,Lx
      class(my_sweeper_t), intent(inout) :: this
      type(pf_pfasst_t),   intent(inout),target :: pf
      integer, intent(in) :: level_index
      
      integer ::  nx,nnodes
      class(pf_level_t), pointer :: lev

      lev => pf%levels(level_index)
      this%nx=lev%lev_shape(1)
      nx=this%nx
      nnodes=lev%nnodes

      ! call superclass initialize
      call this%imex_initialize(pf,level_index)
      
      !>  Set variables for explicit and implicit parts
      this%implicit=.TRUE.
      this%explicit=.TRUE.
      
      
      ! allocate fft & differentiation matrices
      allocate(this%fft_tool)
      call this%fft_tool%fft_setup([nx,nx],2,Lx)
      allocate(this%u(nx,nx))
      allocate(this%v(nx,nx))
      allocate(this%lap(nx,nx))
      allocate(this%ilap(nx,nx))
      allocate(this%ddx(nx,nx))
      allocate(this%ddy(nx,nx))
      allocate(this%tmp(nx,nx))
      allocate(this%tmphat(nx,nx))

      
      call this%fft_tool%make_deriv(this%ddx,1) !  First derivative
      call this%fft_tool%make_deriv(this%ddy,2) !  First derivative      
      call this%fft_tool%make_lap(this%lap)  !  Second derivative
            
    end subroutine initialize
      

     ! DESTROY: calls superclass destroy and deallocates all array type parameters
     subroutine destroy(this,pf,level_index)
       class(my_sweeper_t), intent(inout) :: this
       type(pf_pfasst_t),   intent(inout),target :: pf
       integer,             intent(in)    :: level_index

       call this%imex_destroy(pf,level_index)
       
       ! deallocate arrays
       deallocate(this%u)
       deallocate(this%v)
       deallocate(this%lap)
       deallocate(this%ilap)
       deallocate(this%ddx)
       deallocate(this%ddy)
       deallocate(this%tmp)
       deallocate(this%tmphat)
       
       call this%fft_tool%fft_destroy()
       deallocate(this%fft_tool)
       
    end subroutine destroy


     !F_EVAL: evaluates the nonlinear function N(t,y(t)) at y, t.
     subroutine f_eval(this, y, t, level_index, f,piece)
       use probin, only: nu,kappa,grav,v0,splitting
        ! arguments
        class(my_sweeper_t), intent(inout) :: this
        class(pf_encap_t),   intent(in)    :: y
        real(pfdp),          intent(in)    :: t
        integer,             intent(in)    :: level_index
        class(pf_encap_t),   intent(inout) :: f
        integer,             intent(in   ) :: piece
        
        ! local variables
        complex(pfdp),  pointer :: omega(:,:),rho(:,:)     
        complex(pfdp),  pointer :: f_omega(:,:),f_rho(:,:) 
        type(pf_fft_t),     pointer :: fft

        !  Assign local pointers
        omega    => get_array2d(y,1)
        rho  => get_array2d(y,2)
        f_omega    => get_array2d(f,1)
        f_rho  => get_array2d(f,2)
        
        fft => this%fft_tool

        !  Get velocities in spectral space  by Biot-Savart then ifft
        this%tmphat=this%ddy*this%ilap*omega   !  u in spectral space        
        call fft%ifft(this%tmphat,this%u)         !  u in real space        
        this%tmphat=-this%ddx*this%ilap*omega  !  v in spectral space        
        call fft%ifft(this%tmphat,this%v)         !  v in real space        
        if (piece .eq. 1) then  !  Explicit
           this%tmphat=this%ddx*omega            !  omega_x in spectral space
           call fft%ifft(this%tmphat,this%tmp)   !  omega_x in real space
           if (splitting .eq. 1) then
              this%tmp = -(this%u+v0(1))*this%tmp !  -(u+v0(1)).omega_x in real space
           else
              this%tmp = -(this%u)*this%tmp       !  -u.omega_x in real space
           end if
           
           call fft%fft(this%tmp,f_omega)        !  -(u+v0(1)).omega_x in spectral space
           
           this%tmphat=this%ddy*omega            !  omega_y in spectral space
           call fft%ifft(this%tmphat,this%tmp)   !  omega_y in spectral
           if (splitting .eq. 1) then           
              this%tmp = (this%v+v0(2))*this%tmp   !  -(v+v0(2)).omega_y in real space
           else
              this%tmp = (this%v)*this%tmp      !  -v.omega_y in real space
           end if
           
           call fft%fft(this%tmp,this%tmphat)    !  v.omega_y in spectral space
           f_omega=f_omega-this%tmphat           !  f = -u.omega_x-v.omega_y

           f_omega=f_omega-grav*this%ddy*rho     !  f = -u.omega_x-v.omega_y-g*rho_y
           
           this%tmphat=this%ddx*rho              !  rho_x in spectral space
           call fft%ifft(this%tmphat,this%tmp)   !  rho_x in real space
           if (splitting .eq. 1) then           
              this%tmp = (this%u+v0(1))*this%tmp    !  u.rho_x in real space
           else
              this%tmp = (this%u)*this%tmp    !  u.rho_x in real space
           end if
           
           call fft%fft(this%tmp,this%tmphat)    !  u.rho_x in spectral space
           f_rho=  -this%tmphat                  !  f = -u.rho_x
           
           this%tmphat=this%ddy*rho              !  rho_y in spectral space
           call fft%ifft(this%tmphat,this%tmp)   !  rho_y in real space
           if (splitting .eq. 1) then           
              this%tmp = (this%v+v0(2))*this%tmp    !  v.rho_y in real space
           else
              this%tmp = (this%v)*this%tmp    !  v.rho_y in real space
           end if
           
           call fft%fft(this%tmp,this%tmphat)    !  v.rho_y in spectral space
           f_rho = f_rho - this%tmphat           !  f = -u.rho_x-v.rho_y
        else  !  Implicit
           if (splitting .eq. 1) then                      
              f_omega = nu*this%lap*omega           !  f = nu*Lap*omega
              f_rho = kappa*this%lap*rho            !  f = kappa*Lap*rho
           else
              f_omega = (nu*this%lap-v0(1)*this%ddx-v0(2)*this%ddy)*omega 
              f_rho = (kappa*this%lap-v0(1)*this%ddx-v0(2)*this%ddy)*rho 
           end if

        endif

     end subroutine f_eval

     ! Solve for y and return f2 also.
     subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
       use probin, only:  nu,kappa,splitting,v0
       class(my_sweeper_t), intent(inout) :: this
       class(pf_encap_t),   intent(inout) :: y
       real(pfdp),          intent(in   ) :: t
       real(pfdp),          intent(in   ) :: dtq
       class(pf_encap_t),   intent(in   ) :: rhs
       integer,             intent(in   ) :: level_index
       class(pf_encap_t),   intent(inout) :: f
       integer,             intent(in   ) :: piece

       ! local variables
       complex(pfdp),  pointer :: omega(:,:),rho(:,:)     
       complex(pfdp),  pointer :: f_omega(:,:),f_rho(:,:) 
       complex(pfdp),  pointer :: rhs_omega(:,:),rhs_rho(:,:) 
       type(pf_fft_t),     pointer :: fft

       !  Assign local pointers
       omega    => get_array2d(y,1)
       rho  => get_array2d(y,2)
       f_omega    => get_array2d(f,1)
       f_rho  => get_array2d(f,2)
       rhs_omega    => get_array2d(rhs,1)
       rhs_rho  => get_array2d(rhs,2)
       
       fft => this%fft_tool
       if (piece == 2) then
          ! Apply the inverse operator in spectral space
          if (splitting .eq. 1) then
             omega=rhs_omega/(1.0_pfdp - nu*dtq*this%lap)
             rho=rhs_rho/(1.0_pfdp - kappa*dtq*this%lap)
          else
             omega=rhs_omega/(1.0_pfdp - dtq*(nu*this%lap-v0(1)*this%ddx-v0(2)*this%ddy))
             rho=rhs_rho/(1.0_pfdp - dtq*(kappa*this%lap-v0(1)*this%ddx-v0(2)*this%ddy))
          end if
       

          f_omega = (omega - rhs_omega) / dtq
          f_rho = (rho - rhs_rho) / dtq
       else
          print *,'Bad piece in f_comp ',piece
          call exit(0)
       end if
     end subroutine f_comp



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine set_ic(this,y_0)
    class(my_sweeper_t), intent(inout)    :: this
    type(pf_zndsysarray_t), intent(inout) :: y_0

    ! local variables
    complex(pfdp),  pointer :: omega(:,:),rho(:,:)     
    
    omega => get_array2d(y_0,1)
    rho => get_array2d(y_0,2)
    call this%exact(0.0_pfdp, omega, rho)
    
  end subroutine set_ic

  !> Routine to return the exact solution
  subroutine exact(this,t, omega_ex, rho_ex)
    use probin, only: nx
    class(my_sweeper_t), intent(inout) :: this    
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: omega_ex(:,:)
    complex(pfdp), intent(inout) :: rho_ex(:,:)

    !  Local variables
    type(pf_fft_t),     pointer :: fft

    !  Get exact solution in real space (Borrow the temp space in the sweeper)
    call exact_realspace(t,this%u,this%v)

    ! Get exact solution in spectral space
    fft => this%fft_tool
    call fft%fft(this%u,omega_ex)
    call fft%fft(this%v,rho_ex)

    
  end subroutine exact

  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, omega_ex,rho_ex)
    use probin, only: eq_type,v0,nu, kappa, t00, kfreq,Lx,ic_type
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(out) :: omega_ex(:,:)
    complex(pfdp), intent(out) :: rho_ex(:,:)

    complex(pfdp) ::  fct

    select case (ic_type)
    case (0) 
       call exact_ad_exp(t,rho_ex,kappa,v0,Lx)
       omega_ex=cmplx(0.0,0.0,pfdp)
    case (1) 
       call exact_ad_cos(t,rho_ex,kappa,v0,kfreq,Lx)
       omega_ex=cmplx(0.0,0.0,pfdp)
    case (2)  !  Advection diffusion on rho with Taylor-Green on u,v (grav should be 0)
       call exact_tg_vortex(t,omega_ex)
       call exact_ad_exp(t,rho_ex,kappa,v0,Lx)
    case DEFAULT
       print *,'Bad case for ic_type in exact_realspace ', ic_type
       call exit(0)
    end select
    
  end subroutine exact_realspace

  subroutine exact_tg_vortex(t, omega_ex)
    use probin, only: nu, v0,Lx
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(inout) :: omega_ex(:,:)

    integer    :: nx(2), i, j
    real(pfdp) ::  x,y,diff_fac

    nx = SHAPE(omega_ex)

    diff_fac=exp(-2.0_pfdp*(two_pi/Lx(1))*(two_pi/Lx(2))*nu*t)
    do j = 1, nx(2)
       y = dble(j-1)/dble(nx(2))*Lx(2) - t*v0(2) 
       do i = 1, nx(1)    
          x = dble(i-1)/dble(nx(1))*Lx(1) - t*v0(1) 
          omega_ex(i,j) = -2.0_pfdp*cos(two_pi/Lx(1)*x)*cos(two_pi/Lx(2)*y)*diff_fac
       end do
    end do

  end subroutine exact_tg_vortex



end module pf_my_sweeper


