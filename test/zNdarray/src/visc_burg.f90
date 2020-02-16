  subroutine ex_visc_burg(this,yreal,t,fft,nx)
    use probin, only: nu, Lx
    class(my_sweeper_t), intent(inout) :: this
    complex(pfdp),intent(inout) :: yreal(:)
    type(pf_fft_t), intent(in), pointer :: fft    
    real(pfdp), intent(in)  :: t
    integer, intent(in)  :: nx

    real(pfdp) :: x,uint    
    complex(pfdp),pointer :: phi(:)
    complex(pfdp),pointer :: phihat(:)
    complex(pfdp),pointer  :: phix(:)
    integer    ::  i
        
    allocate(phi(nx))
    allocate(phihat(nx))
    allocate(phix(nx))

    
    !  Form phi(x,0)
    do i = 1, nx
       x = Lx*real(i-1-nx/2,pfdp)/real(nx,pfdp) 
       !          phi(i) = exp(-(cos(two_pi*x)-1.0_pfdp)/(two_pi*2.0_pfdp*nu))
       uint = -(cos(two_pi*x)-1.0_pfdp)/two_pi
       phi(i) = exp(-uint/(2.0_pfdp*nu))          
    end do

    !  Now solve the heat equation for phi
    call fft%fft(phi,phihat)
    if (t .gt. 0.0) then
       phihat=phihat*exp(this%lap*nu*t)
    end if
    
       
    ! Get phi back in real space
    call fft%ifft(phihat,phi)
    
    ! Get phi_x back in real space
    phihat=phihat*this%ddx
    call fft%ifft(phihat,phix)
    
    ! now the solution is the quotient
    do i = 1, nx
       yreal(i) = -2.0_pfdp*nu*phix(i)/phi(i)
    end do
    deallocate(phi,phihat,phix)       
  end subroutine ex_visc_burg
