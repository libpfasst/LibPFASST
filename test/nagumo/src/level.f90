module pf_my_level
  use pf_mod_dtype
  use pf_mod_ndarray_oc
  use pf_mod_restrict
  use pf_mod_imexQ_oc
  use probin
  use solutions
  use pf_mod_fftpackage
  
  implicit none


  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

contains

  subroutine restrict_control(sG, sF)
    use feval, only: ad_sweeper_t, as_ad_sweeper
    class(pf_sweeper_t), intent(inout) :: sG, sF

    real(pfdp), pointer  :: uF(:,:), uG(:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m

    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    uF => sweeperF%u
    uG => sweeperG%u

    nnodesF = size(uF,1)
    nnodesG = size(uG,1)
    nvarF = size(uF,2)
    nvarG = size(uG,2)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict u', xrat, trat

    !do m=1,nnodesG
       uG(:,:) = uF(::trat,::xrat)
    !end do
  end subroutine restrict_control


  subroutine restrict_ydesired(sG, sF)
    use feval, only: ad_sweeper_t, as_ad_sweeper
    class(pf_sweeper_t), intent(inout) :: sG, sF

    real(pfdp), pointer  :: ydesiredF(:,:), ydesiredG(:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m

    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    ydesiredF => sweeperF%ydesired
    ydesiredG => sweeperG%ydesired

    nnodesF = size(ydesiredF,1)
    nnodesG = size(ydesiredG,1)
    nvarF = size(ydesiredF,2)
    nvarG = size(ydesiredG,2)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict ydesired', xrat, trat

    !do m=1,nnodesG
       ydesiredG(:,:) = ydesiredF(::trat,::xrat)
    !end do
  end subroutine restrict_ydesired


  subroutine restrict_for_adjoint(pf, which)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: which
    real(pfdp), pointer :: tF(:), tG(:) !zF(:), zG(:)
    integer :: l, m !, nnodesF, nnodesG, nvarsF, nvarsG

    do l = pf%nlevels, 2, -1
      allocate(tF(pf%levels(l)%nnodes))
      allocate(tG(pf%levels(l-1)%nnodes))
      tF = pf%state%t0 + pf%state%dt*pf%levels(l)%nodes
      tG = pf%state%t0 + pf%state%dt*pf%levels(l-1)%nodes
      call restrict_ts(pf%levels(l), pf%levels(l-1), pf%levels(l)%Q, pf%levels(l-1)%Q,tF, which)

        call pf%levels(l-1)%ulevel%sweeper%evaluate_all(pf,l-1, tG, which)
      deallocate(tF)
      deallocate(tG)
    end do
  end subroutine restrict_for_adjoint


  subroutine interp1(qF, qC, adF, adC)
    use feval, only: ad_sweeper_t, as_ad_sweeper
    class(ad_sweeper_t), pointer :: adF, adC
    real(pfdp),  pointer, intent(inout) :: qF(:), qC(:)

    class(ad_sweeper_t), pointer :: sweeperF, sweeperC
    complex(pfdp),       pointer :: wkF(:), wkC(:)
    type(pf_fft_t),      pointer :: fftF,fftC

    integer      :: NxF, NxC, xrat,i

    sweeperC => as_ad_sweeper(adC)
    sweeperF => as_ad_sweeper(adF)
    fftC => sweeperC%fft_tool
    fftF => sweeperF%fft_tool
    
    NxF = size(qF)
    NxC = size(qC)
    xrat  = NxF / NxC

    if (xrat == 1) then
       qF = qC
       return
    endif

    call fftC%interp(qC,fftF,qF)

  end subroutine interp1


  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use feval, only: ad_sweeper_t, as_ad_sweeper
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    integer :: which
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),          pointer :: f(:), g(:)

    which = 0
    if(present(flags)) which = flags

    adG => as_ad_sweeper(c_lev%ulevel%sweeper)
    adF => as_ad_sweeper(f_lev%ulevel%sweeper)

    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array1d_oc(f_vec,1)
      g => get_array1d_oc(c_vec,1)
      call interp1(f, g, adF, adG)
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array1d_oc(f_vec,2)
      g => get_array1d_oc(c_vec,2)
      call interp1(f, g, adF, adG)
    end if
  end subroutine interpolate


  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    real(pfdp), pointer :: f(:), g(:)

    integer :: nvarF, nvarG, xrat, which
    which = 0
    if(present(flags)) which = flags

    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array1d_oc(f_vec,1)
      g => get_array1d_oc(c_vec,1)
      nvarF = size(f)
      nvarG = size(g)
      xrat  = nvarF / nvarG
      g = f(::xrat)
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array1d_oc(f_vec,2)
      g => get_array1d_oc(c_vec,2)
      nvarF = size(f)
      nvarG = size(g)
      xrat  = nvarF / nvarG
      g = f(::xrat)
    end if
  end subroutine restrict

end module pf_my_level

