module initial
  use encap
  implicit none
  include 'fftw3.f03'

  type :: amplitude
     double precision :: ax, ay, az
     integer :: kx, ky, kz
  end type amplitude
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vortex_sheets(q0)
    type(carray4), intent(inout) :: q0

    double precision :: x, y, z
    integer :: ii, i, j, k, n

    type(c_ptr) :: ffft, wkp

    complex(c_double), pointer :: u(:,:,:), v(:,:,:), w(:,:,:), wk(:,:,:)
    double precision, parameter :: delta = 0.1d0, rho=50.0d0

    n = q0%shape(1)

    allocate(u(n,n,n), v(n,n,n), w(n,n,n))

    u = -1.0d0
    w =  0.0d0

    do i = 1, n
       x = dble(i) / n
       do j = 1, n
          y = dble(j) / n
          do k = 1, n
             z = dble(k) / n

             do ii = -4, 4
                u(k, j, i) = u(k, j, i) + tanh(rho*(ii + y - 0.25d0))
                u(k, j, i) = u(k, j, i) + tanh(rho*(ii + 0.75d0 - y))
             end do

             v(k, j, i) = delta * sin(6.28318530718d0*(x + 0.25d0))

          end do
       end do
    end do

    wkp  = fftw_alloc_complex(int(n**3, c_size_t))
    ffft = fftw_plan_dft_3d(n, n, n, wk, wk, FFTW_FORWARD, FFTW_ESTIMATE)

    call c_f_pointer(wkp, wk, [ n, n, n ])

    wk = u
    call fftw_execute_dft(ffft, wk, wk)
    q0%array(:,:,:,1) = wk

    wk = v
    call fftw_execute_dft(ffft, wk, wk)
    q0%array(:,:,:,2) = wk

    wk = w
    call fftw_execute_dft(ffft, wk, wk)
    q0%array(:,:,:,3) = wk

    q0%array = q0%array / n**3

    deallocate(u,v,w,wk)

    call fftw_destroy_plan(ffft)

  end subroutine vortex_sheets

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine random_full(q0)
    type(carray4), intent(inout) :: q0

    type(amplitude), allocatable :: amps(:)
    integer :: namps

    double precision :: kappa, scale
    integer :: kx, ky, kz, i, j, k, n, a

    double precision :: x, y, z, r(3)
    double precision, parameter :: sigma = 99.0d0

    complex(c_double), pointer :: u(:,:,:), v(:,:,:), w(:,:,:), wk(:,:,:)

    type(c_ptr) :: ffft, wkp

    n = q0%shape(1)

    allocate(u(n,n,n), v(n,n,n), w(n,n,n))
    u = 0; v = 0; w = 0

    namps = 0
    allocate(amps(n*n*n/8))

    do kz = 1, n/2
       do ky = 1, n/2
          do kx = 1, n/2
             if (kx + ky + kz > n/2) cycle

             kappa = dble(kx**2 + ky**2 + kz**2)**0.5
             scale = kappa**(-5.0/3) * exp(-kappa**2/sigma)

             call random_number(r)

             namps = namps + 1
             amps(namps)%kx = kx
             amps(namps)%ky = ky
             amps(namps)%kz = kz
             amps(namps)%ax = r(1) * scale
             amps(namps)%ay = r(2) * scale
             amps(namps)%az = r(3) * scale
          end do
       end do
    end do

    !$omp parallel do private(i,j,k)
    do i = 1, n
       x = 6.28318530718d0 * dble(i-1) / n
       do j = 1, n
          y = 6.28318530718d0 * dble(j-1) / n
          do k = 1, n
             z = 6.28318530718d0 * dble(k-1) / n
             do a = 1, namps
                u(k,j,i) = u(k,j,i) + amps(a)%ax * cos(amps(a)%ky * y) * cos(amps(a)%kz * z)
                v(k,j,i) = v(k,j,i) + amps(a)%ay * cos(amps(a)%kx * x) * cos(amps(a)%kz * z)
                w(k,j,i) = w(k,j,i) + amps(a)%az * cos(amps(a)%kx * x) * cos(amps(a)%ky * y)
             end do
          end do
       end do
    end do
    !$omp end parallel do

    wkp  = fftw_alloc_complex(int(n**3, c_size_t))
    ffft = fftw_plan_dft_3d(n, n, n, wk, wk, FFTW_FORWARD, FFTW_ESTIMATE)

    call c_f_pointer(wkp, wk, [ n, n, n ])

    wk = u
    call fftw_execute_dft(ffft, wk, wk)
    q0%array(:,:,:,1) = wk

    wk = v
    call fftw_execute_dft(ffft, wk, wk)
    q0%array(:,:,:,2) = wk

    wk = w
    call fftw_execute_dft(ffft, wk, wk)
    q0%array(:,:,:,3) = wk

    q0%array = q0%array / n**3

    deallocate(u,v,w,wk,amps)

    call fftw_destroy_plan(ffft)

  end subroutine random_full

  subroutine load(q0, fname)
    use hdf5
    implicit none
    type(carray4),    intent(inout) :: q0
    character(len=*), intent(in   ) :: fname

    integer(hid_t)   :: plist, ctype, file, group, dataset
    integer(size_t)  :: re_size
    integer(hsize_t) :: dims(4)
    integer          :: err

    real(8), pointer    :: buf(:)
    complex(8), pointer :: cbuf(:,:,:,:)

    type complex_t
       double precision :: re
       double precision :: im
    end type complex_t

    call h5open_f(err)

    ! create 'complex' compound type
    call h5tget_size_f(H5T_NATIVE_DOUBLE, re_size, err)
    call h5tcreate_f(H5T_COMPOUND_F, 2*re_size, ctype, err)
    call h5tinsert_f(ctype, "r", int(0, size_t), H5T_NATIVE_DOUBLE, err)
    call h5tinsert_f(ctype, "i", int(re_size, size_t), H5T_NATIVE_DOUBLE, err)

    ! open the initial condition file and root group
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist, err)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file, err, access_prp=plist)
    call h5gopen_f(file, "/", group, err)

    ! read 'q0' dataset
    dims = q0%shape
    allocate(buf(2*product(q0%shape)))

    call h5dopen_f(group, "q0", dataset, err)
    call h5dread_f(dataset, ctype, buf, dims, err)
    call h5dclose_f(dataset, err)

    call c_f_pointer(c_loc(buf(1)), cbuf, q0%shape)
    q0%array = cbuf

    deallocate(buf)

    ! tidy up
    call h5gclose_f(group, err)
    call h5fclose_f(file, err)
    call h5tclose_f(ctype, err)
    call h5pclose_f(plist, err)
    call h5close_f(err)
  end subroutine load


end module initial
