!
! Copyright (c) 2015, Andreas Kreienbuehl and Michael Minion. All rights reserved.
!

! ------------------------------------------------------ module `encap`: start
! PFASST: Encapsulation for Dahlquist equation
module encap
	! For, e.g., C-pointers
	use iso_c_binding 

	! For LIBPFASST
	use pfasst

	! Variables starting with `i`, `j`, `k`, `l`, `m`, or `n` represent integers
	implicit none

	! Encapsulated data (e.g. that of the solution)
	type :: dat
		! Number of variables required to define the encapsulated data array
		integer :: nvars

		! Distribution of number of variables over the encapsulated data array
		integer :: shape(4)

		! Data array (e.g. the unknown solution array)
		real(pfdp), pointer :: u(:, :, :, :)

		! Field variables
		integer :: nfields

		! Computational domain
		integer :: nx, ny, nz
	end type dat

contains
	! ------------------------------------------------------ function `get_nvars`: start
	! Number of variables required to define the encapsulated data
	function get_nvars(Cptr2_dat) result(dat_nvars)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_nvars

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_nvars = Fptr2_dat%nvars
	end function get_nvars
	! ------------------------------------------------------ function `get_nvars`: stop

	! ------------------------------------------------------ function `get_shape`: start
	! Distribution of number of variables over the encapsulated data array
	function get_shape(Cptr2_dat) result(dat_shape)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_shape(4)

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_shape = Fptr2_dat%shape
	end function get_shape
	! ------------------------------------------------------ function `get_shape`: stop

	! ------------------------------------------------------ function `get_u`: start
	! Data array (e.g. the unknown solution array)
	function get_u(Cptr2_dat) result(dat_u)
		type(c_ptr), intent(in), value :: Cptr2_dat
		real(pfdp), pointer :: dat_u(:, :, :, :)

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_u => Fptr2_dat%u
	end function get_u
	! ------------------------------------------------------ function `get_u`: stop

	! ------------------------------------------------------ function `get_nfields`: start
	! Computational domain
	function get_nfields(Cptr2_dat) result(dat_nfields)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_nfields

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_nfields = Fptr2_dat%nfields
	end function get_nfields
	! ------------------------------------------------------ function `get_nfields`: stop

	! ------------------------------------------------------ function `get_nx`: start
	! Computational domain
	function get_nx(Cptr2_dat) result(dat_nx)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_nx

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_nx = Fptr2_dat%nx
	end function get_nx
	! ------------------------------------------------------ function `get_nx`: stop

	! ------------------------------------------------------ function `get_ny`: start
	! Computational domain
	function get_ny(Cptr2_dat) result(dat_ny)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_ny

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_ny = Fptr2_dat%ny
	end function get_ny
	! ------------------------------------------------------ function `get_ny`: stop

	! ------------------------------------------------------ function `get_nz`: start
	! Computational domain
	function get_nz(Cptr2_dat) result(dat_nz)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_nz

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_nz = Fptr2_dat%nz
	end function get_nz
	! ------------------------------------------------------ function `get_nz`: stop

	! ------------------------------------------------------ subroutine `create_dat`: start
	! Creation of encapsulated data
	subroutine create_dat(Cptr2_dat, lvl, kind, nvars, shape, Cptr2_ctx)
		type(c_ptr), intent(inout) :: Cptr2_dat
		integer, intent(in) :: lvl, kind, nvars, shape(:)
		type(c_ptr), intent(in), value :: Cptr2_ctx

		type(dat), pointer :: Fptr2_dat
		integer :: ierr

		allocate(Fptr2_dat, STAT = ierr)
		if(ierr .ne. 0) stop 'Allocation of `Fptr2_dat` in `create_dat` failed'

		Fptr2_dat%nvars = nvars

		Fptr2_dat%shape = shape

		allocate(Fptr2_dat%u(Fptr2_dat%shape(1), Fptr2_dat%shape(2), Fptr2_dat%shape(3), Fptr2_dat%shape(4)), STAT = ierr)
		if(ierr .ne. 0) stop 'Allocation of `Fptr2_dat%u` in `create_dat` failed'

		Fptr2_dat%nfields = Fptr2_dat%shape(1)

		Fptr2_dat%nx = Fptr2_dat%shape(2)
		Fptr2_dat%ny = Fptr2_dat%shape(3)
		Fptr2_dat%nz = Fptr2_dat%shape(4)

		Cptr2_dat = c_loc(Fptr2_dat)
	end subroutine create_dat
	! ------------------------------------------------------ subroutine `create_dat`: stop

	! ------------------------------------------------------ subroutine `destroy_dat`: start
	! Destruction of encapsulated data
	subroutine destroy_dat(Cptr2_dat)
		type(c_ptr), intent(in), value :: Cptr2_dat

		type(dat), pointer :: Fptr2_dat
		integer :: ierr

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		deallocate(Fptr2_dat%u, STAT = ierr)
		if(ierr .ne. 0) stop 'Deallocation of `Fptr2_dat%u` in `destroy_dat` failed'

		deallocate(Fptr2_dat, STAT = ierr)
		if(ierr .ne. 0) stop 'Deallocation of `Fptr2_dat` in `destroy_dat` failed'
	end subroutine destroy_dat
	! ------------------------------------------------------ subroutine `destroy_dat`: stop

	! ------------------------------------------------------ subroutine `setval_dat`: start
	! Assignment of scalar value to all elements in encapsulated data array
	subroutine setval_dat(Cptr2_dst_dat, src_dat_u, flags)
		type(c_ptr), intent(in), value :: Cptr2_dst_dat
		real(pfdp), intent(in) :: src_dat_u
		integer, intent(in), optional :: flags

		real(pfdp), pointer :: dst_dat_u(:, :, :, :)
		integer :: which

		which = 0
		if(present(flags)) which = flags

		dst_dat_u => get_u(Cptr2_dst_dat)

		dst_dat_u = src_dat_u
	end subroutine setval_dat
	! ------------------------------------------------------ subroutine `setval_dat`: stop

	! ------------------------------------------------------ subroutine `copy_dat`: start
	! Copy of encapsulated data array
	subroutine copy_dat(Cptr2_dst_dat, Cptr2_src_dat, flags)
		type(c_ptr), intent(in), value :: Cptr2_dst_dat, Cptr2_src_dat
		integer, intent(in), optional :: flags

		real(pfdp), pointer :: dst_dat_u(:, :, :, :), src_dat_u(:, :, :, :)
		integer :: which

		which = 0
		if (present(flags)) which = flags

		dst_dat_u => get_u(Cptr2_dst_dat)
		src_dat_u => get_u(Cptr2_src_dat)

		dst_dat_u = src_dat_u
	end subroutine copy_dat
	! ------------------------------------------------------ subroutine `copy_dat`: stop

	! ------------------------------------------------------ function `norm_dat`: start
	! Norm of encapsulated data array
	function norm_dat(Cptr2_dat) result(norm)
		type(c_ptr), intent(in), value :: Cptr2_dat
		real(pfdp) :: norm

		real(pfdp), pointer :: dat_u(:, :, :, :)

		dat_u => get_u(Cptr2_dat)
		
		norm = maxval(abs(dat_u))
	end function norm_dat
	! ------------------------------------------------------ function `norm_dat`: stop

	! ------------------------------------------------------ subroutine `pack_dat`: start
	! Packing of encapsulated data array into a flat one-dimensional array
	subroutine pack_dat(dim1_dat_u, Cptr2_dat)
		real(pfdp), intent(out) :: dim1_dat_u(:)
		type(c_ptr), intent(in), value :: Cptr2_dat

		integer :: dat_nvars
		real(pfdp), pointer :: dat_u(:, :, :, :)

		dat_nvars = get_nvars(Cptr2_dat)
		dat_u => get_u(Cptr2_dat)

		dim1_dat_u = reshape(dat_u, (/dat_nvars/))
	end subroutine pack_dat
	! ------------------------------------------------------ subroutine `pack_dat`: stop

	! ------------------------------------------------------ subroutine `unpack_dat`: start
	! Unacking of an encapsulated, flat one-dimensional data array into a non-flat data array
	subroutine unpack_dat(Cptr2_dat, dim1_dat_u)
		type(c_ptr), intent(in), value :: Cptr2_dat
		real(pfdp), intent(in) :: dim1_dat_u(:)

		integer :: dat_shape(4)
		real(pfdp), pointer :: dat_u(:, :, :, :)

		dat_shape = get_shape(Cptr2_dat)
		dat_u => get_u(Cptr2_dat)

		dat_u = reshape(dim1_dat_u, dat_shape)
	end subroutine unpack_dat
	! ------------------------------------------------------ subroutine `unpack_dat`: stop

	! ------------------------------------------------------ subroutine `axpy_dat`: start
	! Operation `a times x plus y` for a scalar `a` and data arrays `x` and `y`
	subroutine axpy_dat(Cptr2_y_dat, a, Cptr2_x_dat, flags)
		type(c_ptr), intent(in), value :: Cptr2_y_dat
		real(pfdp), intent(in) :: a
		type(c_ptr), intent(in), value :: Cptr2_x_dat
		integer, intent(in), optional :: flags

		integer :: which
		real(pfdp), pointer :: y_dat_u(:, :, :, :), x_dat_u(:, :, :, :)

		which = 0
		if(present(flags)) which = flags

		y_dat_u => get_u(Cptr2_y_dat)
		x_dat_u => get_u(Cptr2_x_dat)

		y_dat_u = a*x_dat_u+y_dat_u
	end subroutine axpy_dat
	! ------------------------------------------------------ subroutine `axpy_dat`: stop

	! ------------------------------------------------------ subroutine `print_dat`: start
	! Printing of encapsulation's data array
	subroutine print_dat(Cptr2_dat)
		type(c_ptr), intent(in), value :: Cptr2_dat

		real(pfdp), pointer :: dat_u(:, :, :, :)

		dat_u => get_u(Cptr2_dat)

		print *, 'u = ', dat_u
	end subroutine print_dat
	! ------------------------------------------------------ subroutine `print_dat`: stop

	! ------------------------------------------------------ subroutine `create_encap`: start
	! Creation of encapsulation
	subroutine create_encap(pf_encap)
		type(pf_encap_t), intent(out) :: pf_encap

		pf_encap%create => create_dat
		pf_encap%destroy => destroy_dat
		pf_encap%setval => setval_dat
		pf_encap%copy => copy_dat
		pf_encap%norm => norm_dat
		pf_encap%pack => pack_dat
		pf_encap%unpack => unpack_dat
		pf_encap%axpy => axpy_dat
		pf_encap%eprint => print_dat
	end subroutine create_encap
	! ------------------------------------------------------ subroutine `create_encap`: stop

	! ------------------------------------------------------ subroutine `destroy_encap`: start
	! Destruction of encapsulation
	subroutine destroy_encap(pf_encap)
		type(pf_encap_t), intent(inout) :: pf_encap
	end subroutine destroy_encap
	! ------------------------------------------------------ subroutine `destroy_encap`: stop
end module encap
! ------------------------------------------------------ module `encap`: stop
