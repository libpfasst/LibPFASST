!
! Copyright (c) 2015, Michael Minion and Andreas Kreienbuehl. All rights reserved.
!

! ------------------------------------------------------ module *encap*: start
! PFASST: Encapsulation for compressible Navier-Stokes equation (Cart)
module encap
	! For, e.g., C-pointers
	use iso_c_binding 

	! For LIBPFASST
	use pfasst

	! Variables starting with *i*, *j*, *k*, *l*, *m*, or *n* represent integers
	implicit none

	! We declare the encapsulation's data (i.e. `dat`)
	type :: dat
		! We have five field variables (i.e. `nflds = 5`): `rho` (for density) , `P = rho*(u, v, w)` (for momentum density), and `E` (for energy)
		integer :: nflds
		integer :: nvars

		! For *convenience* we include an F-array carrying the three integers `nx`, `ny`, `nz`
		integer :: shape(3)

		! We live on a three-dimensional grid
		integer :: nx, ny, nz

		! The encapsulated data are `nflds` fields on a three-dimensional grid or `arr = [rho, rho*(u, v, w), E](:, :, :)`
		real(pfdp), pointer :: arr(:, :, :, :)
	end type dat

contains
	! ------------------------------------------------------ function *get_nflds*: start
	! We get an encapsulation's number of field variables
	function get_nflds(Cptr2_dat) result(dat_nflds)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_nflds

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_nflds = Fptr2_dat%nflds
	end function get_nflds
	! ------------------------------------------------------ function *get_nflds*: stop

	! ------------------------------------------------------ function *get_shape*: start
	! We get an encapsulation's data array's size, i.e. `nx`, `ny`, and `nz` value
	function get_shape(Cptr2_dat) result(dat_shape)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_shape(3)

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_shape = Fptr2_dat%shape
	end function get_shape
	! ------------------------------------------------------ function *get_shape*: stop

	! ------------------------------------------------------ function *get_nx*: start
	! We get an encapsulation's number of elements in `x`-direction
	function get_nx(Cptr2_dat) result(dat_nx)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_nx

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_nx = Fptr2_dat%nx
	end function get_nx
	! ------------------------------------------------------ function *get_nx*: stop

	! ------------------------------------------------------ function *get_ny*: start
	! We get an encapsulation's number of elements in `y`-direction
	function get_ny(Cptr2_dat) result(dat_ny)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_ny

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_ny = Fptr2_dat%ny
	end function get_ny
	! ------------------------------------------------------ function *get_ny*: stop

	! ------------------------------------------------------ function *get_nz*: start
	! We get an encapsulation's number of elements in `z`-direction
	function get_nz(Cptr2_dat) result(dat_nz)
		type(c_ptr), intent(in), value :: Cptr2_dat
		integer :: dat_nz

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_nz = Fptr2_dat%nz
	end function get_nz
	! ------------------------------------------------------ function *get_nz*: stop

	! ------------------------------------------------------ function *get_arr*: start
	! We get an encapsulation's data array
	function get_arr(Cptr2_dat) result(dat_arr)
		type(c_ptr), intent(in), value :: Cptr2_dat
		real(pfdp), pointer :: dat_arr(:, :, :, :)

		type(dat), pointer :: Fptr2_dat

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		dat_arr => Fptr2_dat%arr
	end function get_arr
	! ------------------------------------------------------ function *get_arr*: stop

	! ------------------------------------------------------ subroutine *create_dat*: start
	! We create the encapsulation's data
	subroutine create_dat(Cptr2_dat, lvl, kind, nvars, shape, Cptr2_ctx)
		type(c_ptr), intent(inout) :: Cptr2_dat
		integer, intent(in) :: lvl, kind, nvars, shape(:)
		type(c_ptr), intent(in), value :: Cptr2_ctx

		type(dat), pointer :: Fptr2_dat
		integer :: ierr

		allocate(Fptr2_dat, STAT = ierr)
		if(ierr .ne. 0) stop 'Allocation of `Fptr2_dat` in `def_dat` failed'

		Fptr2_dat%shape = shape
Fptr2_dat%nflds = 5

		Fptr2_dat%nx = Fptr2_dat%shape(1)
		Fptr2_dat%ny = Fptr2_dat%shape(2)
		Fptr2_dat%nz = Fptr2_dat%shape(3)

		Fptr2_dat%nvars = Fptr2_dat%nflds*Fptr2_dat%nx*Fptr2_dat%ny*Fptr2_dat%nz

		allocate(Fptr2_dat%arr(Fptr2_dat%nvars, Fptr2_dat%nx, Fptr2_dat%ny, Fptr2_dat%nz), STAT = ierr)
		if(ierr .ne. 0) stop 'Allocation of `Fptr2_dat%arr` in `def_dat` failed'

		Cptr2_dat = c_loc(Fptr2_dat)
	end subroutine create_dat
	! ------------------------------------------------------ subroutine *create_dat*: stop

	! ------------------------------------------------------ subroutine *destroy_dat*: start
	! We destroy the encapsulation's data
	subroutine destroy_dat(Cptr2_dat)
		type(c_ptr), intent(in), value :: Cptr2_dat

		type(dat), pointer :: Fptr2_dat
		integer :: ierr

		call c_f_pointer(Cptr2_dat, Fptr2_dat)

		deallocate(Fptr2_dat%arr, STAT = ierr)
		if(ierr .ne. 0) stop 'Deallocation of `Fptr2_dat%arr` in `destroy_dat` failed'

		deallocate(Fptr2_dat, STAT = ierr)
		if(ierr .ne. 0) stop 'Deallocation of `Fptr2_dat` in `destroy_dat` failed'
	end subroutine destroy_dat
	! ------------------------------------------------------ subroutine *destroy_dat*: stop

	! ------------------------------------------------------ subroutine *setval_dat*: start
	! We set the encapsulation's data values
	subroutine setval_dat(Cptr2_dst_dat, src_dat_arr, flags)
		type(c_ptr), intent(in), value :: Cptr2_dst_dat
		real(pfdp), intent(in) :: src_dat_arr
		integer, intent(in), optional :: flags

		real(pfdp), pointer :: dst_dat_arr(:, :, :, :)
		integer :: which

		which = 0
		if(present(flags)) which = flags

		dst_dat_arr => get_arr(Cptr2_dst_dat)

		dst_dat_arr = src_dat_arr
	end subroutine setval_dat
	! ------------------------------------------------------ subroutine *setval_dat*: stop

	! ------------------------------------------------------ subroutine *copy_dat*: start
	! We copy the encapsulation's data array (from source `src` to destination `dst`)
	subroutine copy_dat(Cptr2_dst_dat, Cptr2_src_dat, flags)
		type(c_ptr), intent(in), value :: Cptr2_dst_dat, Cptr2_src_dat
		integer, intent(in), optional :: flags

		real(pfdp), pointer :: dst_dat_arr(:, :, :, :), src_dat_arr(:, :, :, :)
		integer :: which

		which = 0
		if (present(flags)) which = flags

		dst_dat_arr => get_arr(Cptr2_dst_dat)
		src_dat_arr => get_arr(Cptr2_src_dat)

		dst_dat_arr = src_dat_arr
	end subroutine copy_dat
	! ------------------------------------------------------ subroutine *copy_dat*: stop

	! ------------------------------------------------------ function *norm_dat*: start
	! We determine the norm of the encapsulation's data array
	function norm_dat(Cptr2_dat) result(norm)
		type(c_ptr), intent(in), value :: Cptr2_dat
		real(pfdp) :: norm

		real(pfdp), pointer :: dat_arr(:, :, :, :)

		dat_arr => get_arr(Cptr2_dat)
		
		norm = maxval(abs(dat_arr))
	end function norm_dat
	! ------------------------------------------------------ function *norm_dat*: stop

	! ------------------------------------------------------ subroutine *pack_dat*: start
	! We pack the encapsulation's data array into a flat array
	subroutine pack_dat(dim1_dat_arr, Cptr2_dat)
		real(pfdp), intent(out) :: dim1_dat_arr(:)
		type(c_ptr), intent(in), value :: Cptr2_dat

		integer :: dat_nflds
		integer :: dat_nx, dat_ny, dat_nz
		real(pfdp), pointer :: dat_arr(:, :, :, :)

		integer :: dims(1)

		dat_nflds = get_nflds(Cptr2_dat)
		dat_nx = get_nx(Cptr2_dat)
		dat_ny = get_ny(Cptr2_dat)
		dat_nz = get_nz(Cptr2_dat)
		dat_arr => get_arr(Cptr2_dat)

dims = (/dat_nflds*dat_nx*dat_ny*dat_nz/)

print *, "HELLO YOU 1"

		dim1_dat_arr = reshape(dat_arr, dims)
	end subroutine pack_dat
	! ------------------------------------------------------ subroutine *pack_dat*: stop

	! ------------------------------------------------------ subroutine *unpack_dat*: start
	! We unpack a flat array into the encapsulation's data array
	subroutine unpack_dat(Cptr2_dat, dim1_dat_arr)
		type(c_ptr), intent(in), value :: Cptr2_dat
		real(pfdp), intent(in) :: dim1_dat_arr(:)

		integer :: dat_nflds
		integer :: dat_nx, dat_ny, dat_nz
		real(pfdp), pointer :: dat_arr(:, :, :, :)


		integer :: dims(4)

		dat_nflds = get_nflds(Cptr2_dat)
		dat_nx = get_nx(Cptr2_dat)
		dat_ny = get_ny(Cptr2_dat)
		dat_nz = get_nz(Cptr2_dat)
		dat_arr => get_arr(Cptr2_dat)

		dims = (/dat_nflds, dat_nx, dat_ny, dat_nz/)

print *, "HELLO YOU 2", dims, dat_nflds

		dat_arr = reshape(dim1_dat_arr, dims)
	end subroutine unpack_dat
	! ------------------------------------------------------ subroutine *unpack_dat*: stop

	! ------------------------------------------------------ subroutine *axpy_dat*: start
	! We define the operation `a times x plus y` for the encapsulation's data array
	subroutine axpy_dat(Cptr2_y_dat, a, Cptr2_x_dat, flags)
		type(c_ptr), intent(in), value :: Cptr2_y_dat
		real(pfdp), intent(in) :: a
		type(c_ptr), intent(in), value :: Cptr2_x_dat
		integer, intent(in), optional :: flags

		integer :: which
		real(pfdp), pointer :: y_dat_arr(:, :, :, :), x_dat_arr(:, :, :, :)

		which = 0
		if(present(flags)) which = flags

		y_dat_arr => get_arr(Cptr2_y_dat)
		x_dat_arr => get_arr(Cptr2_x_dat)

		y_dat_arr = a * x_dat_arr + y_dat_arr
	end subroutine axpy_dat
	! ------------------------------------------------------ subroutine *axpy_dat*: stop

	! ------------------------------------------------------ subroutine *print_dat*: start
	! We print the encapsulation's data
	subroutine print_dat(Cptr2_dat)
		type(c_ptr), intent(in), value :: Cptr2_dat

		real(pfdp), pointer :: dat_arr(:, :, :, :)

		dat_arr => get_arr(Cptr2_dat)

		print *, 'arr = ', dat_arr
	end subroutine print_dat
	! ------------------------------------------------------ subroutine *print_dat*: stop

	! ------------------------------------------------------ subroutine *create_encap*: start
	! We create the capsulation's data and context
	subroutine create_encap(encap_dat)
		type(pf_encap_t), intent(out) :: encap_dat

		encap_dat%create => create_dat
		encap_dat%destroy => destroy_dat
		encap_dat%setval => setval_dat
		encap_dat%copy => copy_dat
		encap_dat%norm => norm_dat
		encap_dat%pack => pack_dat
		encap_dat%unpack => unpack_dat
		encap_dat%axpy => axpy_dat
		encap_dat%eprint => print_dat
	end subroutine create_encap
	! ------------------------------------------------------ subroutine *create_encap*: stop

	! ------------------------------------------------------ subroutine *destroy_encap*: start
	! We destroy the encapsulation's data and context
	subroutine destroy_encap(encap_dat)
		type(pf_encap_t), intent(inout) :: encap_dat
	end subroutine destroy_encap
	! ------------------------------------------------------ subroutine *destroy_encap*: stop
end module encap
! ------------------------------------------------------ module *encap*: stop
