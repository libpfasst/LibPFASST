!
! Copyright (c) 2015, Michael Minion and Andreas Kreienbuehl. All rights reserved.
!

! ------------------------------------------------------ Program *main*: start
! PFASST: Main program for the Cart application
program main
	! For, e.g., C-pointers
	use iso_c_binding 

	! For LIBPFASST
	use pfasst
	! For MPI
  use pf_mod_mpi, only: MPI_COMM_WORLD

	! Problem parameters
	use probin
	! Encapsulation of data structure
	use encap
	! Function evaluations
	use feval
	! Restriction and interpolation
	use transfer
	! Diagnostic routines
	use hooks

	! Variables starting with `i`, `j`, `k`, `l`, `m`, or `n` represent integers
	implicit none

	! Input file name
  character(256) :: probin_fname

	! MPI error
	integer :: err

	! MPI rank
	integer :: rank

	! LIBPFASST data container defined in *libpfasst/src/pf_dtype_t.f90*
  type(pf_pfasst_t) :: pf

	! Communicator defined in *libpfasst/src/pf_dtype_t.f90*
  type(pf_comm_t) :: comm

	! Encapsulated data defined in *pfasst-cart/src/encap.f90* (here, ultimately, a pointer destination)
  type(pf_encap_t), target :: encap_dat

	! Timing
	real t_start, t_stop




	integer :: l





	real(pfdp), parameter :: pi_pfdp = 3.141592653589793_pfdp

	real(pfdp), parameter :: t_char_x_pfdp = 1.0_pfdp
	real(pfdp), parameter :: t_char_y_pfdp = 2.0_pfdp
	real(pfdp), parameter :: t_char_z_pfdp = 3.0_pfdp

	type(c_ptr) :: Cptr2_ini_dat
	type(dat), pointer :: Fptr2_ini_dat
	real(pfdp), pointer :: ini_dat_arr(:, :, :, :)


	type(c_ptr) :: Cptr2_end_dat

	integer :: inx, iny, inz







	!
	! Init. input
	!

  if(command_argument_count() >= 1) then
     call get_command_argument(1, value = probin_fname)
  else
     probin_fname = "probin.nml"
  end if

  call probin_init(probin_fname)

	!
	! Init. MPI
	!

  call mpi_init(err)
  if(err .ne. 0) stop "ERROR: Cannot initialize MPI"
  call mpi_comm_rank(MPI_COMM_WORLD, rank, err)

	!
	! Init. LIBPFASST
	!

	! See *libpfasst/src/pf_mpi.f90*
  call pf_mpi_create(comm, MPI_COMM_WORLD)

	! See *libpfasst/src/pf_pfasst.f90*
  call pf_pfasst_create(pf, comm, fname = probin_fname)

	!
	! Init. user data encapsulation
	!

	call create_encap(encap_dat)

	!
	! Start ticking
	!

  call cpu_time(t_start)


  !  Set some level parameters
  do l = 1, pf%nlevels
     allocate(pf%levels(l)%shape(3))
     pf%levels(l)%shape(1) = nx(l)
     pf%levels(l)%shape(2) = ny(l)
     pf%levels(l)%shape(3) = nz(l)

		 ! This is the total number of variables in flatarray
     pf%levels(l)%nvars = nflds*nx(l)*ny(l)*nz(l)

     pf%levels(l)%nnodes = nnodes(l)
     pf%levels(l)%nsweeps = nsweeps(l)

     pf%levels(l)%nsweeps_pred = nsweeps_pred(l)


     pf%levels(l)%interpolate => interpolate

     pf%levels(l)%restrict => restrict

     pf%levels(l)%encap => encap_dat

		 call pf_imexQ_create(pf%levels(l)%sweeper, f1eval, f2eval, f2comp)
  end do







	! Set initial data
	call c_f_pointer(Cptr2_ini_dat, Fptr2_ini_dat)

	call create_dat(Cptr2_ini_dat, pf%nlevels, 1, pf%levels(pf%nlevels)%nvars, pf%levels(pf%nlevels)%shape, pf%levels(pf%nlevels)%ctx)
	call create_dat(Cptr2_end_dat, pf%nlevels, 1, pf%levels(pf%nlevels)%nvars, pf%levels(pf%nlevels)%shape, pf%levels(pf%nlevels)%ctx)

	ini_dat_arr => get_arr(Cptr2_ini_dat)

	do inx = 1, nx(pf%nlevels)
		do iny = 1, ny(pf%nlevels)
			do inz = 1, nz(pf%nlevels)
				ini_dat_arr(:, inx, iny, inz) = sin(pi_pfdp*inx/t_char_x_pfdp)*sin(pi_pfdp*iny/t_char_y_pfdp)*sin(pi_pfdp*inz/t_char_z_pfdp)
			end do
		end do
	end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  print *,'Calling pfasst ...'
  if (nsteps < 0) then
     nsteps = comm%nproc
  end if
  call mpi_barrier(MPI_COMM_WORLD, err)
  call pf_pfasst_run(pf, Cptr2_ini_dat, dt, 0.0_pfdp, nsteps, Cptr2_end_dat)

  call cpu_time(t_stop)

  call mpi_barrier(MPI_COMM_WORLD, err)

!	call setval_dat(Cptr2_ini_dat, ini_dat_arr, 0)
end program main
! ------------------------------------------------------ Program *main*: stop
