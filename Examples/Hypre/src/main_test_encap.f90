program main
   use pfasst
   use encap
   use pf_my_sweeper
   use pf_mod_mpi
   use pf_my_level
   use probin
   use pf_space_comm
   use pfasst_hypre

   type(pf_pfasst_t) :: pf
   type(hypre_vector_encap) :: x, y, z, v_f, v_c
   real(c_double) :: yval, xval, zval, vval, a, norm
   real(pfdp), pointer :: z_pack(:), v(:)
   type(hypre_vector_factory) :: sf
   integer :: level_index, l
   integer, allocatable :: lev_shape(:,:)
   class(pf_encap_t), allocatable :: y_base, x_base, z_base, v_f_base, v_c_base
   integer :: n
   integer :: ierror
   type(my_sweeper_t) :: s, s_finest
   type(my_level_t) :: my_lev
   type(pf_comm_t) :: comm 
   integer :: space_comm, time_comm, space_color, time_color
   character(256) :: pf_fname
   integer :: piece
   real(pfdp) :: t, dtq
   integer :: f_level, c_level

   integer :: nproc, rank, error

   call mpi_init(ierror)
   call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
   call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

   call probin_init(pf_fname)

   call create_simple_communicators(nspace, ntime, space_comm, time_comm, space_color, time_color, space_dim)

   call pf_mpi_create(comm, time_comm)

   call pf_pfasst_create(pf, comm, fname=pf_fname)

   call PfasstHypreInit(pf, lev_shape, space_color, time_color)

   n = num_grid_points**2 
   level_index = 1
   yval = -1.6
   xval = -.43
   a = 2.8
   piece = 2
   t = .1
   dtq = .01


   print *,'Initialize vectors:'
   call sf%create_single(x_base, level_index, lev_shape(level_index,:))
   x = cast_as_hypre_vector(x_base)
   call sf%create_single(y_base, level_index, lev_shape(level_index,:))
   y = cast_as_hypre_vector(y_base)
   call sf%create_single(z_base, level_index, lev_shape(level_index,:))
   z = cast_as_hypre_vector(z_base)

   print *,'Set sin init cond:'
   call initial(y)
   call y%eprint()

   print *,'Setval:'
   call y%setval(yval)
   call y%eprint()

   call x%setval(xval)
   call x%eprint()

   norm = y%norm()
   print *,'Norm:'
   print *,norm
   print *,''

   print *,'Copy:'
   call y%copy(x)
   call y%eprint()

   print *,'Axpy:'
   call y%axpy(a, x)
   call y%eprint()

   print *,'Pack then unpack:'
   allocate(v(n))
   call x%setval(a)
   call x%pack(v)
   call y%unpack(v)
   call y%eprint()

   xval = 1.0
   call x%setval(xval)
   zval = 2.0
   call z%setval(zval)
   s = cast_as_my_sweeper_t(pf%levels(level_index)%ulevel%sweeper)
   print *,'FEval:'
   call s%f_eval(x, t, level_index, y, piece)
   call y%eprint()

   xval = 1.0
   call x%setval(xval)
   zval = 2.0
   call z%setval(zval)
   print *,'FComp:'
   call s%f_comp(y, t, dtq, x, level_index, z, piece)
   call y%eprint()   

   if (pf%nlevels > 1) then
      f_level = level_index + 1
      c_level = level_index
      call sf%create_single(v_f_base, f_level, lev_shape(f_level,:))
      v_f = cast_as_hypre_vector(v_f_base)
      call sf%create_single(v_c_base, c_level, lev_shape(c_level,:))
      v_c = cast_as_hypre_vector(v_c_base)

      my_lev = cast_as_my_level_t(pf%levels(level_index)%ulevel);
      vval = 1.0
      call v_f%setval(vval)
      vval = 0.0
      call v_c%setval(vval)
      call my_lev%restrict(pf%levels(f_level), pf%levels(c_level), v_f, v_c, t);
      print *,'Restrict:'
      call v_c%eprint()

      my_lev = cast_as_my_level_t(pf%levels(level_index)%ulevel);
      vval = 1.0
      call v_c%setval(vval)
      vval = 0.0
      call v_f%setval(vval)
      call my_lev%interpolate(pf%levels(f_level), pf%levels(c_level), v_f, v_c, t);
      print *,'Interpolate:'
      call v_f%eprint()
   end if

   !if (rank == 0) then
   !end if
   !call mpi_barrier(MPI_COMM_WORLD, error)
   !if (rank == 1) then
   !   call y%eprint()
   !end if
   
   call s%destroy(pf, level_index)
   call sf%destroy_single(y_base)
   call sf%destroy_single(x_base)

   call mpi_finalize(ierror)

end program
