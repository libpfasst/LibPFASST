program main
   use pfasst
   use encap

   type(scalar_encap) :: y, x
   real(c_double) :: yval, xval, a, norm
   real(pfdp), pointer :: z_pack(:)
   type(scalar_factory) :: sf
   integer :: level_index
   integer, pointer :: lev_shape(:)
   class(pf_encap_t), allocatable :: y_base, x_base

   level_index = 0
   lev_shape(1) = 0
   yval = -1.6
   xval = -.43
   a = 2.8

   call sf%create_single(y_base, level_index, lev_shape)
   y = cast_as_scalar(y_base)
   call sf%create_single(x_base, level_index, lev_shape)
   x = cast_as_scalar(x_base)

   call y%setval(yval)
   call y%eprint()
   norm = y%norm()
   print *, norm

   call x%setval(xval)
   call x%eprint()

   call y%axpy(a, x)
   call y%eprint()
   call y%copy(x)
   call y%eprint()

   call sf%destroy_single(y_base)
   call sf%destroy_single(x_base)

end program
