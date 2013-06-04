module output
  use iso_c_binding

  interface
     subroutine dump_mkdir(dname, dlen) bind(c)
       use iso_c_binding
       character(c_char), intent(in) :: dname
       integer(c_int),    intent(in), value :: dlen
     end subroutine dump_mkdir

     subroutine dump_numpy(fname, flen, endian, dim, shape, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in) :: fname, endian(4)
       integer(c_int),    intent(in), value :: flen, dim
       integer(c_int),    intent(in) :: shape(dim)
       real(c_double),    intent(in) :: array
     end subroutine dump_numpy
  end interface

end module output
