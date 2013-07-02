!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none

  interface
     subroutine dump_mkdir(dname, dlen) bind(c)
       use iso_c_binding
       character(c_char), intent(in) :: dname
       integer(c_int),    intent(in), value :: dlen
     end subroutine dump_mkdir

     subroutine dump_numpy(dname, fname, endian, dim, shape, nvars, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in) :: dname, fname, endian(4)
       integer(c_int),    intent(in), value :: dim, nvars
       integer(c_int),    intent(in) :: shape(dim)
       real(c_double),    intent(in) :: array(nvars)
     end subroutine dump_numpy
  end interface

contains

  subroutine echo_error_hook(pf, level, state, ctx)
    use solutions, only: exact
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: ctx

    real(c_double) :: yexact(level%nvars)
    real(pfdp), pointer :: qend(:)

    qend => array1(level%qend)

    call exact(state%t0+state%dt, level%nvars, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(qend-yexact))

  end subroutine echo_error_hook


  subroutine dump_hook(pf, level, state, ctx)
    use probin, only: output
    type(pf_pfasst_t),   intent(inout) :: pf
    type(pf_level_t),    intent(inout) :: level
    type(pf_state_t),    intent(in)    :: state
    type(c_ptr),         intent(in)    :: ctx

    character(len=256)  :: fname
    real(pfdp), pointer :: qend(:)

    qend => array1(level%qend)

    write(fname, "('s',i0.5,'i',i0.3,'l',i0.2,'.npy')") &
         state%step, state%iter, level%level

    call dump_numpy(trim(output)//c_null_char, trim(fname)//c_null_char, '<f8'//c_null_char, &
         1, [ size(qend) ], size(qend), qend)

  end subroutine dump_hook

end module hooks
