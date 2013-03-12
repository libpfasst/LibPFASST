!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! XXX: this module needs to be documented

module pf_mod_dnpy
  use pf_mod_dtype
  use iso_c_binding
  implicit none

  type :: pf_dnpy_t
     integer :: interval = -1
     character(len=64) :: dirname
     character(len=4)  :: endian = "<f8" // c_null_char
  end type pf_dnpy_t

  interface create
     module procedure pf_dnpy_create
  end interface create

  interface attach
     module procedure pf_dnpy_attach
  end interface attach

  interface
     subroutine pf_dnpy_solution_npy(dirname, len, endian, q, nvars, nlevel, nstep, ncycle, niter) &
          bind(c, name='pf_dnpy_solution_npy')
       use iso_c_binding
       character(c_char), intent(in)        :: dirname, endian(4)
       type(c_ptr),       intent(in), value :: q
       integer(c_int),    intent(in), value :: len, nvars, nlevel, nstep, ncycle, niter
     end subroutine pf_dnpy_solution_npy
  end interface

  interface
     subroutine pf_dnpy_npy(fname, len, endian, arr, nvars) &
          bind(c, name='pf_dnpy_npy')
       use iso_c_binding
       character(c_char), intent(in)        :: fname, endian(4)
       type(c_ptr),       intent(in), value :: arr
       integer(c_int),    intent(in), value :: len, nvars
     end subroutine pf_dnpy_npy
  end interface

  interface
     subroutine pf_dnpy_mkdir(dirname, len) &
          bind(c, name='pf_dnpy_mkdir')
       use iso_c_binding
       character(c_char), intent(in)        :: dirname
       integer(c_int),    intent(in), value :: len
     end subroutine pf_dnpy_mkdir
  end interface

contains

  subroutine pf_dnpy_create(dnpy, dirname, endian)
    type(pf_dnpy_t),  intent(out) :: dnpy
    character(len=*), intent(in)  :: dirname
    character(len=3), intent(in), optional  :: endian

    dnpy%dirname  = dirname
    dnpy%interval = 1
    if (present(endian)) then
       dnpy%endian   = endian // c_null_char
    end if

    call pf_dnpy_mkdir(dnpy%dirname, len_trim(dnpy%dirname))
  end subroutine pf_dnpy_create

  subroutine pf_dnpy_hook(pf, level, state, ctx)
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(pf_dnpy_t), pointer :: dnpy
    real(pfdp),      target  :: buf(level%nvars)

    call c_f_pointer(level%dctx, dnpy)

    call pack(buf, level%qend)

    if (mod(state%step, dnpy%interval) == 0) then
       call pf_dnpy_solution_npy(dnpy%dirname, len_trim(dnpy%dirname), dnpy%endian, &
            c_loc(buf(1)), level%nvars, level%level, state%step, state%cycle, state%iter)
    end if
  end subroutine pf_dnpy_hook

  subroutine pf_dnpy_attach(dnpy, pf, level, hook)
    use pf_mod_hooks
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level
    integer,           intent(in)    :: hook
    type(pf_dnpy_t),   intent(in), target :: dnpy

    pf%levels(level)%dctx = c_loc(dnpy)
    call add_hook(pf, level, hook, pf_dnpy_hook)
  end subroutine pf_dnpy_attach

end module pf_mod_dnpy
