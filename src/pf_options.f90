!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

module pf_mod_options
  use pf_mod_dtype
  implicit none
contains

  subroutine pf_opts_from_cl(pf)
    type(pf_pfasst_t), intent(inout) :: pf
    character(len=128) :: arg
    integer            :: argc, argi
    
    argc = command_argument_count()
    argi = 1
    do while (argi <= argc)
       call get_command_argument(argi, arg)
       select case(arg)
       case ("--pf-ring")
          pf%window = PF_WINDOW_RING
       case ("--pf-niters")
          argi = argi + 1
          call get_command_argument(argi, arg)
          read(arg,*) pf%niters
       case ("--pf-abs-res-tol")
          argi = argi + 1
          call get_command_argument(argi, arg)
          read(arg,*) pf%abs_res_tol
       case ("--pf-rel-res-tol")
          argi = argi + 1
          call get_command_argument(argi, arg)
          read(arg,*) pf%rel_res_tol
       case ('--')
          exit
       case default
       end select
       argi = argi + 1
    end do
  end subroutine pf_opts_from_cl



end module pf_mod_options

