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

  subroutine pf_opts_from_file(pf, un)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: un

    

  end subroutine pf_opts_from_file

  subroutine pf_set_options(pf, fname, unitno, cmdline)
    type(pf_pfasst_t), intent(inout) :: pf
    character(len=*),  intent(in), optional :: fname
    integer,           intent(in), optional :: unitno
    logical,           intent(in), optional :: cmdline

    if (present(fname) .and. len_trim(fname) > 0) then
       open(unit=66, file=fname, status='old', action='read')
       call pf_opts_from_file(pf, 66)
       close(unit=66)
    end if
    if (present(unitno)) then
       call pf_opts_from_file(pf, unitno)
    end if
    if (.not. present(cmdline) .or. cmdline) then
       call pf_opts_from_cl(pf)
    end if
  end subroutine pf_set_options

  subroutine pf_print_options(pf, unitno)
    use pf_mod_version

    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in), optional :: unitno

    integer :: un = 6
    character(8)   :: date
    character(10)  :: time

    if (pf%rank /= 0) return
    if (present(unitno)) un = unitno

    write(un,*) 'PFASST Configuration'
    write(un,*) '===================='

    call date_and_time(date=date, time=time)
    write(un,*) 'date:        ', date
    write(un,*) 'time:        ', time
    write(un,*) 'version:     ', pf_version
    write(un,*) 'git version: ', pf_git_version

    write(un,*) 'nlevels:     ', pf%nlevels, '! number of pfasst levels'
    write(un,*) 'nprocs:      ', pf%comm%nproc, '! number of pfasst "time" processors'
    write(un,*) 'niters:      ', pf%niters, '! maximum number of sdc/pfasst iterations'
    write(un,*) 'nnodes:      ', pf%levels(:)%nnodes, '! number of sdc nodes per level'
    write(un,*) 'nvars:       ', pf%levels(:)%nvars, '! number of degrees of freedom per level'
    write(un,*) 'nsweeps:     ', pf%levels(:)%nsweeps, '! number of sdc sweeps performed per visit to each level'
    if (pf%window == PF_WINDOW_RING) then
       write(un,*) 'window:     ', '      "ring"', ' ! pfasst processors advance through time in a ring'
    else
       write(un,*) 'window:     ', '      "block"', ' ! pfasst processors advance through time as a block'
    end if

  end subroutine pf_print_options

end module pf_mod_options

