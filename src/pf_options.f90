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

  subroutine pf_read_opts(pf, read_cmd, fname)
    type(pf_pfasst_t), intent(inout)           :: pf
    logical,           intent(in   )           :: read_cmd
    character(len=*),  intent(in   ), optional :: fname

    ! local versions of pfasst parameters
    integer          :: niters, nlevels, qtype, window, taui0
    double precision :: abs_res_tol, rel_res_tol
    logical          :: pipeline_g , pfasst_pred, echo_timings

    ! stuff for reading the command line
    integer, parameter :: un = 9
    integer            :: i, ios
    character(len=32)  :: arg
    character(len=255) :: istring  ! stores command line argument
    character(len=255) :: message  ! use for i/o error messages
    character(len=512) :: outdir

    ! define the namelist for reading
    namelist /pf_params/ niters, nlevels, qtype, abs_res_tol, rel_res_tol, window
    namelist /pf_params/ pipeline_g, pfasst_pred, echo_timings, taui0, outdir

    ! set local variables to pf_pfasst defaults
    nlevels      = pf%nlevels
    niters       = pf%niters
    qtype        = pf%qtype
    window       = pf%window
    abs_res_tol  = pf%abs_res_tol
    rel_res_tol  = pf%rel_res_tol
    pipeline_g   = pf%pipeline_g
    pfasst_pred  = pf%pfasst_pred
    echo_timings = pf%echo_timings
    taui0        = pf%taui0
    outdir       = pf%outdir

    ! open the file fname and read the pfasst namelist
    if (present(fname))  then
       open(unit=un, file=fname, status='old', action='read')
       read(unit=un, nml=pf_params)
       close(unit=un)
    end if

    ! overwrite with the command line
    if (read_cmd) then
       i = 0
       do
          call get_command_argument(i, arg)
          if (len_trim(arg) == 0) exit
          if (i > 0) then
             istring="&pf_params " // trim(arg) // " /"
             read(istring, nml=pf_params, iostat=ios, iomsg=message) ! internal read of namelist
          end if
          i = i+1
       end do
    end if

    ! re-assign the pfasst internals
    pf%niters       = niters
    pf%nlevels      = nlevels
    pf%qtype        = qtype
    pf%window       = window
    pf%abs_res_tol  = abs_res_tol
    pf%rel_res_tol  = rel_res_tol
    pf%pipeline_g   = pipeline_g
    pf%pfasst_pred  = pfasst_pred
    pf%echo_timings = echo_timings
    pf%taui0        = taui0
    pf%outdir       = outdir

    if (pf%nlevels < 1) then
       write(*,*) 'Bad specification for nlevels=', pf%nlevels
       stop
    endif
  end subroutine pf_read_opts

  subroutine pf_print_options(pf, unitno, show_mats)
    type(pf_pfasst_t), intent(inout)           :: pf
    integer,           intent(in   ), optional :: unitno
    logical,           intent(in   ), optional :: show_mats

    integer :: un = 6
    integer :: l, i
    character(8)   :: date
    character(10)  :: time

        print *,'Print options', pf%rank
    if (pf%rank /= 0) return
    if (present(unitno)) un = unitno

    write(un,*) 'PFASST Configuration'
    write(un,*) '===================='

    call date_and_time(date=date, time=time)
    write(un,*) 'date:        ', date
    write(un,*) 'time:        ', time

    write(un,*) 'nlevels:     ', pf%nlevels, '! number of pfasst levels'
    write(un,*) 'nprocs:      ', pf%comm%nproc, '! number of pfasst "time" processors'
    if (pf%comm%nproc == 1) then
       write(un,*) '            ', '             ', ' ! since 1 time proc is being used, this is a serial sdc run'
    else
       write(un,*) '            ', '             ', ' ! since >1 time procs are being used, this is a parallel pfasst run'
    end if
    write(un,*) 'niters:      ', pf%niters, '! maximum number of sdc/pfasst iterations'
    write(un,*) 'nnodes:      ', pf%levels(1:pf%nlevels)%nnodes, '! number of sdc nodes per level'
    write(un,*) 'nvars:       ', pf%levels(1:pf%nlevels)%nvars, '! number of degrees of freedom per level'
    write(un,*) 'nsweeps:     ', pf%levels(1:pf%nlevels)%nsweeps, '! number of sdc sweeps performed per visit to each level'
    write(un,*) 'nsweeps_pred:     ', pf%levels(1:pf%nlevels)%nsweeps_pred, '! number of sdc sweeps in predictor'
    write(un,*) 'taui0:     ',   pf%taui0, '! cutoff for tau correction'

    if (pf%comm%nproc > 1) then
       if (pf%window == PF_WINDOW_RING) then
          write(un,*) 'window:     ', '      "ring"', ' ! pfasst processors advance through time in a ring'
       else
          write(un,*) 'window:     ', '      "block"', ' ! pfasst processors advance through time as a block'
       end if
    end if

    if (pf%Pipeline_G) then
       write(un,*) 'Predictor Pipelining is ON    '
    else
       write(un,*) 'Predictor Pipelining is OFF    '
    end if
    if (pf%PFASST_pred) then
       write(un,*) 'PFASST Predictor style  '
    else
       write(un,*) 'Serial Predictor style  '
    end if

    write(un,*) ''

    if (present(show_mats)) then
       if (show_mats) then
          do l = 1, pf%nlevels
             print *, "Level", l
             print *, "-----------------"
             print *, "  nodes"
             print *, pf%levels(l)%nodes
             ! print *, "  Q"
             ! do i = 1, pf%levels(l)%nnodes-1
             !    print *, pf%levels(l)%qmat(i,:)
             ! end do
             print *, "  S"
             do i = 1, pf%levels(l)%nnodes-1
                print *, pf%levels(l)%s0mat(i,:)
             end do
          end do
       end if
    end if


  end subroutine pf_print_options

end module pf_mod_options
