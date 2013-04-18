"""WAF build script.

   To configure and build:

   $ ./waf configure
   $ ./waf build

"""

def options(opt):
  opt.load('compiler_c compiler_fc')

  opt.add_option('--fftw', action='store', default='', dest='fftw',
                 help='location of FFTW3 library, eg, --fftw=$HOME/opt')



def configure(cnf):
  cnf.load('compiler_c compiler_fc')

  if cnf.options.fftw:
    cnf.env.INCLUDES_FFTW = [ cnf.options.fftw + '/include' ]
    cnf.env.LIBPATH_FFTW  = [ cnf.options.fftw + '/lib' ]
    cnf.env.LIB_FFTW      = [ 'fftw3' ]
  else:
    cnf.check(compiler='c',
              lib='fftw3',
              mandatory=False,
              uselib_store='FFTW')

  # cnf.env.FCFLAGS   = [ '-O3' ]
  # cnf.env.CFLAGS    = [ '-O3' ]

  # f90 = '-f90=/home/memmett/gcc-4.5/bin/gfortran'
  f90 = '-f90=/home/memmett/gcc-4.8/bin/gfortran'

  cnf.env.FCFLAGS   += [ '-Wall', '-Wextra', '-g', '-pg', '-fno-strict-aliasing', '-fwrapv' ]
  cnf.env.CFLAGS    += [ '-Wall' ]
  cnf.env.LINKFLAGS += [ '-g', f90 ]

  cnf.env.FCFLAGS   += [ f90 ]
  


def build(bld):

  # libpfasst
  core = bld.path.ant_glob('src/*.c|*.f90')

  bld(features='c fc cshlib fcshlib',
      source=core,
      target='pfasst')

  # test programs
  bld(features='fc c fcprogram', use='pfasst FFTW',
      source=bld.path.ant_glob('examples/mpi-advection/*.c|*.f90'),
      target='mpi-advection')

  bld(features='fc c fcprogram', use='pfasst FFTW',
      source=bld.path.ant_glob('examples/pth-advection/*.c|*.f90'),
      target='pth-advection')

