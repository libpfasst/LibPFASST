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

  cnf.env.FCFLAGS   += [ '-Wall', '-g', '-pg' ]
  cnf.env.CFLAGS    += [ '-Wall' ]
  cnf.env.LINKFLAGS += [ '-g' ]


def build(bld):

  # libpfasst
  core = bld.path.ant_glob('src/*.c|*.f90')

  bld(features='c fc cshlib fcshlib',
      source=core,
      target='pfasst')

  # test programs
  bld(features='fc c fcprogram', use='pfasst FFTW',
      source=bld.path.ant_glob('examples/advection/*.c|*.f90'),
      target='advection')

