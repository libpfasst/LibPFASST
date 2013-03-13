"""WAF build script.

   To configure and build:

   $ ./waf configure
   $ ./waf build

"""

def options(opt):
  opt.load('compiler_c compiler_fc')

def configure(cnf):
  cnf.load('compiler_c compiler_fc')

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
  # ...
