"""Run several regression tests."""

import subprocess
import collections
import re

ErrorTuple = collections.namedtuple('ErrorTuple', [ 'step', 'iter', 'error' ])

def run(exe):
  p = subprocess.Popen(exe, stdout=subprocess.PIPE, shell=True)
  stdout, stderr = p.communicate()

  # print '==== stdout ===='
  # print stdout
  # print '==== stderr ===='
  # print stderr

  assert stderr is None
  return stdout


def errors(out):
  rx   = re.compile(r"error:\s*step:\s*(\d+)\s*iter:\s*(\d+)\s*error:\s*(\S+)")
  cast = [ int, int, float ]

  errors = []
  for line in out.splitlines():
    m = rx.search(line)
    if m:
      errors.append(ErrorTuple(*[ c(x) for c, x in zip(cast, m.groups()) ]))

  return errors
    

def test_mpi_advection():
  out = run('mpiexec -n 4 examples/mpi-advection/main.exe')
  err = errors(out)

  maxstep = max([ x.step for x in err ])
  maxiter = max([ x.iter for x in err ])
  lasterr = max([ x.error for x in err if x.step == maxstep and x.iter == maxiter ])

  assert lasterr < 5.0e-9


if __name__ == '__main__':
  test_mpi_advection()


