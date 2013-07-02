"""Run several regression tests."""

import subprocess

def run(exe, tol=1.0e-14):
  p = subprocess.Popen('''%s | grep "error:" | egrep -v "iter: +-" | sort -n ''' % exe,
                       stdout=subprocess.PIPE,
                       shell=True)

  (stdout, stderr) = p.communicate()

  print '==== stdout ===='
  print stdout

  print '==== stderr ===='
  print stderr

  # from the output, grab the last 'word' of the last line
  last = stdout.splitlines()[-1]
  err  = float(last.split()[-1])

  print '==== checks ===='
  print 'error of last step and iteration:', err

  assert(stderr is None)
  assert(err < tol)


# def test_mpi_advection():
#   run('mpirun -n 4 bin/pfasst-mpi-advection.exe')

# def test_fake_advection():
#   run('bin/pfasst-fake-advection.exe', tol=5.0e-14)

# def test_pth_advection():
#   run('bin/pfasst-pth-advection.exe')

# def test_rk_advection():
#   run('bin/pfasst-rk-advection.exe', tol=5.0e-8)

# def test_mpi_harmV():
#   run('mpirun -n 8 bin/pfasst-mpi-harmV.exe')

def test_t1():
  run('mpiexec -n 4 tests/t1.exe', tol=5e-9)



if __name__ == '__main__':
  test_t1()
  #test_mpi_advection()
  #test_pth_advection()
  #test_rk_advection()
  # test_mpi_harmV()
  #test_fake_advection()

