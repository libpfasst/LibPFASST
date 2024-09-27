"""Run several regression tests."""

import subprocess   #  python native module to run processes from python
import collections  #  python native module of container classes
import re           #  python native module to compare regular expressions
import pytest


ErrorTuple = collections.namedtuple('ErrorTuple', [ 'rank', 'step', 'iter', 'level', 'residual' ])   # This is what is scraped from the output
EXE = 'test/nagumo/main_split.exe'                                                        # The name of the executable to run
NMLFILE = 'test/nagumo/{}.nml'                                                            # The name of the input files
TOL = 1e-11                                                                               # The error tolerance that must be met for a successful test

"""make pfasst tests from the mlsdc.nml input file."""
def make_pfasst():
    pfasst = []
    nml = NMLFILE.format('probin')
    mpi_tasks = 1 #32
    for max_opt_iter in [0, 1]:   # 0: check state solution; 1: check adjoint solution
       pfasst.append((mpi_tasks, nml, max_opt_iter))

    return pfasst

"""scrape the output looking for the error statement"""
def errors(out):
    #rx = re.compile(r"error:\s*step:\s*(\d+)\s*iter:\s*(\d+)\s*level:\s*(\d+)\s*error:\s*(\S+)")
    rx = re.compile(r"rank:\s*(\d+)\s*lev:\s*(\d+)\s*step:\s*(\d+)\s*iter:\s*(\d+)\s*res:\s*(\S+)")
# rank:    0 lev:    3 step:    1 t= 1.562E-01 t0= 0.000E+00 iter: 14 Max_y: 1.866421E+00 Max_p: 5.796460E+00 yres: 7.590695E-11 pres: 3.757933E-11 Res: 2.202122E-12

    cast = [int, int, int, int, float]

    errors = []
    for line in out.splitlines():
        m = rx.search(line)
        if m:
            try:
                errors.append(ErrorTuple(*[c(x) for c, x in zip(cast, m.groups())]))
            except ValueError:
                raise ValueError

    return errors

"""set up the list of tests to do and call the routines to make the list of tests"""
tests = []
tests.extend(make_pfasst())

"""pytest command to call the following routing with the different parameters loaded into tests """
@pytest.mark.parametrize('mpi_tasks, nml, max_opt_iter',
                          tests)
# later: add misdc sweeper as well

def test_nagumo(mpi_tasks, nml, max_opt_iter):
    command = 'mpirun -np {} {} {} max_opt_iter={}'.format(mpi_tasks, EXE, nml, max_opt_iter)

    print(command)
    output = subprocess.check_output(command.split(), text=True)   #  This will run all the tests

    try:
        err = errors(output)   # Try to read the output for error statistics
    except ValueError:
        assert False, \
                'Unexpected poor convergence behavior\nunable to parse output\n {}'.format(line)
    else:   # Find the error from the last step and last iteratoin 
        if max_opt_iter == 0:
          maxrank = max([x.rank for x in err])
          maxstep = max([x.step for x in err if x.rank == maxrank])
          maxiter = max([x.iter for x in err if x.step == maxstep])
          lasterr = max([x.residual for x in err if x.step == maxstep and x.iter == maxiter])
          
          assert lasterr < TOL, "error: {}, tol: {}".format(lasterr, TOL)   # This decides if the test was successful

        if max_opt_iter == 1: # for the adjoint: rank 0 is the last one, with supposedly the largest error
          minrank = min([x.rank for x in err])
          maxstep = max([x.step for x in err if x.rank == minrank])
          maxiter = max([x.iter for x in err if x.step == maxstep and x.rank == minrank])
          lasterr = max([x.residual for x in err if x.step == maxstep and x.iter == maxiter and x.rank == minrank])
          #print(minrank, maxstep, maxiter, lasterr)
          
          assert lasterr < TOL, "error: {}, tol: {}".format(lasterr, TOL)   # This decides if the test was successful


