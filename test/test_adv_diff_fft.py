"""Run several regression tests."""

import subprocess   #  python native module to run processes from python
import collections  #  python native module of container classes
import re           #  python native module to compare regular expressions
import pytest


ErrorTuple = collections.namedtuple('ErrorTuple', [ 'step', 'iter', 'level', 'error' ])   # This is what is scraped from the output
EXE = 'test/adv_diff_fft/1d/main.exe'                                                        # The name of the executable to run
NMLFILE = 'test/adv_diff_fft/{}.nml'                                                      # The name of the input files
TOL = 1e-8                                                                               # The error tolerance that must be met for a successful test

"""make tests from the sdc.nml input file."""
def make_sdc():
    sdc = []
    nml = NMLFILE.format('sdc')
    mpi_tasks = 1
    nsteps = 32
    for nodes in range(3, 9):
        for imex_status in [0, 2]:
            sdc.append((nodes, imex_status, mpi_tasks, nml, nsteps))

    return sdc

"""make serial tests from the mlsdc.nml input file."""
def make_mlsdc():
    mlsdc = []
    nml = NMLFILE.format('multi_level')
    mpi_tasks = 1
    nodes = '[2 3 5]'
    for imex_status in [0, 2]:
        if imex_status == 0:
            nsteps = 256
        else:
            nsteps = 32

        mlsdc.append((nodes, imex_status, mpi_tasks, nml,nsteps))

    return mlsdc

"""make pfasst tests from the mlsdc.nml input file."""
def make_pfasst():
    pfasst = []
    nml = NMLFILE.format('multi_level')
    mpi_tasks = 4
    nodes = '[2 3 5]'
    for imex_status in [0, 2]:
        if imex_status == 0:
            nsteps = 256
        else:
            nsteps = 32

        pfasst.append((nodes, imex_status, mpi_tasks, nml,nsteps))

    return pfasst

"""scrape the output looking for the error statement"""
def errors(out):
    rx = re.compile(r"error:\s*step:\s*(\d+)\s*iter:\s*(\d+)\s*level:\s*(\d+)\s*error:\s*(\S+)")
    cast = [int, int, int, float]

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
tests.extend(make_sdc())
tests.extend(make_mlsdc())
tests.extend(make_pfasst())

"""pytest command to call the following routing with the different parameters loaded into tests """
@pytest.mark.parametrize('nodes, imex_status, mpi_tasks, nml, nsteps',
                          tests)
def test_advdiff(nodes, imex_status, mpi_tasks, nml,nsteps):
    command = 'mpirun -np {} {} {} nnodes={} imex_stat={} nsteps={}'.format(mpi_tasks, EXE, nml, nodes, imex_status,nsteps)

    print command
    output = subprocess.check_output(command.split())   #  This will run all the tests

    try:
        err = errors(output)   # Try to read the output for error statistics
    except ValueError:
        if imex_status == 0:
            assert True, \
                'Expected poor convergence behavior\nunable to parse output\n {}'.format(line)
        else:
            assert False, \
                'Unexpected poor convergence behavior\nunable to parse output\n {}'.format(line)
    else:   # Find the error from the last step and last iteratoin 
        maxstep = max([x.step for x in err])
        maxiter = max([x.iter for x in err if x.step == maxstep])
        lasterr = max([x.error for x in err if x.step == maxstep and x.iter == maxiter])

        assert lasterr < TOL, "error: {}, tol: {}".format(lasterr, TOL)   # This decides if the test was successful

