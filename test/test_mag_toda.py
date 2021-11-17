"""Run several regression tests."""

import subprocess   #  python native module to run processes from python
import collections  #  python native module of container classes
import re           #  python native module to compare regular expressions
import pytest


ErrorTuple = collections.namedtuple('ErrorTuple', [ 'step', 'iter', 'level', 'error' ])   # This is what is scraped from the output
ResidTuple = collections.namedtuple('ResidTuple', ['time', 'step', 'rank','iter', 'level', 'resid' ])   # This is what is scraped from the output
EXE = 'test/magpicard/main.exe'                                                        # The name of the executable to run
NMLFILE = 'test/magpicard/{}.nml'                                                      # The name of the input files
TOL = 1e-5                                                                             # The error tolerance that must be met for a successful test

"""make tests from the sdc.nml input file."""
def make_sdc():
    sdc = []
    nml = NMLFILE.format('toda')
    mpi_tasks = 1
    nsteps = 32
    nodes=3
    nsteps = 32

    sdc.append((nodes, mpi_tasks, nml, nsteps))

    return sdc

"""make serial tests from the mlsdc.nml input file."""
def make_mlsdc():
    mlsdc = []
    nml = NMLFILE.format('toda')
    mpi_tasks = 1
    nodes = '[2 3 5]'
    nsteps = 32

    mlsdc.append((nodes, mpi_tasks, nml,nsteps))

    return mlsdc

"""make pfasst tests from the mlsdc.nml input file."""
def make_pfasst():
    pfasst = []
    nml = NMLFILE.format('multi_level')
    mpi_tasks = 4
    nodes = '[2 3 5]'
    nsteps = 32

    pfasst.append((nodes, mpi_tasks, nml,nsteps))

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

"""scrape the output looking for the resid statement"""
def resids(out):
    #rx = re.compile(r"resid:\s*time:\s*(\S+)\s*step:\s*(\d+)\s*rank:\s*(\d+)\s*iter:\s*(\d+)\s*level:\s*(\d+)\s*resid:\s*(\S+)")
    rx = re.compile(r"time:\s*(\S+)\s*step:\s*(\d+)\s*rank:\s*(\d+)\s*iter:\s*(\d+)\s*level:\s*(\d+)\s*resid:\s*(\S+)")
    cast = [float, int, int, int, int,  float]

    resids = []
    for line in out.splitlines():
        m = rx.search(line)
        if m:
            try:
                resids.append(ResidTuple(*[c(x) for c, x in zip(cast, m.groups())]))
            except ValueError:
                raise ValueError
                
    return resids

"""set up the list of tests to do and call the routines to make the list of tests"""
tests = []
tests.extend(make_sdc())

"""pytest command to call the following routing with the different parameters loaded into tests """
@pytest.mark.parametrize('nodes,  mpi_tasks, nml, nsteps',tests)
def test_advdiff(nodes, mpi_tasks, nml,nsteps):
    command = 'mpirun -np {} {} {} nnodes={} nsteps={}'.format(mpi_tasks, EXE, nml, nodes, nsteps)

    print(command)
    output = subprocess.check_output(command.split())   #  This will run all the tests
    output = output.decode('ascii')

    try:
        err = resids(output)   # Try to read the output for error statistics
    except ValueError:
            assert True, \
                'Expected poor convergence behavior\nunable to parse output\n {}'.format(line)
    else:   # Find the error from the last step and last iteration 
        maxstep = max([x.step for x in err])
        maxiter = max([x.iter for x in err if x.step == maxstep])
        lasterr = max([x.resid for x in err if x.step == maxstep and x.iter == maxiter])
        assert lasterr < TOL, "error: {}, tol: {}".format(lasterr, TOL)   # This decides if the test was successful
        
