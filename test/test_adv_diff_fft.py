"""Run several regression tests."""

import subprocess
import collections
import re
import pytest


ErrorTuple = collections.namedtuple('ErrorTuple', [ 'step', 'iter', 'level', 'error' ])
EXE = 'test/adv_diff_fft/main.exe'
NMLFILE = 'test/adv_diff_fft/{}.nml'
TOL = 1e-4

def make_sdc():
    sdc = []
    nml = NMLFILE.format('sdc')
    tasks = 1
    for nodes in range(4, 8):
        for status in ['', 0, 1]:
            sdc.append((nodes, status, tasks, nml))

    return sdc

def make_mlsdc():
    mlsdc = []
    nml = NMLFILE.format('mlsdc')
    tasks = 1
    for status in [1, 2]:
        if status == 2:
            nodes = '[2 3 5]'
        else:
            nodes = 5

        mlsdc.append((nodes, status, tasks, nml))

    return mlsdc

def make_pfasst():
    pfasst = []
    nml = NMLFILE.format('mlsdc')
    tasks = 4
    for status in [1, 2]:
        if status == 1:
            nodes = 5
        else:
            nodes = '[2 3 5]'

        pfasst.append((nodes, status, tasks, nml))

    return pfasst

tests = []
tests.extend(make_sdc())
tests.extend(make_mlsdc())
tests.extend(make_pfasst())
@pytest.mark.parametrize('nodes, status, tasks, nml',
                          tests)
def test_advdiff(nodes, status, tasks, nml):
    command = 'mpirun -np {} {} {} nnodes={} imex_stat={}'.format(tasks, EXE, nml, nodes, status)

    print command
    output = subprocess.check_output(command.split())

    try:
        err = errors(output)
    except ValueError:
        if status == 0:
            assert True, \
                'Expected poor convergence behavior\nunable to parse output\n {}'.format(line)
        else:
            assert False, \
                'Unexpected poor convergence behavior\nunable to parse output\n {}'.format(line)
    else:
        maxstep = max([x.step for x in err])
        maxiter = max([x.iter for x in err if x.step == maxstep])
        lasterr = max([x.error for x in err if x.step == maxstep and x.iter == maxiter])

        assert lasterr < TOL, "error: {}, tol: {}".format(lasterr, TOL)

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
