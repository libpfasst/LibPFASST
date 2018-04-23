"""Run several regression tests."""

import subprocess
import collections
import re
from nose import with_setup

ErrorTuple = collections.namedtuple('ErrorTuple', [ 'step', 'iter', 'level', 'error' ])
EXE = 'tests/adv_diff_fft/main.exe'
NMLFILE = 'tests/adv_diff_fft/{}.nml'
TOL = 1e-4

def test_sdc():
    nml = NMLFILE.format('sdc')
    tasks = 1
    for nodes in range(4, 8):
        for status in ['', 0, 1]:
            yield advdiff, nodes, status, tasks, nml

def test_mlsdc():
    nml = NMLFILE.format('mlsdc')
    tasks = 1
    for status in [1, 2]:
        if status == 2:
            nodes = '[2 3 5]'
        else:
            nodes = 5

        yield advdiff, nodes, status, tasks, nml

def test_pfasst():
    nml = NMLFILE.format('mlsdc')
    tasks = 4
    for status in [1, 2]:
        if status == 1:
            nodes = 5
        else:
            nodes = '[2 3 5]'

        yield advdiff, nodes, status, tasks, nml

@with_setup()
def advdiff(nodes, status, tasks, nml):
    command = 'mpirun -np {} {} {} nnodes={} imex_stat={}'.format(tasks, EXE, nml, nodes, status)

    print command
    output = subprocess.check_output(command.split())

    err = errors(output)

    print err

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
            errors.append(ErrorTuple(*[c(x) for c, x in zip(cast, m.groups())]))

    return errors
