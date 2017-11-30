"""Run several regression tests."""

import subprocess
import collections
import re
from pf.pfasst import PFASST, Experiment, Params

ErrorTuple = collections.namedtuple('ErrorTuple', [ 'step', 'iter', 'level', 'error' ])

def run(exe):
    output = subprocess.check_output(exe.split())

    return output


def errors(out):
    rx   = re.compile(r"error:\s*step:\s*(\d+)\s*iter:\s*(\d+)\s*level:\s*(\d+)\s*error:\s*(\S+)")
    cast = [ int, int, int, float ]

    errors = []
    for line in out.splitlines():
        m = rx.search(line)
        if m:
            errors.append(ErrorTuple(*[ c(x) for c, x in zip(cast, m.groups()) ]))

    return errors


def check_last_error(exe, tol):
    out = run(exe)
    err = errors(out)

    print out
    print err
    maxstep = max([ x.step for x in err ])
    maxiter = max([ x.iter for x in err if x.step == maxstep ])
    lasterr = max([ x.error for x in err if x.step == maxstep and x.iter == maxiter ])

    print maxstep, maxiter, lasterr

    print "check_last_error:", lasterr, tol

    assert lasterr < tol


def test_adv_diff_fft_serial():
    check_last_error('tests/adv_diff_fft/main.exe tests/adv_diff_fft/pipeline.nml', 5e-8)

def test_adv_diff_fft_m2():
    check_last_error('mpirun -np 2 tests/adv_diff_fft/main.exe tests/adv_diff_fft/pipeline.nml', 5e-8)

