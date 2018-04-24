import glob
import pytest
from os import remove, mkdir
from os.path import abspath
from errno import EEXIST
import numpy as np
from pf.pfasst import PFASST, Params

# defined relative to root of project
home = abspath('.')
base_dir = home+'/test/imk'
exe = base_dir+'/main.exe'
output_dir = base_dir+'/output'

try:
    mkdir(base_dir)
except OSError as exc:
    if exc.errno == EEXIST:
        pass
    else:
        raise OSError, exc.message

params = Params(nodes=[3], nterms=[2], \
                nsteps=128, tfinal=1.0, \
                inttype='imk', nb=False, base_dir=output_dir)

ref_particle0 = 1.6125274564234153
ref_particle8 = -1.7151098584854585
toda = PFASST(exe, params)

def make():
    tests = []
    for vcycle in [True, False]:
        for sdc in [True, False]:
            for nlevels in [1, 3]:
                for nodes in [[3]]:
                    for nterms in [[3]]:
                        tests.append((nlevels, vcycle, sdc, nodes, nterms))
    return tests

tests = make()
@pytest.mark.parametrize('levels, vcycle, sdc, nodes, nterms',
                         tests)
def test_toda(levels, vcycle, sdc, nodes, nterms):
    params = Params(levels=levels, nodes=nodes*levels, nterms=nterms*levels,
                    sweeps=[1]*levels, sweeps_pred=[1]*levels,
                    nsteps=128, tfinal=1.0, inttype='imk',
                    nb=False, base_dir=output_dir)

    toda = PFASST(exe, params)
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    particle0 = final_solution[0, 0]
    particle8 = final_solution[7, 7]

    e0 = particle0 - ref_particle0
    e8 = particle8 - ref_particle8

    print particle0
    print e0
    print particle8
    print e8

    print ref_particle8

    e = max(e0, e8)

    print 'e = {}'.format(e)

    assert abs(e) < 1e-6
