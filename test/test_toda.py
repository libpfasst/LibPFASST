import glob
import pytest
from os import remove, mkdir
from os.path import abspath
from errno import EEXIST
import numpy as np
from pf.pfasst import PFASST, Params

# defined relative to root of project
home = abspath('.')
exe = home+'/test/toda/main.exe'
base_dir = home+'/test/toda/output'

try:
    mkdir(base_dir)
except OSError as exc:
    if exc.errno == EEXIST:
        pass
    else:
        raise OSError, exc.message

params = Params(nodes=[3], magnus=[2], \
                nsteps=128, tfinal=1.0, iterations=15, \
                nb=False, base_dir=base_dir)

ref_particle0 = 1.6125274564234153
ref_particle8 = -1.7151098584854585
toda = PFASST(exe, params)

def cleanup():
    for fname in glob.iglob(base_dir+'/*pkl'):
        remove(fname)

def make():
    tests = []
    for nodes in [[3]]:
        for magnus in [[2], [1]]:
            if nodes[0] == 2 and magnus[0] == 2: continue
            tests.append((nodes, magnus))

    return tests

tests = make()
@pytest.mark.parametrize('nodes, magnus', tests)
def test_toda(nodes, magnus):
    params = Params(nodes=nodes, magnus=magnus,
                    tolerance=1e-12,
                    nsteps=128, tfinal=1.0, iterations=30,
                    nb=False, base_dir=base_dir)

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

    assert abs(e) < 1e-5
