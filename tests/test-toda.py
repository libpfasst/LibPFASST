import subprocess
import numpy as np
from nose import with_setup
from pf.pfasst import PFASST, Experiment, Params

# defined relative to root of project
base_dir = 'tests/toda/output'
exe = 'tests/toda/main.exe'

def setup_periodic():
    params = Params(nodes=[3], magnus=[2], tasks=1, periodic=True, \
                    nsteps=2048, tfinal=10.0, particles=11, iterations=200, \
                    nb=False, base_dir=base_dir, solutions=True)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']
    periodic_ref = final_solution

def setup_nonperiodic():
    params = Params(nodes=[3], magnus=[2], tasks=1, periodic=False, \
                    nsteps=2048, tfinal=10.0, particles=11, iterations=200, \
                    nb=False, base_dir=base_dir, solutions=True)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']
    nonperiodic_ref = final_solution

def cleanup():
    subprocess.call(['rm', basedir+'/*pkl'])
    

@with_setup(setup_nonperiodic, cleanup)
def test_toda11_nonperiodic_n2_m1_mpi1():
    params = Params(nodes=[2], magnus=[1], periodic=False, tasks=1, \
                    nsteps=512, tfinal=10.0, particles=11, iterations=200, \
                    nb=False, base_dir=base_dir, solutions=True)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    assert np.amax(final_solution - nonperiodic_ref) < 1e-4

@with_setup(setup_nonperiodic, cleanup)
def test_toda11_nonperiodic_n3_m2_mpi1():
    params = Params(nodes=[3], magnus=[2], periodic=False, tasks=1, \
                    nsteps=512, tfinal=10.0, particles=11, iterations=200, \
                    nb=False, base_dir=base_dir, solutions=True)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    assert np.amax(final_solution - nonperiodic_ref) < 1e-4

@with_setup(setup_periodic, cleanup)
def test_toda11_periodic_n2_m1_mpi1():
    params = Params(nodes=[2], magnus=[1], periodic=True, tasks=1, \
                    nsteps=256, tfinal=10.0, particles=11, iterations=200, \
                    nb=False, base_dir=base_dir, solutions=True)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    assert np.amax(final_solution - periodic_ref) < 1e-4

@with_setup(setup_periodic, cleanup)
def test_toda11_periodic_n3_m2_mpi1():
    params = Params(nodes=[3], magnus=[2], tasks=1, periodic=True, \
                    nsteps=256, tfinal=10.0, particles=11, iterations=200, \
                    nb=False, base_dir=base_dir, solutions=True)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    print np.amax(final_solution - periodic_ref)
    assert np.amax(final_solution - periodic_ref) < 1e-4
