import numpy as np
from pf.pfasst import PFASST, Experiment, Params

# defined relative to root of project
exe = 'tests/toda/main.exe'

params = Params(nb=False, nodes=[3], magnus=[2], tasks=1, periodic=True, \
                nsteps=2048, tfinal=10.0, particles=11, iterations=200, \
                solutions=True)

toda = PFASST(exe, params) 
results = toda.run()[0]
final_solution = results.loc[len(results)-1, 'solution']
periodic_ref = final_solution

params = Params(nb=False, nodes=[3], magnus=[2], tasks=1, periodic=False, \
                nsteps=2048, tfinal=10.0, particles=11, iterations=200, \
                solutions=True)

toda = PFASST(exe, params) 
results = toda.run()[0]
final_solution = results.loc[len(results)-1, 'solution']
nonperiodic_ref = final_solution

def test_toda11_nonperiodic_n2_m1_mpi1():
    params = Params(nodes=[2], magnus=[1], periodic=False, tasks=1, \
                    nsteps=512, tfinal=10.0, particles=11, iterations=200, \
                    nb=False)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    assert np.amax(final_solution - nonperiodic_ref) < 1e-4

def test_toda11_nonperiodic_n3_m2_mpi1():
    params = Params(nodes=[3], magnus=[2], periodic=False, tasks=1, \
                    nsteps=512, tfinal=10.0, particles=11, iterations=200, \
                    nb=False)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    assert np.amax(final_solution - nonperiodic_ref) < 1e-4

def test_toda11_periodic_n2_m1_mpi1():
    params = Params(nodes=[2], magnus=[1], periodic=True, tasks=1, \
                    nsteps=256, tfinal=10.0, particles=11, iterations=200, \
                    nb=False)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    assert np.amax(final_solution - periodic_ref) < 1e-4

def test_toda11_periodic_n3_m2_mpi1():
    params = Params(nodes=[3], magnus=[2], tasks=1, periodic=True, \
                    nsteps=256, tfinal=10.0, particles=11, iterations=200, \
                    nb=False)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    print np.amax(final_solution - periodic_ref)
    assert np.amax(final_solution - periodic_ref) < 1e-4
