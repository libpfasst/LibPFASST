
import pf.io
import numpy as np

from collections import namedtuple

Speedup = namedtuple('Speedup', [ 'nproc', 'stime', 'ptime', 'speedup', 'efficiency', 'theory' ])


def theoretical_speedup_fixed_block(serial, parallel, nvars, nnodes, C=-1, verbose=False):
    '''Compute the theoretical speedup assuming that PFASST was run
    block mode with a fixed number of iterations per block.

    Note that the communication overhead is taken to be maximum
    communication time over all send/recv timings.

    '''

    # serial cost
    N  = max([ x.step for x in serial ])
    Ks = max([ x.iter for x in serial ])
    Y0 = np.mean([ x.delta for x in serial if x.timer == 'sweep0' ])
    Cs = N * Ks * Y0

    # parallel cost
    L  = len(set([ x.timer for x in parallel if x.timer.startswith('sweep') ]))
    Kp = max([ x.iter for x in parallel ])
    P  = len(set([ x.rank for x in parallel ]))
    B  = N / P

    if C < 0:
        C = np.max([ x.delta for x in parallel if x.timer.startswith('send') or x.timer.startswith('recv') ])

    nnodesF, nvarsF = nnodes[-1], nvars[-1]
    alpha = [ float(nv)/nvarsF * float(nn)/nnodesF for nv, nn in zip(nnodes, nvars) ]

    OI = np.mean([ x.delta for x in parallel if x.timer == 'interp%d' % (L-1) ])
    OR = np.mean([ x.delta for x in parallel if x.timer == 'restrict%d' % (L-1) ])
    O0 = OI + OR

    Y = [ Y0 * alpha[l] for l in range(L) ]
    O = [ O0 * alpha[l] for l in range(L) ]

    Cp =  N * 2 * Y[0] + B * Kp * ( C + sum([ 2*Y[lev] for lev in range(L-1) ])
                                      + sum([   O[lev] for lev in range(L) ]) + Y[L-1] )

    if verbose:
        measCs = max([ x.delta for x in serial if x.timer == 'total' ])
        measCp = max([ x.delta for x in parallel if x.timer == 'total' ])
        measY = []
        for l in range(L):
            measY.append(np.mean([ x.delta for x in parallel if x.timer == 'sweep%d' % l ]))
        measO = []
        for l in range(L):
            OI = np.mean([ x.delta for x in parallel if x.timer == 'interp%d' % (l) ])
            OR = np.mean([ x.delta for x in parallel if x.timer == 'restrict%d' % (l) ])
            measO.append(OI + OR)

        print('')
        print('steps:                      ', N)
        print('serial iterations:          ', Ks)
        print('mean cost of serial sweeps: ', Y0)
        print('estimated serial cost:      ', Cs)
        print('actual serial cost:         ', measCs)
        print('parallel levels:            ', L)
        print('parallel iterations:        ', Kp)
        print('parallel processors:        ', P)
        print('parallel blocks:            ', B)
        print('max cost of any/all comm:   ', C)
        print('alpha ratios:               ', alpha)
        print('estimated cost of sweeps:   ', Y)
        print('actual cost of sweeps:      ', measY)
        print('estimated cost of overhead: ', O)
        print('actual cost of overhead:    ', measO)
        print('estimated parallel cost:    ', Cp)
        print('actual parallel cost:       ', measCp)
        print('estimated speedup:          ', Cs/Cp)
        print('actual speedup:             ', measCs/measCp)

    return Cs / Cp



def speedup(serial, parallel):
    '''Compute speedup/efficiency between timings in serial directory
       'sdname' and parallel directory 'pdname'.'''

    # measured
    stime = max([ x.delta for x in serial if x.timer == 'total' ])

    ptime = {}
    processors = set([ x.rank for x in parallel ])
    for proc in processors:
        ptime[proc] = max([ x.delta for x in parallel if x.timer == 'total' and x.rank == proc ])

    ptime = list(ptime.itervalues())
    nproc = len(processors)

    return Speedup(nproc, stime, max(ptime), stime/max(ptime), stime/max(ptime)/nproc, 0)


def echo_speedup(sdname, pdname):
    '''Compute speedup/efficiency between timings in serial directory
       'sdname' and parallel directory 'pdname'.'''

    serial   = pf.io.read_all_timings(sdname)
    parallel = pf.io.read_all_timings(pdname)

    stime = max([ x.delta for x in serial if x.timer == 'total' ])

    ptime = {}
    processors = set([ x.rank for x in parallel ])
    for proc in processors:
        ptime[proc] = max([ x.delta for x in parallel if x.timer == 'total' and x.rank == proc ])

    ptime = list(ptime.itervalues())
    nproc = len(processors)

    print("processors:           ", nproc)
    print("serial time:          ", stime)
    print("parallel time:        ", max(ptime))
    print("speedup (slowest):    ", stime/max(ptime))
    print("efficiency (slowest): ", stime/max(ptime)/nproc)
    print("speedup (fastest):    ", stime/min(ptime))
    print("efficiency (fastest): ", stime/min(ptime)/nproc)
