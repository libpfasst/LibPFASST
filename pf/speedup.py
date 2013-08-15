
import pf.io
import numpy as np

from collections import namedtuple

Speedup = namedtuple('Speedup', [ 'nproc', 'stime', 'ptime', 'speedup', 'efficiency', 'theory' ])


def theoretical_speedup_fixed_block(serial, parallel, nvars, nnodes, C=-1, verbose=False):
    '''Compute the theoretical speedup...'''

    # N   -  number of time steps
    # Ks  -  number of serial sdc iterations
    # Y0  -  cost of fine sweep
    # Cs  -  serial cost
    # L   -  number of pfasst levels
    # Kp  -  number of pfasst iterations
    # P   -  number of pfasst processors
    # B   -  number of pfasst blocks
    # Cp  -  parallel cost


    # serial cost
    N  = max([ x.step for x in serial ]) + 1
    Ks = max([ x.iter for x in serial ])
    Y0 = np.mean([ x.delta for x in serial if x.timer == 'sweep0' ])
    Cs = N * Ks * Y0

    if verbose:
        print 'N:  ', N
        print 'Ks: ', Ks
        print 'Y0: ', Y0
        print 'Cs: ', Cs

    # parallel cost
    L  = len(set([ x.timer for x in parallel if x.timer.startswith('sweep') ]))
    Kp = max([ x.iter for x in parallel ])
    P  = len(set([ x.rank for x in parallel ]))
    B  = N / P

    if C < 0:
        C = np.mean([ x.delta for x in parallel if x.timer == 'send0' ])

    nnodesF, nvarsF = float(nnodes[-1]), float(nvars[-1])
    alpha = [ nv/nvarsF * nn/nnodesF  for nv, nn in zip(nnodes, nvars) ]

    if verbose:
        print 'L:     ', L
        print 'Kp:    ', Kp
        print 'P:     ', P
        print 'B:     ', B
        print 'C:     ', C
        print 'alpha: ', alpha

    # Y = {}                      # sweep times
    # O = {}                      # overhead: iterp and restrict costs
    # C = np.mean([ x.delta for x in parallel if x.timer == 'send0' ]) # communication overhead
    # for lev in range(nlevs):
    #     Y[lev] = np.mean([ x.delta for x in parallel if x.timer == 'sweep%d' % lev ])
    #     if lev > 0:
    #         O[lev] = np.mean([ x.delta for x in parallel if x.timer == 'interp%d' % lev 
    #                                                      or x.timer == 'restrict%d' % lev ])
    #     else:
    #         O[lev] = 0
    # Cp = N * 2 * Y[0] +  B * Kp * ( C + sum([ 2*Y[lev] for lev in range(nlevs-1) ]) 
    #                                   + sum([   O[lev] for lev in range(nlevs) ]) + Y[nlevs-1] )

    OI = np.mean([ x.delta for x in parallel if x.timer == 'interp%d' % (L-1) ])
    OR = np.mean([ x.delta for x in parallel if x.timer == 'restrict%d' % (L-1) ])
    O0 = OI + OR

    Y = [ Y0 * alpha[l] for l in range(L) ]
    O = [ O0 * alpha[l] for l in range(L) ]

    Cp =  N * 2 * Y[0] + B * Kp * ( C + sum([ 2*Y[lev] for lev in range(L-1) ]) 
                                      + sum([   O[lev] for lev in range(L) ]) + Y[L-1] )

    
    if verbose:
        print 'Y:     ', Y
        print 'O:     ', O
        print 'Cp:    ', Cp
        print 'E:     ', Cs/Cp


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

    print "processors:           ", nproc
    print "serial time:          ", stime
    print "parallel time:        ", max(ptime)
    print "speedup (slowest):    ", stime/max(ptime)
    print "efficiency (slowest): ", stime/max(ptime)/nproc
    print "speedup (fastest):    ", stime/min(ptime)
    print "efficiency (fastest): ", stime/min(ptime)/nproc
