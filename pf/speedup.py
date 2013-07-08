
import pf.io
import numpy as np

from collections import namedtuple

Speedup = namedtuple('Speedup', [ 'nproc', 'stime', 'ptime', 'speedup', 'efficiency', 'theory' ])


def speedup(sdname, pdname):
    '''Compute speedup/efficiency between timings in serial directory
       'sdname' and parallel directory 'pdname'.'''

    serial   = pf.io.read_all_timings(sdname)
    parallel = pf.io.read_all_timings(pdname)

    # measured
    stime = max([ x.delta for x in serial if x.timer == 'total' ])

    ptime = {}
    processors = set([ x.rank for x in parallel ])
    for proc in processors:
        ptime[proc] = max([ x.delta for x in parallel if x.timer == 'total' and x.rank == proc ])

    ptime = list(ptime.itervalues())
    nproc = len(processors)

    # theory: serial cost
    N  = max([ x.step for x in serial ]) # number of time steps
    Ks = max([ x.iter for x in serial ]) # number of serial iterations
    Y0 = np.mean([ x.delta for x in serial if x.timer == 'sweep0' ]) # sweep time
    Cs = N * Ks * Y0

    # theory: parallel cost
    nlevs = len(set([ x.timer for x in parallel if x.timer.startswith('sweep') ])) # number of levels
    Kp    = max([ x.iter for x in parallel ]) # number of parallel iterations
    P     = len(set([ x.rank for x in parallel ])) # number of processors
    B     = N / P                                  # number of time blocks
    
    # here we assume that, during each iteration 2 sweeps are done on each level except the finest

    Y = {}                      # sweep times
    O = {}                      # overhead: iterp and restrict costs
    C = np.mean([ x.delta for x in parallel if x.timer == 'send0' ]) # communication overhead
    for lev in range(nlevs):
        Y[lev] = np.mean([ x.delta for x in parallel if x.timer == 'sweep%d' % lev ])
        if lev > 0:
            O[lev] = np.mean([ x.delta for x in parallel if x.timer == 'interp%d' % lev 
                                                         or x.timer == 'restrict%d' % lev ])
        else:
            O[lev] = 0
    Cp = N * 2 * Y[0] +  B * Kp * ( C + sum([ 2*Y[lev] for lev in range(nlevs-1) ]) 
                                      + sum([   O[lev] for lev in range(nlevs) ]) + Y[nlevs-1] )

    return Speedup(nproc, stime, max(ptime), stime/max(ptime), stime/max(ptime)/nproc, Cs / Cp)


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
