
import pf.io

def echo_speedup(sdname, pdname):
    '''Compute speedup/efficiency between timings in serial directory
       'sdname' and parallel directory 'pdname'.'''

    serial   = pf.io.read_all_timings(sdname)
    parallel = pf.io.real_all_timings(pdname)

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
