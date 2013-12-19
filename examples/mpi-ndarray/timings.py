"""Generate several speedup/efficient plots, timing plots, and a
speedup table for the 'A parallel full approximation scheme in space
and time' paper by Emmett and Minion."""

from itertools import product

#problems   = [ 'heat' , 'burgers', 'ks' ] #'wave', 'ks' ]
problems   = [ 'navier-stokes' ]
processors = [ 4, 8, 16, 32  ]
#trials     = [ 1, 2 ]
#trials     = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
trials     = [ 1, 2, 3, 4, 5 ]
levels     = [ 2, 3 ]
#host       = 'edison.r1'
host       = 'edison'


def load_timings(host='edison'):
    import pf

    timings = {}

    import jobtools.progressbar as pb
    pbar  = pb.ProgressBar(widgets=["Loading timings: ", pb.Percentage(), ' ', pb.Bar()],
                           maxval=len(problems)*len(levels)*len(processors)*len(trials)+len(problems))
    pbar.start()

    for prob, nlev, nproc, trial in product(problems, levels, processors, trials):
        parallel = 'timings/%s/%s_t%02dp%02dl%d' % (host, prob, trial, nproc, nlev)
        timings[prob, nlev, nproc, trial] = pf.io.read_all_timings(parallel)
        pbar.bump()

    for prob in problems:
        serial = 'timings/%s/%s_p01l1' % (host, prob)
        timings[prob, 1, 1, 0] = pf.io.read_all_timings(serial)
        pbar.bump()

    pbar.finish()

    return timings


def compute_speedups(timings):
    import pf
    from pf.speedup import theoretical_speedup_fixed_block as theory
    from fabfile import nnodes, nvars

    speedups = {}
    levels = [ 3 ]
    for prob, nlev, nproc, trial in product(problems, levels, processors, trials):

        serial   = timings[prob, 1, 1, 0]
        parallel = timings[prob, nlev, nproc, trial]
        try:
            sp = pf.speedup.speedup(serial, parallel)
        except ValueError:
            print "warning: unable to compute speedup for:", prob, nlev, nproc, trial
            continue
        sp = sp._replace(theory=theory(serial, parallel,
                                       nnodes[prob][:nlev],
                                       nvars[prob][:nlev], verbose=False))
        print prob, nlev, nproc, trial, sp.efficiency, sp.theory/nproc
        speedups[prob, nlev, nproc, trial] = sp
    return speedups


def echo_speedup_table(speedups):

    import numpy as np

    print r"""
    \begin{table}
      \centering
      \begin{tabular}{llccrrr} \toprule
        Equation & Method & Processors & Iterations & Time & Speedup & Efficiency \\ \midrule
    """

    from collections import defaultdict
    niters = {
      'ad':      defaultdict(lambda: 8, { 1: 12 }),
      'wave':    defaultdict(lambda: 8, { 1: 12 }),
      'heat':    defaultdict(lambda: 8, { 1: 12 }),
      'burgers': defaultdict(lambda: 8, { 1: 12 }),
      'ks':      defaultdict(lambda: 8, { 1: 12 }),
      'navier-stokes': defaultdict(lambda: 5, { 1: 5 }), # XXX
    }
    # from fabfile import niters

    summary = lambda a: (np.mean(a), np.std(a))

    for prob in problems:

        # HEAT & Serial SDC & 1  & 6 & 0.46s & & \\
        #      & PFASST 2   & 64 & 3 & 0.11s & 4.38 & 0.069 \\
        #      & PFASST 3   & 64 & 3 & 0.09s & 5.25 & 0.082 \\

        stime = speedups[prob, 2, 4, 1].stime
        print "{prob} & Serial SDC & 1 & {niters} & {stime:.02f}s & & \\\\".format(
            prob=prob.upper(), niters=niters[prob][1], stime=stime)

        for nprocs, nlevs in product(processors, levels):
            s = [ speedups[prob, nlevs, nprocs, t] for t in trials ]
            ptime = summary([ x.ptime for x in s ])
            spdup = summary([ x.speedup for x in s ])
            eff   = summary([ x.efficiency for x in s])

            print " & PFASST {nlevs} & {nprocs} & {niters} & {pmean:.02f}s \pm {pstd:.02f}s & {smean:.02f} \pm {sstd:0.2f} & {emean:.02f} \pm {estd:.02f} \\\\".format(
                nlevs=nlevs, nprocs=nprocs, niters=niters[prob][nprocs],
                pmean=ptime[0], pstd=ptime[1], smean=spdup[0], sstd=spdup[1], emean=eff[0], estd=eff[1])

    print r"""
        \bottomrule
      \end{tabular}
      \caption{Speedups}
      \label{tab:speedup}
    \end{table}
    """


def dump_speedups(speedups):
    import pandas
    df = []
    for x in speedups:
        d = { 'prob': x[0] }
        d.update(speedups[x]._asdict())
        df.append(d)
    df = pandas.DataFrame(df)
    df.to_csv('speedups.csv', index=False)


def dump_level_timings(timings, timer, prob, nproc, nlev):
    import pandas
    df = []
    for trial in trials:
        for l in range(nlev):
            t = [ x.delta for x in timings[prob, nlev, nproc, trial] if x.timer == timer + str(l) ]
            for delta in t:
                df.append({ 'prob': prob, 'nlev': nlev, 'nproc': nproc, 'trial': trial, 'timer': timer,
                            'level': l, 'delta': delta })
    df = pandas.DataFrame(df)
    df.to_csv('timings_%s_p%02dl%d_%s.csv' % (prob, nproc, nlev, timer), index=False)


if __name__ == '__main__':
    timings = load_timings()
    speedups = compute_speedups(timings)



