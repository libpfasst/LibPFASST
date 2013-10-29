"""Generate several speedup/efficient plots, timing plots, and a
speedup table for the 'A parallel full approximation scheme in space
and time' paper by Emmett and Minion."""

from itertools import product

problems   = [ 'heat' , 'burgers', 'ks' ] #'wave', 'ks' ]
processors = [ 4, 8, 16, 32, 64 ]
#trials     = [ 1, 2 ]
trials     = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
levels     = [ 2, 3 ]
#host       = 'edison.r1'


pens = {
    2: { 'linestyle': 'none', 'color': 'black',  'marker': 's', 'markersize': 8,
         'markerfacecolor': 'white', 'markeredgewidth': 2 },
    3: { 'linestyle': 'none', 'color': 'black',  'marker': 'o', 'markersize': 8 },
    (2, 'theory'): { 'linestyle': '-', 'color': 'black' },
    (3, 'theory'): { 'linestyle': '--', 'color': 'black' },
    'ideal': { 'linestyle': '-', 'color': 'black' },
}


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

    speedups = {}
    for prob, nlev, nproc, trial in product(problems, levels, processors, trials):

        serial   = timings[prob, 1, 1, 0]
        parallel = timings[prob, nlev, nproc, trial]
        sp = pf.speedup.speedup(serial, parallel)
        # XXX
        # sp = sp._replace(theory=theoretical_speedup_fixed_block(serial, parallel,
        #                                                         nnodes[prob][:nlevs],
        #                                                         nvars[prob][:nlevs], C=0.01/64*nprocs, verbose=True))
        speedups[prob, nlev, nproc, trial] = sp
    return speedups


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


def pickle(filename, phatty):
    import cPickle as p
    with open(filename, 'w') as f:
        p.dump(phatty, f)





if __name__ == '__main__':
    timings = load_timings()
    speedups = compute_speedups(timings)



