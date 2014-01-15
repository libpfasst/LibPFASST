"""Generate several speedup/efficient plots, timing plots, and a
speedup table for the 'A parallel full approximation scheme in space
and time' paper by Emmett and Minion."""

from itertools import product
from collections import namedtuple, defaultdict

Trial = namedtuple('Trial', [ 'problem', 'processors', 'levels', 'trial', 'timings' ])

#problems   = [ 'heat' , 'burgers', 'ks' ] #'wave', 'ks' ]
# problems   = [ 'navier-stokes' ]
problems   = [ 'ks' ]
processors = [ 4, 8, 16, 32  ]
#trials     = [ 1, 2 ]
#trials     = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
trials     = [ 1, 2, 3, 4, 5 ]
#levels     = [ 2, 3 ]
levels     = [ 3 ]
#host       = 'edison.r1'
host       = 'edison'

class TimingsContainer:
    def __init__(self):
        self._parallel = defaultdict(lambda: [])
        self._serial   = {}
    def add_parallel_trial(self, trial):
        self._parallel[trial.problem].append(trial)
    def add_serial_trial(self, trial):
        self._serial[trial.problem] = trial
    def iterparallel():
        #prob, nlev, nproc, trial
        for prob in self._parallel.iterkeys():
            for nlev in set([ x.levels for x in self._parallel[prob] ]):
                for nproc in set([ x.nproc for x in self._parallel[prob] if x.levels == nlev ]):
                    pass
                    # for trial in [ x.nproc for x in self._parallel[prob] if x.levels == nlev ]:
    def parallel(prob,

class SpeedupsContainer:
    def __init__(self):
        pass


def load_timings(host='edison'):
    import pf
    import glob
    import jobtools.progressbar as pb
    import re
    import pandas

    ptimings = glob.glob('timings/%s/*_t*' % host)
    stimings = glob.glob('timings/%s/*_p*' % host)

    pparse = re.compile('.*/([-a-z]+)_t(\d+)p(\d+)l(\d+)')
    # sparse = re.compile('([-a-z]+)_t(\d+)p(\d+)l(\d+)')

    pbar  = pb.ProgressBar(widgets=["Loading timings: ", pb.Percentage(), ' ', pb.Bar()],
                           maxval=len(ptimings) + len(stimings))
    pbar.start()

    timings = TimingsContainer()

    for t in ptimings:
        problem, trial, processors, levels = pparse.match(t).group(1, 2, 3, 4)
        timings.add_parallel_trial(
            Trial(problem, int(processors), int(levels), int(trial), pf.io.read_all_timings(t)))
        pbar.bump()

    # for prob in problems:
    #     serial = 'timings/%s/%s_p01l1' % (host, prob)
    #     timings[prob, 1, 1, 0] = pf.io.read_all_timings(serial)
    #     pbar.bump()

    pbar.finish()

    return timings


# def load_timings(host='edison'):
#     import pf

#     timings = {}

#     import jobtools.progressbar as pb
#     pbar  = pb.ProgressBar(widgets=["Loading timings: ", pb.Percentage(), ' ', pb.Bar()],
#                            maxval=len(problems)*len(levels)*len(processors)*len(trials)+len(problems))
#     pbar.start()

#     for prob, nlev, nproc, trial in product(problems, levels, processors, trials):
#         parallel = 'timings/%s/%s_t%02dp%02dl%d' % (host, prob, trial, nproc, nlev)
#         timings[prob, nlev, nproc, trial] = pf.io.read_all_timings(parallel)
#         pbar.bump()

#     for prob in problems:
#         serial = 'timings/%s/%s_p01l1' % (host, prob)
#         timings[prob, 1, 1, 0] = pf.io.read_all_timings(serial)
#         pbar.bump()

#     pbar.finish()

#     return timings


def compute_speedups(timings):
    import pf
    from pf.speedup import theoretical_speedup_fixed_block as theory
    from fabfile import nnodes, nvars

    speedups = {}
    levels = [ 3 ]
    for prob, nlev, nproc, trial in timings.iterparallel():
        serial   = timings.serial(prob)
        parallel = timings.parallel(prob, nlev, nproc, trial)
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
    import jinja2

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

    speedup_table = jinja2.Template(r"""
\begin{table}
  \centering
  \begin{tabular}{llccrrr} \toprule
    Equation & Method & Processors & Iterations & Time & Speedup & Efficiency \\ \midrule
    {% for p in problems %}
      {{p.name}} & Serial SDC & 1 & {{p.serial_niters}} & {{"%.02fs"|format(p.serial_time)}} & & \\
      {% for r in p.parallel_runs %}
                 & PFASST {{r.nlevs}} & {{r.nprocs}} & {{r.niters}}
                   & {{"%.02fs \pm %.02fs"|format(*r.time)}}
                   & {{"%.02f  \pm %.02f "|format(*r.speedup)}}
                   & {{"%.02f  \pm %.02f "|format(*r.efficiency)}}
      {% endfor %}
    {% endfor %}
    \bottomrule
  \end{tabular}
  \caption{Speedups}
  \label{tab:speedup}
\end{table}
""")

    summary = lambda a: (np.mean(a), np.std(a))

    table = []

    levels = [ 3 ]


    for problem in speedups.problems:

        # HEAT & Serial SDC & 1  & 6 & 0.46s & & \\
        #      & PFASST 2   & 64 & 3 & 0.11s & 4.38 & 0.069 \\
        #      & PFASST 3   & 64 & 3 & 0.09s & 5.25 & 0.082 \\

        # stime = speedups[prob, 2, 4, 1].stime
        # print "{prob} & Serial SDC & 1 & {niters} & {stime:.02f}s & & \\\\".format(
        #     prob=prob.upper(), niters=niters[prob][1], stime=stime)

        runs = []
        for nprocs, nlevs in product(problem.processors, problem.levels):
            runs.append({ 'nlevs': nlevs, 'nprocs': nprocs, 'niters': niters[prob][nprocs],
                          'time':       summary([ x.ptime for x in s ]),
                          'speedup':    summary([ x.speedup for x in s ]),
                          'efficiency': summary([ x.efficiency for x in s]) })

        table.append({ 'name': prob.upper(),
                       'serial_niters': niters[prob][1],
                       'serial_time': speedups[prob, 3, 4, 1].stime,
                       'parallel_runs': runs })

    print speedup_table.render(problems=table)


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


def filter_trials(timings, prob, nproc, nlev):
    keys = [ k for k in timings.iterkeys() if k[:3] == (prob, nlev, nproc) ]
    return [ timings[k] for k in keys ]


def timing_boxplot(output, timings, timer, prob, nproc, nlev):
    import numpy as np
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    import rpy2
    from rpy2.robjects import r as R
    from rpy2.robjects import DataFrame

    _trials = []
    _levels = []
    _deltas = []
    for trial in trials:
        for l in range(nlev):
            d = [ x.delta for x in timings[prob, nlev, nproc, trial]
                                if x.timer == timer + str(l) ]
            _trials.extend(len(d) * [trial])
            _levels.extend(len(d) * [l])
            _deltas.extend(d)

    df = DataFrame({ 'trial': np.asarray(_trials),
                     'level': np.asarray(_levels),
                     'delta': np.asarray(_deltas) })

    rpy2.robjects.globalenv["df"] = df
    R("""library(ggplot2)
         library(ggthemes)
         plt = ggplot(df, aes(factor(level), log(delta))) +
               geom_boxplot(outlier.colour=NA) +
               geom_jitter(aes(color=factor(level))) +
               labs(x='level', y='log(elapsed time [s])') +
               facet_wrap(~trial) + theme_few()
         ggsave('{fname}', plt)""".format(fname=output))


def speedup_boxplot(output, speedups, prob, nlev):
    import numpy as np
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    import rpy2
    from rpy2.robjects import r as R
    from rpy2.robjects import DataFrame

    _nprocs = []
    _speedups = []
    _theory = []
    for nproc in processors:
        s = [ speedups[prob, nlev, nproc, trial].speedup for trial in trials ]
        t = [ speedups[prob, nlev, nproc, trial].theory for trial in trials ]
        _nprocs.extend(len(s) * [nproc])
        _speedups.extend(s)
        _theory.append(np.mean(t))

    df = DataFrame({ 'nproc': np.asarray(_nprocs),
                     'speedup': np.asarray(_speedups) })

    dft = DataFrame({ 'nproc': np.asarray(processors),
                      'theory': np.asarray(_theory) })

    rpy2.robjects.globalenv["df"] = df
    rpy2.robjects.globalenv["dft"] = dft
    R("""library(ggplot2)
         library(ggthemes)
         plt = ggplot(df, aes(factor(nproc), speedup)) +
               geom_boxplot() +
               labs(x='no. of processors') +
               geom_line(data=dft, aes(as.numeric(ordered(nproc)), theory), colour="red") +
               theme_few()
         ggsave('{fname}', plt)""".format(fname=output))


if __name__ == '__main__':
    timings = load_timings()
    speedups = compute_speedups(timings)
