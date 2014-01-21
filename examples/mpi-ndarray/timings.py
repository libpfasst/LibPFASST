"""Generate several speedup/efficient plots, timing plots, and a
speedup table for the 'A parallel full approximation scheme in space
and time' paper by Emmett and Minion."""

from itertools import product
from collections import namedtuple, defaultdict

TimingTrial = namedtuple('TimingTrial',
                         [ 'problem', 'processors', 'levels', 'trial', 'timings' ])

SpeedupTrial = namedtuple('SpeedupTrial',
                          [ 'problem', 'processors', 'levels', 'trial', 'speedup' ])

class TimingsContainer:
    def __init__(self):
        self.parallel = defaultdict(lambda: [])
        self.serial   = {}
    def add_parallel_trial(self, trial):
        self.parallel[trial.problem].append(trial)
    def add_serial_trial(self, problem, timings):
        self.serial[problem] = timings
    def parallel_trials(self):
        for prob in self.parallel.iterkeys():
            for t in self.parallel[prob]:
                yield prob, t


def load_timings(host='edison'):
    import pf
    import glob
    import jobtools.progressbar as pb
    import re
    import pandas

    ptimings = glob.glob('timings/%s/*_t*' % host)
    stimings = glob.glob('timings/%s/*_p*' % host)

    pparse = re.compile('.*/([-a-z]+)_t(\d+)p(\d+)l(\d+)')
    sparse = re.compile('.*/([-a-z]+)_p01')

    pbar  = pb.ProgressBar(widgets=["Loading timings: ", pb.Percentage(), ' ', pb.Bar()],
                           maxval=len(ptimings) + len(stimings))
    pbar.start()

    timings = TimingsContainer()

    for t in ptimings:
        problem, trial, processors, levels = pparse.match(t).group(1, 2, 3, 4)
        timings.add_parallel_trial(
            TimingTrial(problem, int(processors), int(levels), int(trial), pf.io.read_all_timings(t)))
        pbar.bump()

    for t in stimings:
        problem = sparse.match(t).group(1)
        timings.add_serial_trial(problem, pf.io.read_all_timings(t))
        pbar.bump()

    pbar.finish()

    return timings


def compute_speedups(timings, verbose=False):
    import pf
    from pf.speedup import theoretical_speedup_fixed_block as theory
    from fabfile import nnodes, nvars

    speedups = defaultdict(lambda: [])

    for prob, trial in timings.parallel_trials():
        serial   = timings.serial[prob]
        parallel = trial.timings
        try:
            sp = pf.speedup.speedup(serial, parallel)
        except ValueError:
            print "warning: unable to compute speedup for:", prob, trial.levels, trial.processors, trial.trial
            continue
        sp = sp._replace(theory=theory(serial, parallel,
                                       nnodes[prob][:trial.levels],
                                       nvars[prob][:trial.levels], verbose=verbose))
        speedups[prob].append(SpeedupTrial(prob, trial.processors, trial.levels, trial.trial, sp))

    return speedups



def speedup_table(speedups):
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
                   & {{"%.02fs $\pm$ %.02fs"|format(*r.time)}}
                   & {{"%.02f  $\pm$ %.02f "|format(*r.speedup)}}
                   & {{"%.02f  $\pm$ %.02f "|format(*r.efficiency)}} \\
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
    for prob in speedups:

        processors = sorted(set([ x.processors for x in speedups[prob] ]))
        levels = sorted(set([ x.levels for x in speedups[prob] ]))

        runs = []
        for nprocs, nlevs in product(processors, levels):
            s = [ x.speedup for x in speedups[prob] if x.processors == nprocs and x.levels == nlevs ]
            runs.append({ 'nlevs': nlevs, 'nprocs': nprocs, 'niters': niters[prob][nprocs],
                          'time':       summary([ x.ptime for x in s ]),
                          'speedup':    summary([ x.speedup for x in s ]),
                          'efficiency': summary([ x.efficiency for x in s]) })

        table.append({ 'name': prob.upper(),
                       'serial_niters': niters[prob][1],
                       'serial_time': speedups[prob][0].speedup.stime,
                       'parallel_runs': runs })

    return speedup_table.render(problems=table)


def timing_boxplot(output, timings, timer, prob, nproc, nlev):
    import numpy as np
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    import rpy2
    from rpy2.robjects import r as R
    from rpy2.robjects import DataFrame

    trials = sorted(set([ x.trial for x in timings.parallel[prob] if x.levels == nlev]))

    _trials = []
    _levels = []
    _deltas = []
    for trial in trials:
        t = [ x.timings for x in timings.parallel[prob]
                              if x.processors == nproc
                              and x.levels == nlev
                              and x.trial == trial ][0]
        for l in range(nlev):
            d = [ x.delta for x in t if x.timer == timer + str(l) ]
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

    processors = sorted(set([ x.processors for x in speedups[prob] if x.levels == nlev]))

    _nprocs = []
    _speedups = []
    _theory = []
    for nproc in processors:
        spds = [ x.speedup for x in speedups[prob]
                                 if x.processors == nproc and x.levels == nlev ]
        s = [ x.speedup for x in spds ]
        t = [ x.theory for x in spds ]
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


def timing_plot(output, timings, prob, nproc, nlev, trial):

    import pylab as pl

    pl.clf()

    pens = {
        'total': { 'label': 'total', 'marker': 'o', 'linestyle': '-', 'color': 'k' },
        'iteration': { 'label': 'iteration', 'marker': 's', 'linestyle': '-', 'color': 'k' },
        'predictor': { 'label': 'predictor', 'marker': 's', 'linestyle': '-', 'color': 'r' },
        ('recv', 0): { 'label': 'recv 0', 'marker': '^', 'linestyle': '-.', 'color': 'b' },
        ('recv', 1): { 'label': 'recv 1', 'marker': '^', 'linestyle': '-.', 'color': 'g' },
        ('recv', 2): { 'label': 'recv 2', 'marker': '^', 'linestyle': '-.', 'color': 'c' },
        ('send', 0): { 'label': 'send 0', 'marker': 'x', 'linestyle': ':', 'color': 'b' },
        ('send', 1): { 'label': 'send 1', 'marker': 'x', 'linestyle': ':', 'color': 'g' },
        ('send', 2): { 'label': 'send 2', 'marker': 'x', 'linestyle': ':', 'color': 'c' },
        ('sweep', 0): { 'label': 'sweep 0', 'marker': 'o', 'linestyle': '-', 'color': 'b' },
        ('sweep', 1): { 'label': 'sweep 1', 'marker': 'o', 'linestyle': '-', 'color': 'g' },
        ('sweep', 2): { 'label': 'sweep 2', 'marker': 'o', 'linestyle': '-', 'color': 'c' },
        }

    t = [ x.timings for x in timings.parallel[prob]
                          if x.processors == nproc
                          and x.levels == nlev
                          and x.trial == trial ][0]

    xmax = max([ x.rank for x in t ])

    iterations = sorted(set([ x.iter for x in t if x.iter > 0 ]))
    levels     = range(nlev)

    unpack = lambda z: zip(*sorted(z))
    select = lambda timer, k: unpack([ (x.step, x.delta) for x in t
                                       if x.timer == timer and x.iter == k ])

    ax = pl.axes([0.125,0.725,0.75-0.125,0.9-0.725])
    pl.axes(ax)

    pl.title("PFASST timing")
    pl.ylabel("total elapsed time")
    pl.xlim([0, xmax])

    x, y = unpack([ (x.rank, x.delta) for x in t if x.timer == 'total' ])
    pl.plot(x, y, 's-b', **pens['total'])


    ax = pl.axes([0.125,0.1,0.75-0.125,0.65-0.1])
    pl.axes(ax)

    x, y = unpack([ (x.step, x.delta) for x in t if x.timer == 'predictor' ])
    pl.plot(x, y, 's-b', **pens['predictor'])

    for k in iterations:
        x, y = select('iteration', k)
        pl.plot(x, y, **pens['iteration'])

        for l in levels:
            x, y = select('sweep%d' % l, k)
            pl.plot(x, y, **pens['sweep', l])

            x, y = select('recv%d' % l, k)
            pl.plot(x, y, **pens['recv', l])

            x, y = select('send%d' % l, k)
            pl.plot(x, y, **pens['send', l])

    from collections import OrderedDict
    handles, labels = pl.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    pl.figlegend(by_label.values(), by_label.keys(), loc='center right')


    pl.ylabel("elapsed time")
    pl.xlabel("processor")
    pl.xlim([0, xmax])

    pl.savefig(output)
    # pl.show()





if __name__ == '__main__':
    timings = load_timings()
    speedups = compute_speedups(timings)
