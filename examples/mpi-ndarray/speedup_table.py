from cPickle import load

from collections import defaultdict

with open('speedups.pkl', 'r') as f:
    speedups = load(f)

# taken from fabfile.py...
niters = {
  'ad':      defaultdict(lambda: 8, { 1: 12 }),
  'wave':    defaultdict(lambda: 8, { 1: 12 }),
  'heat':    defaultdict(lambda: 8, { 1: 12 }),
  'burgers': defaultdict(lambda: 8, { 1: 12 }),
  'ks':      defaultdict(lambda: 8, { 1: 12 }),
}

problems = [ 'heat', 'burgers', 'wave', 'ks' ]

print r"""
\begin{table}
  \centering
  \begin{tabular}{llccrrr} \toprule
    Equation & Method & Processors & Iterations & Time & Speedup & Efficiency \\ \midrule
"""

for prob in problems:

    # HEAT & Serial SDC & 1  & 6 & 0.46s & & \\
    #      & PFASST 2   & 64 & 3 & 0.11s & 4.38 & 0.069 \\
    #      & PFASST 3   & 64 & 3 & 0.09s & 5.25 & 0.082 \\

    stime = max([ x.stime for x in speedups[prob, 2] ])

    print "{prob} & Serial SDC & 1 & {niters} & {stime:.02f}s & & \\\\".format(
        prob=prob.upper(), niters=niters[prob][1], stime=stime)

    for nlevs in [ 2, 3 ]:

        nprocs = max([ x.nproc for x in speedups[prob, nlevs] ])
        assert len([ x.ptime for x in speedups[prob, nlevs] if x.nproc == nprocs ]) == 1
        ptime  = [ x.ptime for x in speedups[prob, nlevs] if x.nproc == nprocs ][0]
        spdup  = [ x.speedup for x in speedups[prob, nlevs] if x.nproc == nprocs ][0]
        eff    = [ x.efficiency for x in speedups[prob, nlevs] if x.nproc == nprocs ][0]

        print " & PFASST {nlevs} & {nprocs} & {niters} & {ptime:.02f}s & {spdup:.02f} & {eff:.02f} \\\\".format(
            nlevs=nlevs, nprocs=nprocs, niters=niters[prob][nprocs], ptime=ptime, spdup=spdup, eff=eff)

print r"""
    \bottomrule
  \end{tabular}
  \caption{Speedups}
  \label{tab:speedup}
\end{table}
"""


