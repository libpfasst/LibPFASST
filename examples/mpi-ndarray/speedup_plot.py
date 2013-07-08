from cPickle import load
import sys
sys.path.append('/home/memmett/projects/libpfasst')

import pylab as pl

with open('speedups.pkl', 'r') as f:
    speedups = load(f)


problems = [ 'heat', 'burgers', 'wave', 'ks' ]

pens = {
    2: { 'linestyle': 'none', 'color': 'blue', 'marker': 's', 'markersize': 8 },
    3: { 'linestyle': 'none', 'color': 'black', 'marker': 'o', 'markersize': 8 },
    (2, 'theory'): { 'linestyle': '-', 'color': 'blue' },
    (3, 'theory'): { 'linestyle': '-', 'color': 'black' },
    'ideal': { 'linestyle': '-', 'color': 'black' },
}

for prob in problems:

    pl.figure()

    for nlevs in [ 2, 3 ]:

        xy = sorted([ (s.nproc, s.speedup) for s in speedups[prob, nlevs] ])
        x, y = zip(*xy)

        pl.loglog(x, y, label='%d levels' % nlevs, **pens[nlevs])

        xy = sorted([ (s.nproc, s.theory) for s in speedups[prob, nlevs] ])
        x, y = zip(*xy)

        pl.loglog(x, y, label='%d levels, theory' % nlevs, **pens[nlevs, 'theory'])

    pl.legend(loc='upper left')
    pl.title(prob)
    pl.ylabel('speedup')
    pl.xlabel('no. of processors')

pl.show()
