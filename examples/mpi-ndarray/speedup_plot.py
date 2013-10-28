from cPickle import load
import sys
sys.path.append('/home/memmett/projects/libpfasst')

import pylab as pl

with open('speedups.pkl', 'r') as f:
    speedups = load(f)


problems = [ 'heat', 'burgers', 'wave', 'ks' ]

pens = {
    2: { 'linestyle': 'none', 'color': 'black',  'marker': 's', 'markersize': 8, 'markerfacecolor': 'white', 'markeredgewidth': 2 },
    3: { 'linestyle': 'none', 'color': 'black',  'marker': 'o', 'markersize': 8 },
    (2, 'theory'): { 'linestyle': '-', 'color': 'black' },
    (3, 'theory'): { 'linestyle': '--', 'color': 'black' },
    'ideal': { 'linestyle': '-', 'color': 'black' },
}

for prob in problems:

    fig, ax = pl.subplots(2, sharex=True)

    for nlevs in [ 2, 3 ]:

        xy = sorted([ (s.nproc, s.speedup) for s in speedups[prob, nlevs] ])
        x, y = zip(*xy)
        ax[0].plot(x, y, label='%d levels' % nlevs, **pens[nlevs])

        xy = sorted([ (s.nproc, s.efficiency) for s in speedups[prob, nlevs] ])
        x, y = zip(*xy)
        ax[1].plot(x, y, label='%d levels' % nlevs, **pens[nlevs])


        xy = sorted([ (s.nproc, s.theory) for s in speedups[prob, nlevs] ])
        x, y = zip(*xy)
        ax[0].plot(x, y, label='%d levels, theory' % nlevs, **pens[nlevs, 'theory'])

        xy = sorted([ (s.nproc, s.theory/s.nproc) for s in speedups[prob, nlevs] ])
        x, y = zip(*xy)
        ax[1].plot(x, y, label='%d levels, theory' % nlevs, **pens[nlevs, 'theory'])


    # pl.plot([0, 70], [0, 70], '--k', label='ideal')
    #ax[0].legend(loc='best', fontsize=8)
    ax[0].legend(loc='best')
    ax[0].set_ylabel('speedup')
    ax[0].set_ylim(0, 24)

    ax[1].set_ylabel('efficiency')
    ax[1].set_xlabel('no. of processors')
    ax[1].set_ylim(0, 0.5)

    fig.suptitle(prob)
    # fig.set_tight_layout(True)

pl.show()
