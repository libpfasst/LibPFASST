import sys
sys.path.append('/home/memmett/projects/libpfasst')

from pf.speedup import echo_speedup, speedup

import pylab as pl

#problems = [ 'heat', 'burgers', 'wave', 'ks' ]
problems = [ 'ks' ]

pens = {
    2: { 'linestyle': 'none', 'color': 'black', 'marker': 's', 'markersize': 4 },
    3: { 'linestyle': 'none', 'color': 'black', 'marker': 'o', 'markersize': 4 },
}

for prob in problems:

    pl.figure()

    for nlevs in [ 2, 3 ]:
        speedups = []

        for nprocs in [4, 8, 16, 32, 64]:
            serial   = 'speed/%s_p01l1' % prob
            parallel = 'speed/%s_p%02dl%d' % (prob, nprocs, nlevs)
            speedups.append(speedup(serial, parallel))

            print '====>', prob, nprocs, nlevs
            echo_speedup(serial, parallel)
            print ''

        xy = [ (s.nproc, s.speedup) for s in speedups ]
        x, y = zip(*xy)
        print x
        print y

        pl.loglog(x, y, label='%d levels' % nlevs, **pens[nlevs])

    pl.legend()

pl.show()
