
import sys
sys.path.append('/home/memmett/projects/libpfasst')

import pf
import pylab as pl
import pprint
import numpy as np

#problems = [ 'heat', 'burgers', 'wave', 'ks' ]
problems = [ 'ks' ]

pens = {
    2: { 'linestyle': 'none', 'color': 'black', 'marker': 's', 'markersize': 4 },
    3: { 'linestyle': 'none', 'color': 'black', 'marker': 'o', 'markersize': 4 },
}

for prob in problems:

    for nlevs in [ 2, 3 ]:

        send = []
        recv = []

        for nprocs in [4, 8, 16, 32, 64]:
            parallel = 'speed.out/%s_p%02dl%d' % (prob, nprocs, nlevs)

            print '====>', prob, nprocs, nlevs
            timings = pf.io.read_all_timings(parallel)

            for lev in range(nlevs):
                sends = [ x.delta for x in timings if x.timer == 'send%d' % lev ]
                recvs = [ x.delta for x in timings if x.timer == 'recv%d' % lev ]

                print 'send: mean: %.2e, std: %.2e' % (np.mean(sends), np.std(sends))
                print 'recv: mean: %.2e, std: %.2e' % (np.mean(recvs), np.std(recvs))

                send.append(sends)
                recv.append(recvs)

        pl.figure()
        pl.boxplot(send)
        pl.title('send ' + prob + ' nlevs ' + str(nlevs))

        pl.figure()
        pl.boxplot(recv)
        pl.title('recv ' + prob + ' nlevs ' + str(nlevs))

        # xy = [ (s.nproc, s.speedup) for s in speedups ]
        # x, y = zip(*xy)
        # print x
        # print y

        # pl.loglog(x, y, label='%d levels' % nlevs, **pens[nlevs])

    # pl.legend()

pl.show()
