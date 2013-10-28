import sys
sys.path.append('/home/memmett/projects/libpfasst')

import pf.io
from itertools import product

import pylab
#from fabfile import nnodes, nvars

for prob, nlevs in product([ 'heat' , 'burgers', 'wave' ], #, 'ks' ],
                           [ 2, 3 ]):       

    for nprocs in [ 4, 8, 16, 32, 64 ]:
        parallel = 'speed/%s_p%02dl%d' % (prob, nprocs, nlevs)
        parallel = pf.io.read_all_timings(parallel)

        send = []
        recv = []
        for lev in range(nlevs):
            send.append([ x.delta for x in parallel if x.timer == 'send%d' % lev ])
            recv.append([ x.delta for x in parallel if x.timer == 'recv%d' % lev ])

        pylab.figure()
        pylab.boxplot(send)
        pylab.title("send: %s (%d:%d)" % (prob, nprocs, nlevs))

        pylab.figure()
        pylab.boxplot(recv)
        pylab.title("recv: %s (%d:%d)" % (prob, nprocs, nlevs))

pylab.show()
