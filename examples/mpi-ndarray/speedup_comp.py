import sys
sys.path.append('/home/memmett/projects/libpfasst')

import pf.io
from pf.speedup import echo_speedup, speedup, theoretical_speedup_fixed_block
from itertools import product
from cPickle import dump

speedups = {}

from fabfile import nnodes, nvars

for prob, nlevs in product([ 'heat', 'burgers', 'wave', 'ks' ],
                           [ 2, 3 ]):       

    speedups[prob, nlevs] = []

    for nprocs in [ 4, 8, 16, 32, 64 ]:
        serial   = 'speed/%s_p01l1' % prob
        parallel = 'speed/%s_p%02dl%d' % (prob, nprocs, nlevs)

        print parallel

        serial   = pf.io.read_all_timings(serial)
        parallel = pf.io.read_all_timings(parallel)

        sp = speedup(serial, parallel)
        sp = sp._replace(theory=theoretical_speedup_fixed_block(serial, parallel, 
                                                                nnodes[prob][:nlevs], 
                                                                nvars[prob][:nlevs], C=0, verbose=True))
        speedups[prob, nlevs].append(sp)

        print sp.speedup

with open('speedups.pkl', 'w') as f:
    dump(speedups, f)

    

