import sys
sys.path.append('/home/memmett/projects/libpfasst')

from pf.speedup import echo_speedup, speedup
from itertools import product
from cPickle import dump

speedups = {}

for prob, nlevs in product([ 'heat', 'burgers', 'wave', 'ks' ],
                           [ 2, 3 ]):       

    speedups[prob, nlevs] = []

    for nprocs in [ 4, 8, 16, 32, 64 ]:
        serial   = 'speed.out/%s_p01l1' % prob
        parallel = 'speed.out/%s_p%02dl%d' % (prob, nprocs, nlevs)
        print serial, parallel
        speedups[prob, nlevs].append(speedup(serial, parallel, True))

with open('speedups.pkl', 'w') as f:
    dump(speedups, f)

    

