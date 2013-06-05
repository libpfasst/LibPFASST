import os
import re
import numpy as np
from collections import namedtuple
from itertools import product

solution = namedtuple('solution', [ 'step', 'iter', 'level', 'base' ])

def read_avail(dname):
  """Read directory *dname* and return available solutions."""

  solutions = []
  for fname in os.listdir(dname):
    m = re.search('(s(\d+)i(\d+)l(\d+)).npy', fname)
    if m:
      step, iteration, level = map(int, m.groups()[1:])
      solutions.append(solution(step, iteration, level, 
                                os.path.join(dname, m.group(1))))

  return solutions


def read(base):
  return np.load(base + '.npy')




ref = read_avail('ref')
out = read_avail('out')

nsteps = max([ x.step for x in ref ])

refiter  = max([ x.iter for x in ref ])
reflevel = max([ x.level for x in ref ])

outiter  = max([ x.iter for x in out ])
outlevel = max([ x.level for x in out ])

for s, i, l in product(range(nsteps+1), range(1, outiter+1), range(2, 3)):
    r = [ x.base for x in ref if x.step == s 
          and x.iter == refiter and x.level == reflevel ]
    r = read(r[0])

    a = [ x.base for x in out if x.step == s 
          and x.iter == i and x.level == l ]
    a = read(a[0])

    print 'step: %d, iter: %d, lev: %d, err: %e' % (s, i, l, abs(r - a).max())

