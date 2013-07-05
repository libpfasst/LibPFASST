
import glob
import re
import os
from collections import namedtuple


Timing = namedtuple('Timing', ['timer', 'rank', 'block', 'step', 'iter', 'cycle', 
                               'delta', 'start', 'end'])

Solution = namedtuple('Solution', [ 'step', 'iter', 'level', 'fname' ])


def read_timings(fname):

    prog = re.compile("timer:(.*), rank:(.*), block:(.*), step:(.*), iter:(.*), cycle:(.*), "
                      + "time .rate(.*)Hz.:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)")

    with open(fname, 'r') as f:
        output = f.read()

    timings = []
    for line in output.split('\n'):
        match = prog.search(line)
        if match:
            timer = str(match.group(1)).strip()
            rate  = float(match.group(7))

            rank, block, step, iteration, cycle = map(int, match.group(2, 3, 4, 5, 6))
            delta, start, end = map(lambda x: float(x)/rate, match.group(8, 10, 11))

            timings.append(Timing(timer, rank, block, step, iteration, cycle, delta, start, end))

    return timings


def read_all_timings(dname):

    timings = []
    for fname in glob.glob(dname + "/fort.*"):
        timings.extend(read_timings(fname))
    return timings



def read_avail(dname):
  """Read output directory *dname* and return list of available
  solutions.

  Note that this does not read the solutions.
  """

  prog = re.compile('(s(\d+)i(\d+)l(\d+)).npy')

  solutions = []
  for fname in os.listdir(dname):
    m = prog.search(fname)
    if m:
      step, iteration, level = map(int, m.groups()[1:])
      solutions.append(Solution(step, iteration, level, 
                                os.path.join(dname, m.group(0))))

  return solutions
