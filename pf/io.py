import glob
import re
import os
from collections import namedtuple

Timing = namedtuple('Timing', [
    'timer', 'rank', 'step', 'level', 'iter', 'delta', 'start', 'end'
])

Solution = namedtuple('Solution', ['step', 'iter', 'level', 'fname'])

timers_of_interest = ['exp', 'feval', 'omega']


def read_timings(fname):

    prog = re.compile(
        "timer:(.*), rank:(.*), step:(.*), level:(.*), iter:(.*), cycle:(.*), "
        + "time .rate(.*)Hz.:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)")

    with open(fname, 'r') as f:
        output = f.read()

        timings = []
        for line in output.split('\n'):
            match = prog.search(line)
            if match:
                timer = str(match.group(1)).strip()
                rate = float(match.group(7))

                if timer not in timers_of_interest:
                    continue

                rank, step, iteration = map(int, match.group(2, 3, 5))
                try:
                    level = map(int, match.group(4))
                except ValueError:
                    pass  # currently unable to print level correctly in exe
                    level = 0

                delta, start, end = map(lambda x: float(x) / rate,
                                        match.group(8, 9, 10))

                timing = Timing(timer, rank, step + 1, level, iteration, delta,
                                start, end)
                timings.append(timing)

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
            solutions.append(
                Solution(step, iteration, level,
                         os.path.join(dname, m.group(0))))

    return solutions


def read_final(dname):
    """Read output directory *dname* and return list of final solutions.

  Note that this does not read the solutions.
  """

    avail = read_avail(dname)

    solutions = []
    for step in set([x.step for x in avail]):
        tmp = [x for x in avail if x.step == step]
        max_iter = max([x.iter for x in tmp])
        max_level = max([x.level for x in tmp])

        solutions.extend(
            [x for x in tmp if x.iter == max_iter and x.level == max_level])

    return sorted(solutions)
