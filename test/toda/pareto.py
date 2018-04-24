#! /usr/bin/env python

import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from pf.pfasst import PFASST, Experiment, Params
plt.switch_backend('agg')

sns.set()
sns.set_style('whitegrid',
              {'legend.frameon': True})
sns.set_context('paper')
colors = sns.color_palette()

tfinal = 5.0
params = Params(
    tfinal=tfinal,
    tasks=1,
    sweeps_pred=[2],
    sweeps=[1],
    nodes=[3],
    magnus=[3],
    particles=11,
    iterations=200,
    tolerance=1e-12,
    qtype='gauss',
    periodic=True,
    solutions=False,
    timings=False,
    nersc=True)

markers = ['o', 'o', '^', 's', 'p', 'h', 'd']
labels = ['Leg-6', 'Leg-4-3', 'Leg-2']

if __name__ == '__main__':

    exe = '/global/homes/b/bkrull/apps/pfasst-nwchem/libpfasst/tests/toda/main.exe'
    print exe
    toda = PFASST(exe, params)
    exp = Experiment()

    try:
        with open('times{}.pkl'.format(tfinal)) as pkl:
            method_times = pickle.load(pkl)
        TIMES_RECOVERED = True
        print 'Timings were recovered'
    except IOError:
        method_times = []
        TIMES_RECOVERED = False
        print 'Timings were not recovered'

    fig, ax = plt.subplots(dpi=200)

    print 'Starting loops'

    for i, order in enumerate([3, 2, 1]):
        c = colors[i]
        toda.p.magnus = [order]
        label = labels[i]

        print "Starting order {} loop".format(order)
        mpi_times = []
        for j in range(6):
            times = []
            tasks = 2**j
            print "\t Starting MPI{}".format(tasks)
            toda.p.tasks = tasks
            nsteps = range(j, 12)

            r = exp.convergence_exp(toda, steps=nsteps)

            if TIMES_RECOVERED:
                times = method_times[i][j]
            else:
                for time in r.total_times:
                    times.append(sum(time.values()) / tasks)

            fc = c if tasks == 1 else 'none'
            m = markers[j]
            if tasks > 1:
                label = ''

            ax.plot(
                times,
                r.error,
                color=c,
                marker=m,
                markeredgewidth=1,
                markeredgecolor=colors[i],
                markerfacecolor=fc,
                label=label)

            mpi_times.append(times)
        method_times.append(mpi_times)

    with open('times{}.pkl'.format(tfinal), 'w') as pkl:
        pickle.dump(method_times, pkl)

    ax.grid(True)
    ax.legend()
    ax.set_xlabel('Time to solution')
    ax.set_ylabel('Error')
    ax.set_xscale('log')
    ax.set_yscale('log')

    fig.savefig('pareto.png')
