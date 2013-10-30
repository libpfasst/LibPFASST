"""Fabric (fabfile.org) tasks for mpi-ndarray.

Please see README for more info.
"""

import os.path
import numpy as np

from fabric.api import *
from jobtools import JobQueue, Job
from itertools import product
from collections import defaultdict


nnodes = defaultdict(lambda: [ 2, 3, 5 ], { 'ks': [ 3, 5, 9 ] })
nvars  = defaultdict(lambda: [ 128, 256, 512 ], { 'ks': [ 256, 512, 1024 ] })
niters = {
  'ad':      defaultdict(lambda: 8, { 1: 12 }),
  'wave':    defaultdict(lambda: 8, { 1: 12 }),
  'heat':    defaultdict(lambda: 8, { 1: 12 }),
  'burgers': defaultdict(lambda: 8, { 1: 12 }),
  'ks':      defaultdict(lambda: 8, { 1: 12 }),
}
sigma = defaultdict(lambda: 0.004, { 'wave': 0.001 })
dt    = defaultdict(lambda: 0.01, { 'wave': 0.5/512, 'ks': 1.0, })


# XXX: I broke the wave eqn recently...


@task
def timings():
  """Speedup/timing tests."""

  setenv()

  jobs       = JobQueue(rwd=os.path.join(env.scratch, 'timings', env.host_nick), queue='regular')
  problems   = [ 'heat', 'burgers', 'ks' ]
  # problems   = [ 'navier-stokes' ]
  processors = [ 4, 8, 16, 32, 64 ]
  trials     = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
  levels     = [ 2, 3 ]

  # serial reference run
  for prob in problems:
    name = '%s_p%02dl%d' % (prob, 1, 1)
    job  = Job(name=name, rwd=name, width=1, walltime="00:10:00")
    if prob == 'navier-stokes':
      job.exe = os.path.join(env.libpfasst, 'examples/mpi-navier-stokes/main.exe')
    else:
      job.param_file = 'probin.nml.in'
      job.update_params(
        problem=prob, rwd=name, output="", nsteps=64, dt=dt[prob], nlevs=1,
        nnodes=','.join(map(str, nnodes[prob][-1:])), nvars=','.join(map(str, nvars[prob][-1:])),
        niters=niters[prob][1], nu=0.005, sigma=sigma[prob], abs_tol=0)
    jobs.add(job)

  # parallel runs
  for prob, trial, nprocs, nlevs in product(problems, trials, processors, levels):
    name = '%s_t%02dp%02dl%d' % (prob, trial, nprocs, nlevs)
    job  = Job(name=name, rwd=name, width=nprocs, walltime="00:10:00")
    if prob == 'navier-stokes':
      job.depth       = 10
      job.specialized = 2
      job.pernode     = 2
      job.exe         = os.path.join(env.libpfasst, 'examples/mpi-navier-stokes/main.exe')
    else:
      job.param_file = 'probin.nml.in'
      job.update_params(
        problem=prob, rwd=name, output="", nsteps=64, dt=dt[prob], nlevs=nlevs,
        nnodes=','.join(map(str, nnodes[prob][-nlevs:])), nvars=','.join(map(str, nvars[prob][-nlevs:])),
        niters=niters[prob][nprocs], nu=0.005, sigma=sigma[prob], abs_tol=0)
    jobs.add(job)

  jobs.submit_all()


@task
def pull(rwd):
  setenv()

  if rwd is None:
    print 'need to specify a directory'
    return

  local("rsync -aFvz {host}:{scratch}/{rwd} .".format(
    host=env.host_rsync, scratch=env.scratch, rwd=rwd))


@task
def build(example, target=''):
  """Build given example on the remote host."""

  setenv()
  with cd(os.path.join(env.libpfasst, 'examples', example)):
    run('git pull')
    run('make %s' % target)


def setenv():
  """Setup Fabric and jobtools environment."""

  projects = '/home/memmett/projects/'

  if env.host[:6] == 'edison':
    # local
    env.host_nick   = 'edison'
    env.scratch     = '/global/scratch2/sd/memmett/PFASST'
    env.libpfasst   = '/global/homes/m/memmett/projects/libpfasst'
    env.exe         = '/global/homes/m/memmett/projects/libpfasst/examples/mpi-ndarray/main.exe'
    # jobtools
    env.host_rsync  = 'edison-s'
    env.scheduler   = 'edison'
    env.depth       = 1
    # # XXX: not sure if this is still needed
    # env.pbs_cmds    = [
    #   'export MPICH_MAX_THREAD_SAFETY=multiple',
    #   'export MPICH_NEMESIS_ASYNC_PROGRESS=method2',
    #   'export MPICH_GNI_USE_UNASSIGNED_CPUS=enabled',
    #   ]

    # fabric
    env.host_string = 'edison.nersc.gov'

  elif env.host[:5] == 'gigan':
    env.host_nick   = 'gigan'
    env.scratch     = '/scratch/memmett/'
    env.scheduler   = 'serial'
    env.host_string = 'gigan.lbl.gov'
    env.host_rsync  = 'gigan-s'
    env.exe         = 'main.exe'
    env.width       = 1
    env.depth       = 16

  elif env.host[:7] == 'juqueen':
    env.use_ssh_config = True
    env.host_nick   = 'juqueen'
    env.scratch     = '/homea/hwu12/hwu125/scratch/'
    env.scheduler   = 'juqueen'
    env.host_string = 'juqueen'
    env.host_rsync  = 'juqueen'
    env.exe         = 'main.exe'
    env.libpfasst   = '/homea/hwu12/hwu125/projects/libpfasst'
    env.width       = 1
    env.depth       = 1

  else:
    env.host_nick   = 'localhost'
    env.scratch     = '/home/memmett/scratch/'
    env.scheduler   = 'serial'
    env.host_string = 'localhost'
    env.host_rsync  = 'localhost'
    env.exe         = '/home/memmett/projects/libpfasst/examples/mpi-ndarray/main.exe'
    env.width       = 1
    env.depth       = 2

  env.rsync = [ (projects + 'libpfasst', env.scratch + 'libpfasst'), ]
