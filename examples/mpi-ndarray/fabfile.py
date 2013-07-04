"""Fabric (fabfile.org) tasks for mpi-ndarray."""

import numpy as np

from fabric.api import *
from jobtools import JobQueue, Job
from itertools import product
from collections import defaultdict

nnodes = defaultdict(lambda: [ 5, 3, 2 ], { 
  'ks': [ 9, 5, 3 ]
})

nvars  = defaultdict(lambda: [ 512, 256, 128 ], { 
  'ks': [ 1024, 512, 256 ] 
})

niters = {
  'ad':      defaultdict(lambda: 8, { 1: 12 }),
  'heat':    defaultdict(lambda: 8, { 1: 12 }),
  'burgers': defaultdict(lambda: 8, { 1: 12 }),
  'ks':      defaultdict(lambda: 8, { 1: 12 }),
}

sigma = defaultdict(lambda: 0.004, { 
  'wave': 0.001 
})

dt = defaultdict(lambda: 0.01, { 
  'wave': 0.5/512,
  'ks':   1.0,
})


@task
def speed():
  """Speedup/timing tests."""

  setenv()

  jobs = JobQueue(rwd=env.scratch + 'speed', queue='regular')

  for prob, nprocs, nlevs in product( [ 'ad' ], #, 'burgers' ], #, 'ks' ],
                                      [ 4, 8, 16, 32 ],
#                                      [ 1, 4, 8, 16, 32, 64 ],
                                      [ 2, 3 ] ):

    name = '%sp%02dl%d' % (prob, nprocs, nlevs)
    job = Job(name=name, 
              param_file='probin.nml.in', 
              rwd=name, 
              width=nprocs, 
              walltime="00:10:00")

    job.update_params(
      problem=prob, rwd=name, output="", nsteps=64, dt=dt[prob], nlevs=nlevs,
      nnodes=','.join(map(str, nnodes[prob][:nlevs][::-1])), nvars=','.join(map(str, nvars[prob][:nlevs][::-1])),
      niters=niters[prob][nprocs], nu=0.005, sigma=0.004,
      )

    jobs.add(job)

  jobs.submit_all()


@task
def build(target=''):
  """Build mpi-ndarray on the remote host."""

  setenv()
  with cd(env.libpfasst):
    run('git pull')
    run('make %s' % target)


def setenv():
  """Setup Fabric and jobtools environment."""

  projects = '/home/memmett/projects/'

  if env.host[:6] == 'edison':
    env.scratch     = '/scratch/scratchdirs/memmett/'
    env.scheduler   = 'edison'
    env.host_string = 'edison.nersc.gov'
    env.host_rsync  = 'edison-s'
    env.exe         = 'main.exe'

    env.depth   = 6
    env.pernode = 4

    # XXX
    env.aprun_opts = [ '-cc numa_node' ]

  elif env.host[:5] == 'gigan':
    env.scratch     = '/scratch/memmett/'
    env.scheduler   = 'serial'
    env.host_string = 'gigan.lbl.gov'
    env.host_rsync  = 'gigan-s'
    env.exe         = 'main.exe'

    env.width = 1
    env.depth = 16

  elif env.host[:7] == 'juqueen':
    env.use_ssh_config = True
    env.scratch     = '/homea/hwu12/hwu125/scratch/'
    env.scheduler   = 'juqueen'
    env.host_string = 'juqueen'
    env.host_rsync  = 'juqueen'
    env.exe         = 'main.exe'
    env.libpfasst   = '/homea/hwu12/hwu125/projects/libpfasst/examples/mpi-ndarray'

    env.width = 1
    env.depth = 1

  else:
    env.scratch     = '/home/memmett/scratch/'
    env.scheduler   = 'serial'
    env.host_string = 'localhost'
    env.host_rsync  = 'localhost'
    env.exe         = '/home/memmett/projects/libpfasst/examples/mpi-ndarray/main.exe'

    env.width = 1
    env.depth = 2
    
  env.rsync = [ (projects + 'libpfasst', env.scratch + 'libpfasst'), ]

