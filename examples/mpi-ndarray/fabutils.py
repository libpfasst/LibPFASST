"""Fabric (fabfile.org) utilities for launching jobs."""

import glob
import os

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import *


env.mpirun   = 'mpirun'
env.nprocs   = 1
env.nthreads = 1


def setenv():
  if env.host == 'hopper':
    env.host_string = 'hopper.nersc.gov'
    env.host_rsync  = 'hopper-s'
    # env.exe = 'main.Linux.Cray.mpi.omp.exe'
    env.exe = 'main.Linux.gfortran.mpi.omp.exe'
  elif env.host == 'gigan':
    env.exe = 'main.Linux.gfortran.mpi.omp.exe'

@task
def put():
  local('/home/memmett/bin/put')


@task
def rsync():
  """Push (rsync) directories in env.rsync to env.host."""

  if env.host == 'localhost':
    return

  for src, dst in env.rsync:
      command = "rsync -avz -F {src}/ {host}:{dst}".format(
          host=env.host_rsync, src=src, dst=dst)
      local(command)



@task
def make(target=''):
  """Run *make* in the env.bin directory."""

  setenv()

  with cd(env.bin):
    run('make %s' % target)


class stage(object):
    """Context hander: create a staging directory upon entry and push
    to env.host:env.bin upon exit.
    """

    def __init__(self, stage='staging'):
        self.stage = stage

    def __enter__(self):
        self.mkdir()
        return self

    def __exit__(self, type, value, tb):
        import shutil
        puts(green('syncing stage'))
        local('rsync -auz {src}/* {host}:{dst}'.format(src=self.stage, host=env.host_rsync, dst=env.bin))
        shutil.rmtree(self.stage)

    def mkdir(self, *dirnames):
        path = os.path.join(self.stage, *dirnames)
        os.makedirs(path)
        return path + os.sep
